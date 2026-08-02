#!/usr/bin/env python
"""Missing-aware context attention MIL with CV and final-model export.

This model keeps a full fixed-tissue context schema, but each sample carries its
own observed-marker mask and context-present mask.  At inference time, a sample
can therefore provide one or more available time points while absent contexts
are excluded from attention rather than treated as zero signal.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset

from run_context_attention_mil_cv import (
    LABEL_ORDER,
    MARKER_ORDER,
    fit_log_cpm_standardizer,
    load_count_matrix,
    make_class_weights,
    metric_row,
    scheduled_learning_rate,
    set_optimizer_learning_rate,
    set_seed,
    split_train_val,
    write_metadata,
)
from sklearn.model_selection import StratifiedKFold


class MissingAwareContextBagDataset(Dataset):
    def __init__(
        self,
        counts: np.ndarray,
        labels: np.ndarray,
        indices: Sequence[int],
        library_sizes: np.ndarray,
        mean: np.ndarray,
        std: np.ndarray,
        flat_context_marker_ids: np.ndarray,
        context_marker_column_counts: np.ndarray,
        marker_mask: np.ndarray,
        n_contexts: int,
        n_markers: int,
        allowed_context_mask: Optional[np.ndarray] = None,
        augment_contexts: bool = False,
        context_dropout_prob: float = 0.0,
        single_context_prob: float = 0.0,
    ) -> None:
        self.counts = counts
        self.labels = labels
        self.indices = np.asarray(indices, dtype=np.int64)
        self.scale = (1_000_000.0 / library_sizes.astype(np.float32)).astype(np.float32)
        self.mean = mean.astype(np.float32)
        self.std = std.astype(np.float32)
        self.flat_context_marker_ids = flat_context_marker_ids.astype(np.int64)
        self.context_marker_column_counts = context_marker_column_counts.astype(np.float32)
        self.marker_mask = marker_mask.astype(np.float32)
        self.n_contexts = n_contexts
        self.n_markers = n_markers
        self.n_context_marker = n_contexts * n_markers
        self.flat_denominator = self.context_marker_column_counts.reshape(-1)
        self.present_flat = self.flat_denominator > 0
        self.augment_contexts = augment_contexts
        self.context_dropout_prob = context_dropout_prob
        self.single_context_prob = single_context_prob

        if allowed_context_mask is None:
            allowed_context_mask = self.marker_mask.sum(axis=1) > 0
        self.allowed_context_mask = np.asarray(allowed_context_mask, dtype=bool)
        if not self.allowed_context_mask.any():
            raise ValueError("allowed_context_mask contains no available context")

    def __len__(self) -> int:
        return len(self.indices)

    def _sample_context_mask(self) -> np.ndarray:
        keep = self.allowed_context_mask.copy()
        present_ids = np.flatnonzero(keep)
        if not self.augment_contexts or len(present_ids) <= 1:
            return keep

        if self.single_context_prob > 0.0 and np.random.random() < self.single_context_prob:
            chosen = int(np.random.choice(present_ids))
            keep[:] = False
            keep[chosen] = True
            return keep

        if self.context_dropout_prob > 0.0:
            random_keep = np.random.random(len(present_ids)) >= self.context_dropout_prob
            if not random_keep.any():
                random_keep[np.random.randint(0, len(present_ids))] = True
            keep[:] = False
            keep[present_ids[random_keep]] = True
        return keep

    def __getitem__(self, item: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        row_idx = int(self.indices[item])
        x = np.log2(self.counts[row_idx] * self.scale + 1.0).astype(np.float32, copy=False)
        x = ((x - self.mean) / self.std).astype(np.float32, copy=False)
        signal_flat = np.bincount(
            self.flat_context_marker_ids,
            weights=x,
            minlength=self.n_context_marker,
        ).astype(np.float32, copy=False)
        signal_flat[self.present_flat] /= self.flat_denominator[self.present_flat]
        signal = signal_flat.reshape(self.n_contexts, self.n_markers)

        keep_context = self._sample_context_mask()
        observed_marker_mask = self.marker_mask * keep_context[:, None].astype(np.float32)
        signal = signal * observed_marker_mask
        context_present = observed_marker_mask.sum(axis=1) > 0
        y = int(self.labels[row_idx])
        return (
            torch.from_numpy(signal.astype(np.float32, copy=False)),
            torch.from_numpy(observed_marker_mask.astype(np.float32, copy=False)),
            torch.from_numpy(context_present.astype(np.bool_)),
            torch.tensor(y, dtype=torch.long),
            torch.tensor(row_idx),
        )


class MissingAwareContextGatedAttentionMIL(nn.Module):
    def __init__(
        self,
        n_contexts: int,
        n_markers: int,
        n_classes: int,
        context_embedding_dim: int,
        hidden_dim: int,
        attention_dim: int,
        dropout: float,
        instance_dropout: float,
    ) -> None:
        super().__init__()
        self.n_contexts = n_contexts
        self.n_markers = n_markers
        self.instance_dropout = instance_dropout
        self.context_embedding = nn.Embedding(n_contexts, context_embedding_dim)
        self.instance_encoder = nn.Sequential(
            nn.Linear(n_markers * 2 + context_embedding_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        self.attention_v = nn.Linear(hidden_dim, attention_dim)
        self.attention_u = nn.Linear(hidden_dim, attention_dim)
        self.attention_w = nn.Linear(attention_dim, 1)
        self.classifier = nn.Linear(hidden_dim, n_classes)
        self.register_buffer("context_ids", torch.arange(n_contexts, dtype=torch.long))

    def forward(
        self,
        signal: torch.Tensor,
        observed_marker_mask: torch.Tensor,
        context_present: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        batch_size = signal.shape[0]
        context_emb = self.context_embedding(self.context_ids)
        context_emb = context_emb.unsqueeze(0).expand(batch_size, -1, -1)
        instance_input = torch.cat([signal, observed_marker_mask, context_emb], dim=-1)
        h = self.instance_encoder(instance_input)

        gated = torch.tanh(self.attention_v(h)) * torch.sigmoid(self.attention_u(h))
        attention_logits = self.attention_w(gated).squeeze(-1)
        context_present = context_present.bool()
        attention_logits = attention_logits.masked_fill(
            ~context_present,
            torch.finfo(attention_logits.dtype).min,
        )
        if self.training and self.instance_dropout > 0.0:
            keep = (torch.rand_like(attention_logits) > self.instance_dropout) & context_present
            all_dropped = ~keep.any(dim=1, keepdim=True)
            keep = keep | (all_dropped & context_present)
            attention_logits = attention_logits.masked_fill(
                ~keep,
                torch.finfo(attention_logits.dtype).min,
            )
        attention = torch.softmax(attention_logits, dim=1)
        bag_embedding = torch.sum(h * attention.unsqueeze(-1), dim=1)
        logits = self.classifier(bag_embedding)
        return logits, attention


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run missing-aware 5-fold CV and save one final model."
    )
    parser.add_argument("--count-matrix", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--label-order", nargs="+", default=LABEL_ORDER)
    parser.add_argument("--marker-order", nargs="+", default=MARKER_ORDER)
    parser.add_argument("--include-context-regex", default=None)
    parser.add_argument("--min-observed-markers", type=int, default=3)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--val-fraction", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=20260613)
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--patience", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--eval-batch-size", type=int, default=512)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--learning-rate", type=float, default=3e-3)
    parser.add_argument(
        "--lr-scheduler",
        choices=["none", "cosine", "warmup-cosine"],
        default="warmup-cosine",
    )
    parser.add_argument("--warmup-epochs", type=int, default=5)
    parser.add_argument("--min-learning-rate", type=float, default=1e-5)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--context-embedding-dim", type=int, default=32)
    parser.add_argument("--hidden-dim", type=int, default=128)
    parser.add_argument("--attention-dim", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--instance-dropout", type=float, default=0.0)
    parser.add_argument(
        "--context-dropout-prob",
        type=float,
        default=0.25,
        help="Training-time probability of dropping each present context.",
    )
    parser.add_argument(
        "--single-context-prob",
        type=float,
        default=0.25,
        help="Training-time probability of keeping exactly one context.",
    )
    parser.add_argument("--normalizer-chunk-size", type=int, default=4096)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--amp", action="store_true")
    parser.add_argument("--write-predictions", action="store_true")
    parser.add_argument("--save-fold-models", action="store_true")
    parser.add_argument("--no-final-model", action="store_true")
    parser.add_argument(
        "--eval-mode",
        action="append",
        default=None,
        help="Optional evaluation mode as name=regex. Default is full=.* only.",
    )
    return parser.parse_args()


def parse_eval_modes(mode_args: Optional[Sequence[str]]) -> List[Tuple[str, str]]:
    if not mode_args:
        return [("full", ".*")]
    modes = []
    for item in mode_args:
        if "=" not in item:
            raise ValueError(f"--eval-mode must be name=regex, got: {item}")
        name, pattern = item.split("=", 1)
        if not name:
            raise ValueError(f"Empty eval mode name in: {item}")
        modes.append((name, pattern))
    return modes


def build_eval_context_masks(contexts: Sequence[object], eval_modes: Sequence[Tuple[str, str]]) -> Dict[str, np.ndarray]:
    masks: Dict[str, np.ndarray] = {}
    context_names = [info.context for info in contexts]
    for name, pattern in eval_modes:
        regex = re.compile(pattern)
        mask = np.asarray([bool(regex.search(context)) for context in context_names], dtype=bool)
        if not mask.any():
            raise ValueError(f"Evaluation mode {name} matched no contexts with regex: {pattern}")
        masks[name] = mask
    return masks


def make_loader(
    counts: np.ndarray,
    labels: np.ndarray,
    indices: np.ndarray,
    library_sizes: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    flat_context_marker_ids: np.ndarray,
    context_marker_column_counts: np.ndarray,
    marker_mask: np.ndarray,
    n_contexts: int,
    n_markers: int,
    batch_size: int,
    shuffle: bool,
    num_workers: int,
    device: torch.device,
    allowed_context_mask: Optional[np.ndarray],
    augment_contexts: bool,
    context_dropout_prob: float,
    single_context_prob: float,
) -> DataLoader:
    dataset = MissingAwareContextBagDataset(
        counts,
        labels,
        indices,
        library_sizes,
        mean,
        std,
        flat_context_marker_ids,
        context_marker_column_counts,
        marker_mask,
        n_contexts,
        n_markers,
        allowed_context_mask=allowed_context_mask,
        augment_contexts=augment_contexts,
        context_dropout_prob=context_dropout_prob,
        single_context_prob=single_context_prob,
    )
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        pin_memory=device.type == "cuda",
    )


def new_model(
    args: argparse.Namespace,
    n_contexts: int,
    n_markers: int,
    n_classes: int,
    device: torch.device,
) -> MissingAwareContextGatedAttentionMIL:
    return MissingAwareContextGatedAttentionMIL(
        n_contexts=n_contexts,
        n_markers=n_markers,
        n_classes=n_classes,
        context_embedding_dim=args.context_embedding_dim,
        hidden_dim=args.hidden_dim,
        attention_dim=args.attention_dim,
        dropout=args.dropout,
        instance_dropout=args.instance_dropout,
    ).to(device)


def train_one_epoch_missing_aware(
    model: nn.Module,
    loader: DataLoader,
    criterion: nn.Module,
    optimizer: torch.optim.Optimizer,
    device: torch.device,
    amp: bool,
) -> float:
    model.train()
    total_loss = 0.0
    total_n = 0
    scaler = torch.amp.GradScaler(device.type, enabled=amp)
    for signal, observed_marker_mask, context_present, y, _row_idx in loader:
        signal = signal.to(device, non_blocking=True)
        observed_marker_mask = observed_marker_mask.to(device, non_blocking=True)
        context_present = context_present.to(device, non_blocking=True)
        y = y.to(device, non_blocking=True)
        optimizer.zero_grad(set_to_none=True)
        with torch.amp.autocast(device_type=device.type, enabled=amp):
            logits, _attention = model(signal, observed_marker_mask, context_present)
            loss = criterion(logits, y)
        scaler.scale(loss).backward()
        scaler.step(optimizer)
        scaler.update()
        total_loss += float(loss.detach().cpu()) * signal.size(0)
        total_n += signal.size(0)
    return total_loss / max(total_n, 1)


@torch.no_grad()
def predict_missing_aware(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
    amp: bool,
    n_contexts: int,
    collect_attention: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
    model.eval()
    row_parts: List[np.ndarray] = []
    true_parts: List[np.ndarray] = []
    prob_parts: List[np.ndarray] = []
    attention_sum = np.zeros(n_contexts, dtype=np.float64) if collect_attention else None
    attention_n = 0

    for signal, observed_marker_mask, context_present, y, row_idx in loader:
        signal = signal.to(device, non_blocking=True)
        observed_marker_mask = observed_marker_mask.to(device, non_blocking=True)
        context_present = context_present.to(device, non_blocking=True)
        with torch.amp.autocast(device_type=device.type, enabled=amp):
            logits, attention = model(signal, observed_marker_mask, context_present)
            prob = torch.softmax(logits, dim=1)
        row_parts.append(row_idx.numpy())
        true_parts.append(y.numpy())
        prob_parts.append(prob.cpu().numpy())
        if collect_attention and attention_sum is not None:
            attention_sum += attention.sum(dim=0).cpu().numpy()
            attention_n += attention.shape[0]

    rows = np.concatenate(row_parts)
    y_true = np.concatenate(true_parts)
    probs = np.concatenate(prob_parts)
    mean_attention = None
    if collect_attention and attention_sum is not None:
        mean_attention = attention_sum / max(attention_n, 1)
    return rows, y_true, probs, mean_attention


def save_predictions(
    path: Path,
    fold: int,
    eval_mode: str,
    region_ids: Sequence[str],
    label_order: Sequence[str],
    rows: np.ndarray,
    y_true: np.ndarray,
    probs: np.ndarray,
) -> None:
    y_pred = probs.argmax(axis=1)
    data = {
        "fold": np.full(len(rows), fold, dtype=np.int16),
        "eval_mode": np.full(len(rows), eval_mode),
        "region_id": [region_ids[int(idx)] for idx in rows],
        "true_label": [label_order[int(idx)] for idx in y_true],
        "predicted_label": [label_order[int(idx)] for idx in y_pred],
        "max_probability": probs.max(axis=1),
    }
    for class_idx, label in enumerate(label_order):
        data[f"prob_{label}"] = probs[:, class_idx]
    df = pd.DataFrame(data)
    header = not path.exists()
    df.to_csv(path, sep="\t", index=False, mode="a", header=header)


def evaluate_modes(
    args: argparse.Namespace,
    model: nn.Module,
    counts: np.ndarray,
    labels: np.ndarray,
    eval_idx: np.ndarray,
    library_sizes: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    flat_context_marker_ids: np.ndarray,
    context_marker_column_counts: np.ndarray,
    marker_mask: np.ndarray,
    n_contexts: int,
    n_markers: int,
    device: torch.device,
    eval_context_masks: Dict[str, np.ndarray],
) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    outputs = {}
    for eval_mode, allowed_mask in eval_context_masks.items():
        loader = make_loader(
            counts,
            labels,
            eval_idx,
            library_sizes,
            mean,
            std,
            flat_context_marker_ids,
            context_marker_column_counts,
            marker_mask,
            n_contexts,
            n_markers,
            args.eval_batch_size,
            False,
            args.num_workers,
            device,
            allowed_mask,
            False,
            0.0,
            0.0,
        )
        rows, y_true, probs, _attn = predict_missing_aware(
            model,
            loader,
            device,
            args.amp,
            n_contexts,
            collect_attention=False,
        )
        outputs[eval_mode] = (rows, y_true, probs)
    return outputs


def run_selection_training(
    args: argparse.Namespace,
    counts: np.ndarray,
    labels: np.ndarray,
    train_idx: np.ndarray,
    val_idx: Optional[np.ndarray],
    library_sizes: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    flat_context_marker_ids: np.ndarray,
    context_marker_column_counts: np.ndarray,
    marker_mask: np.ndarray,
    n_contexts: int,
    n_markers: int,
    n_classes: int,
    device: torch.device,
    fold: Optional[int] = None,
) -> Tuple[nn.Module, List[dict], dict]:
    train_loader = make_loader(
        counts,
        labels,
        train_idx,
        library_sizes,
        mean,
        std,
        flat_context_marker_ids,
        context_marker_column_counts,
        marker_mask,
        n_contexts,
        n_markers,
        args.batch_size,
        True,
        args.num_workers,
        device,
        allowed_context_mask=None,
        augment_contexts=True,
        context_dropout_prob=args.context_dropout_prob,
        single_context_prob=args.single_context_prob,
    )
    val_loader = None
    if val_idx is not None:
        val_loader = make_loader(
            counts,
            labels,
            val_idx,
            library_sizes,
            mean,
            std,
            flat_context_marker_ids,
            context_marker_column_counts,
            marker_mask,
            n_contexts,
            n_markers,
            args.eval_batch_size,
            False,
            args.num_workers,
            device,
            allowed_context_mask=None,
            augment_contexts=False,
            context_dropout_prob=0.0,
            single_context_prob=0.0,
        )

    model = new_model(args, n_contexts, n_markers, n_classes, device)
    class_weights = make_class_weights(labels, train_idx, n_classes).to(device)
    criterion = nn.CrossEntropyLoss(weight=class_weights)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=args.learning_rate,
        weight_decay=args.weight_decay,
    )

    best_state = None
    best_score = -math.inf
    best_epoch = 0
    stale_epochs = 0
    log_rows: List[Dict[str, float]] = []
    best_val_metrics: Dict[str, float] = {}
    for epoch in range(1, args.epochs + 1):
        current_lr = scheduled_learning_rate(args, epoch)
        set_optimizer_learning_rate(optimizer, current_lr)
        train_loss = train_one_epoch_missing_aware(
            model,
            train_loader,
            criterion,
            optimizer,
            device,
            args.amp,
        )
        row: Dict[str, float] = {
            "epoch": epoch,
            "train_loss": train_loss,
            "learning_rate": current_lr,
            "n_train": len(train_idx),
            "n_val": 0 if val_idx is None else len(val_idx),
        }
        if fold is not None:
            row["fold"] = fold

        if val_loader is not None:
            _rows, y_val, val_probs, _attn = predict_missing_aware(
                model,
                val_loader,
                device,
                args.amp,
                n_contexts,
                collect_attention=False,
            )
            y_val_pred = val_probs.argmax(axis=1)
            val_metrics = metric_row(y_val, y_val_pred)
            row.update({f"val_{key}": value for key, value in val_metrics.items()})
            score = val_metrics["macro_f1"]
        else:
            score = -train_loss

        log_rows.append(row)
        tag = f"[fold {fold}]" if fold is not None else f"[{args.tissue} final]"
        print(
            f"{tag} epoch {epoch} train_loss={train_loss:.5f} "
            f"lr={current_lr:.6g} score={score:.5f}",
            flush=True,
        )

        if score > best_score:
            best_score = score
            best_epoch = epoch
            best_val_metrics = dict(val_metrics) if val_loader is not None else {}
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
            stale_epochs = 0
        else:
            stale_epochs += 1
        if val_loader is not None and stale_epochs >= args.patience:
            print(f"{tag} early stopping at epoch {epoch}", flush=True)
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    selection_summary = {
        "best_epoch": best_epoch,
        "best_selection_score": best_score,
        **{f"best_val_{key}": value for key, value in best_val_metrics.items()},
    }
    return model, log_rows, selection_summary


def train_refit_model(
    args: argparse.Namespace,
    counts: np.ndarray,
    labels: np.ndarray,
    library_sizes: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    flat_context_marker_ids: np.ndarray,
    context_marker_column_counts: np.ndarray,
    marker_mask: np.ndarray,
    n_contexts: int,
    n_markers: int,
    n_classes: int,
    best_epoch: int,
    device: torch.device,
) -> Tuple[nn.Module, List[dict]]:
    all_idx = np.arange(len(labels))
    train_loader = make_loader(
        counts,
        labels,
        all_idx,
        library_sizes,
        mean,
        std,
        flat_context_marker_ids,
        context_marker_column_counts,
        marker_mask,
        n_contexts,
        n_markers,
        args.batch_size,
        True,
        args.num_workers,
        device,
        allowed_context_mask=None,
        augment_contexts=True,
        context_dropout_prob=args.context_dropout_prob,
        single_context_prob=args.single_context_prob,
    )
    model = new_model(args, n_contexts, n_markers, n_classes, device)
    class_weights = make_class_weights(labels, all_idx, n_classes).to(device)
    criterion = nn.CrossEntropyLoss(weight=class_weights)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=args.learning_rate,
        weight_decay=args.weight_decay,
    )
    log_rows: List[Dict[str, float]] = []
    for epoch in range(1, max(best_epoch, 1) + 1):
        current_lr = scheduled_learning_rate(args, epoch)
        set_optimizer_learning_rate(optimizer, current_lr)
        train_loss = train_one_epoch_missing_aware(
            model,
            train_loader,
            criterion,
            optimizer,
            device,
            args.amp,
        )
        log_rows.append(
            {
                "phase": "refit_all",
                "epoch": epoch,
                "learning_rate": current_lr,
                "train_loss": train_loss,
                "n_train": len(all_idx),
            }
        )
        print(
            f"[{args.tissue} final] refit epoch {epoch}/{best_epoch} "
            f"train_loss={train_loss:.5f} lr={current_lr:.6g}",
            flush=True,
        )
    return model, log_rows


def save_checkpoint(
    path: Path,
    model: nn.Module,
    args: argparse.Namespace,
    selected_cols: Sequence[str],
    column_infos: Sequence[object],
    contexts: Sequence[object],
    marker_mask: np.ndarray,
    context_marker_column_counts: np.ndarray,
    flat_context_marker_ids: np.ndarray,
    library_sizes: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    n_contexts: int,
    n_markers: int,
    n_classes: int,
    best_epoch: int,
    best_selection_score: float,
    final_phase: str,
) -> None:
    checkpoint = {
        "model_type": "missing_aware_context_attention_mil",
        "model_state_dict": model.state_dict(),
        "label_order": args.label_order,
        "marker_order": args.marker_order,
        "selected_cols": list(selected_cols),
        "column_infos": [asdict(info) for info in column_infos],
        "contexts": [asdict(info) for info in contexts],
        "marker_mask": marker_mask,
        "context_marker_column_counts": context_marker_column_counts,
        "flat_context_marker_ids": flat_context_marker_ids,
        "library_sizes": library_sizes,
        "mean": mean,
        "std": std,
        "args": vars(args),
        "model_config": {
            "n_contexts": n_contexts,
            "n_markers": n_markers,
            "n_classes": n_classes,
            "context_embedding_dim": args.context_embedding_dim,
            "hidden_dim": args.hidden_dim,
            "attention_dim": args.attention_dim,
            "dropout": args.dropout,
            "instance_dropout": args.instance_dropout,
        },
        "missing_aware_config": {
            "uses_per_sample_observed_marker_mask": True,
            "uses_per_sample_context_present_mask": True,
            "context_dropout_prob": args.context_dropout_prob,
            "single_context_prob": args.single_context_prob,
        },
        "best_epoch": int(best_epoch),
        "best_selection_score": float(best_selection_score),
        "final_phase": final_phase,
    }
    torch.save(checkpoint, path)


def main() -> None:
    args = parse_args()
    set_seed(args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    device = torch.device(args.device)
    if args.amp and device.type != "cuda":
        raise ValueError("--amp requires a CUDA device")

    (
        region_ids,
        labels,
        counts,
        selected_cols,
        column_infos,
        contexts,
        marker_mask,
        context_marker_column_counts,
        flat_context_marker_ids,
    ) = load_count_matrix(
        args.count_matrix,
        args.label_order,
        args.marker_order,
        args.include_context_regex,
        args.min_observed_markers,
    )
    n_classes = len(args.label_order)
    n_contexts = len(contexts)
    n_markers = len(args.marker_order)
    eval_context_masks = build_eval_context_masks(contexts, parse_eval_modes(args.eval_mode))

    write_metadata(
        args.out_dir,
        args.label_order,
        args.marker_order,
        column_infos,
        contexts,
        marker_mask,
    )
    with open(args.out_dir / "run_config.json", "w", encoding="utf-8") as handle:
        json.dump({key: str(value) for key, value in vars(args).items()}, handle, indent=2)
    pd.DataFrame(
        {
            "label": args.label_order,
            "count": np.bincount(labels, minlength=n_classes),
        }
    ).to_csv(args.out_dir / "input_label_counts.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "eval_mode": name,
                "n_contexts": int(mask.sum()),
                "contexts": ",".join([contexts[idx].context for idx in np.flatnonzero(mask)]),
            }
            for name, mask in eval_context_masks.items()
        ]
    ).to_csv(args.out_dir / "eval_modes.tsv", sep="\t", index=False)

    splitter = StratifiedKFold(n_splits=args.n_splits, shuffle=True, random_state=args.seed)
    all_indices = np.arange(len(labels))
    cv_metric_rows: List[Dict[str, float]] = []
    cv_log_rows: List[Dict[str, float]] = []
    prediction_path = args.out_dir / "cv_predictions.tsv"
    if prediction_path.exists() and args.write_predictions:
        prediction_path.unlink()

    for fold, (train_outer, test_idx) in enumerate(splitter.split(all_indices, labels), 1):
        print(f"[fold {fold}] preparing train/val/test splits", flush=True)
        train_idx, val_idx = split_train_val(
            train_outer,
            labels,
            args.val_fraction,
            args.seed + fold,
        )
        library_sizes, mean, std = fit_log_cpm_standardizer(
            counts,
            train_idx,
            args.normalizer_chunk_size,
        )
        model, log_rows, selection_summary = run_selection_training(
            args,
            counts,
            labels,
            train_idx,
            val_idx,
            library_sizes,
            mean,
            std,
            flat_context_marker_ids,
            context_marker_column_counts,
            marker_mask,
            n_contexts,
            n_markers,
            n_classes,
            device,
            fold=fold,
        )
        cv_log_rows.extend(log_rows)

        mode_outputs = evaluate_modes(
            args,
            model,
            counts,
            labels,
            test_idx,
            library_sizes,
            mean,
            std,
            flat_context_marker_ids,
            context_marker_column_counts,
            marker_mask,
            n_contexts,
            n_markers,
            device,
            eval_context_masks,
        )
        for eval_mode, (rows, y_true, probs) in mode_outputs.items():
            y_pred = probs.argmax(axis=1)
            test_metrics = metric_row(y_true, y_pred)
            fold_row = {
                "fold": fold,
                "eval_mode": eval_mode,
                "best_epoch": int(selection_summary["best_epoch"]),
                "best_selection_score": float(selection_summary["best_selection_score"]),
                "n_train": len(train_idx),
                "n_val": 0 if val_idx is None else len(val_idx),
                "n_test": len(test_idx),
                **test_metrics,
            }
            cv_metric_rows.append(fold_row)
            print(
                f"[fold {fold} {eval_mode}] test accuracy={test_metrics['accuracy']:.5f} "
                f"balanced_accuracy={test_metrics['balanced_accuracy']:.5f} "
                f"macro_f1={test_metrics['macro_f1']:.5f}",
                flush=True,
            )
            if args.write_predictions:
                save_predictions(
                    prediction_path,
                    fold,
                    eval_mode,
                    region_ids,
                    args.label_order,
                    rows,
                    y_true,
                    probs,
                )

        if args.save_fold_models:
            save_checkpoint(
                args.out_dir / f"fold{fold}.best_model.pt",
                model,
                args,
                selected_cols,
                column_infos,
                contexts,
                marker_mask,
                context_marker_column_counts,
                flat_context_marker_ids,
                library_sizes,
                mean,
                std,
                n_contexts,
                n_markers,
                n_classes,
                int(selection_summary["best_epoch"]),
                float(selection_summary["best_selection_score"]),
                "cv_best_model",
            )

    cv_metrics = pd.DataFrame(cv_metric_rows)
    cv_metrics.to_csv(args.out_dir / "cv_metrics_by_fold.tsv", sep="\t", index=False)
    pd.DataFrame(cv_log_rows).to_csv(args.out_dir / "cv_training_log.tsv", sep="\t", index=False)

    summary_rows = []
    for eval_mode, group in cv_metrics.groupby("eval_mode", sort=False):
        for metric in ["accuracy", "balanced_accuracy", "macro_f1"]:
            values = group[metric]
            summary_rows.append(
                {
                    "eval_mode": eval_mode,
                    "metric": metric,
                    "mean": values.mean(),
                    "sd": values.std(ddof=1),
                    "n_folds": len(values),
                }
            )
    pd.DataFrame(summary_rows).to_csv(args.out_dir / "cv_metric_summary.tsv", sep="\t", index=False)

    if not args.no_final_model:
        print(f"[{args.tissue} final] selecting best epoch for final model", flush=True)
        train_idx, val_idx = split_train_val(
            all_indices,
            labels,
            args.val_fraction,
            args.seed,
        )
        selection_library_sizes, selection_mean, selection_std = fit_log_cpm_standardizer(
            counts,
            train_idx,
            args.normalizer_chunk_size,
        )
        selection_model, selection_log_rows, selection_summary = run_selection_training(
            args,
            counts,
            labels,
            train_idx,
            val_idx,
            selection_library_sizes,
            selection_mean,
            selection_std,
            flat_context_marker_ids,
            context_marker_column_counts,
            marker_mask,
            n_contexts,
            n_markers,
            n_classes,
            device,
            fold=None,
        )
        final_library_sizes, final_mean, final_std = fit_log_cpm_standardizer(
            counts,
            all_indices,
            args.normalizer_chunk_size,
        )
        set_seed(args.seed)
        final_model, refit_log_rows = train_refit_model(
            args,
            counts,
            labels,
            final_library_sizes,
            final_mean,
            final_std,
            flat_context_marker_ids,
            context_marker_column_counts,
            marker_mask,
            n_contexts,
            n_markers,
            n_classes,
            int(selection_summary["best_epoch"]),
            device,
        )
        pd.DataFrame(selection_log_rows + refit_log_rows).to_csv(
            args.out_dir / "final_training_log.tsv",
            sep="\t",
            index=False,
        )
        final_summary = {
            "tissue": args.tissue,
            "model_phase": "refit_all",
            "n_rows": len(labels),
            "n_contexts": n_contexts,
            "n_markers": n_markers,
            "n_selected_columns": len(selected_cols),
            **selection_summary,
        }
        pd.DataFrame([final_summary]).to_csv(
            args.out_dir / "final_model_selection_summary.tsv",
            sep="\t",
            index=False,
        )
        save_checkpoint(
            args.out_dir / "final_model.pt",
            final_model,
            args,
            selected_cols,
            column_infos,
            contexts,
            marker_mask,
            context_marker_column_counts,
            flat_context_marker_ids,
            final_library_sizes,
            final_mean,
            final_std,
            n_contexts,
            n_markers,
            n_classes,
            int(selection_summary["best_epoch"]),
            float(selection_summary["best_selection_score"]),
            "refit_all",
        )
        del selection_model
        print(f"Wrote final model to {args.out_dir / 'final_model.pt'}", flush=True)

    print(f"Wrote missing-aware CV and final-model outputs to {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
