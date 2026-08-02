#!/usr/bin/env python
"""Predict cCRE labels with a missing-aware full-context MIL checkpoint.

Expected input is a long-format TSV with one row per element-timepoint:

region_id    stage    ATAC    CTCF    H3K27ac    ...    H3K9me3
region_1     E130     ...
region_1     E110     ...
region_2     E130     ...

Alternatively, provide --context-column if the table already contains full
context names such as E130-Brain.  Missing timepoints are naturally absent from
the table and will be masked out of attention.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset

from run_context_attention_mil_cv import metric_row
from run_missing_aware_context_attention_mil_cv_final import (
    MissingAwareContextGatedAttentionMIL,
    predict_missing_aware,
)


class MissingAwareInferenceDataset(Dataset):
    def __init__(
        self,
        counts: np.ndarray,
        observed_columns: np.ndarray,
        labels: np.ndarray,
        library_sizes: np.ndarray,
        mean: np.ndarray,
        std: np.ndarray,
        flat_context_marker_ids: np.ndarray,
        n_contexts: int,
        n_markers: int,
    ) -> None:
        self.counts = counts.astype(np.float32)
        self.observed_columns = observed_columns.astype(bool)
        self.labels = labels.astype(np.int64)
        self.scale = (1_000_000.0 / library_sizes.astype(np.float32)).astype(np.float32)
        self.mean = mean.astype(np.float32)
        self.std = std.astype(np.float32)
        self.flat_context_marker_ids = flat_context_marker_ids.astype(np.int64)
        self.n_contexts = n_contexts
        self.n_markers = n_markers
        self.n_context_marker = n_contexts * n_markers

    def __len__(self) -> int:
        return self.counts.shape[0]

    def __getitem__(self, item: int):
        observed = self.observed_columns[item].astype(np.float32)
        x = np.log2(self.counts[item] * self.scale + 1.0).astype(np.float32, copy=False)
        x = ((x - self.mean) / self.std).astype(np.float32, copy=False)
        weights = x * observed
        signal_flat = np.bincount(
            self.flat_context_marker_ids,
            weights=weights,
            minlength=self.n_context_marker,
        ).astype(np.float32, copy=False)
        denominator = np.bincount(
            self.flat_context_marker_ids,
            weights=observed,
            minlength=self.n_context_marker,
        ).astype(np.float32, copy=False)
        present = denominator > 0
        signal_flat[present] /= denominator[present]
        signal = signal_flat.reshape(self.n_contexts, self.n_markers)
        observed_marker_mask = present.reshape(self.n_contexts, self.n_markers).astype(np.float32)
        context_present = observed_marker_mask.sum(axis=1) > 0
        return (
            torch.from_numpy(signal),
            torch.from_numpy(observed_marker_mask),
            torch.from_numpy(context_present.astype(np.bool_)),
            torch.tensor(int(self.labels[item]), dtype=torch.long),
            torch.tensor(item),
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Predict with a missing-aware MIL final_model.pt")
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--id-column", default=None)
    parser.add_argument("--stage-column", default="stage")
    parser.add_argument("--context-column", default=None)
    parser.add_argument("--label-column", default=None)
    parser.add_argument("--eval-batch-size", type=int, default=512)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--amp", action="store_true")
    return parser.parse_args()


def load_checkpoint(path: Path, device: torch.device) -> dict:
    try:
        return torch.load(path, map_location=device, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=device)


def choose_id_column(df: pd.DataFrame, requested: Optional[str]) -> str:
    if requested:
        if requested not in df.columns:
            raise ValueError(f"Requested id column not found: {requested}")
        return requested
    for candidate in ["region_id", "cCRE", "id", "ID", "name"]:
        if candidate in df.columns:
            return candidate
    return df.columns[0]


def normalize_marker_lookup(df: pd.DataFrame, marker_order: Sequence[str]) -> Dict[str, str]:
    lookup = {}
    for column in df.columns:
        for key in [column, column.strip(), column.strip().upper(), column.strip().lower()]:
            lookup.setdefault(key, column)
    for marker in marker_order:
        for column in df.columns:
            if column == marker or column.upper() == marker.upper():
                lookup.setdefault(marker, column)
                lookup.setdefault(marker.upper(), column)
                lookup.setdefault(marker.lower(), column)
    return lookup


def context_id_for_record(
    row: pd.Series,
    context_column: Optional[str],
    stage_column: str,
    context_to_id: Dict[str, int],
    stage_to_context_ids: Dict[str, List[int]],
) -> int:
    if context_column:
        context = str(row[context_column])
        if context in context_to_id:
            return context_to_id[context]
        if stage_column in row.index and not pd.isna(row[stage_column]):
            stage = str(row[stage_column])
        else:
            stage = context.split("-", 1)[0]
        matches = stage_to_context_ids.get(stage, [])
        if len(matches) == 1:
            return matches[0]
        available = ", ".join(list(context_to_id.keys())[:20])
        if len(context_to_id) > 20:
            available += ", ..."
        if not matches:
            raise ValueError(
                f"Input context is not in model schema and stage fallback also failed: "
                f"context={context}, stage={stage}. Available model contexts: {available}"
            )
        raise ValueError(
            f"Input context is not in model schema: {context}. Stage {stage} maps to "
            f"multiple model contexts; use a context value from the checkpoint schema. "
            f"Available model contexts: {available}"
        )
    stage = str(row[stage_column])
    matches = stage_to_context_ids.get(stage, [])
    if not matches:
        available = ", ".join(list(context_to_id.keys())[:20])
        if len(context_to_id) > 20:
            available += ", ..."
        raise ValueError(
            f"Input stage is not in model schema: {stage}. "
            f"Available model contexts: {available}"
        )
    if len(matches) > 1:
        raise ValueError(
            f"Stage {stage} maps to multiple model contexts; use --context-column."
        )
    return matches[0]


def build_feature_arrays(
    df: pd.DataFrame,
    id_col: str,
    label_col: Optional[str],
    context_column: Optional[str],
    stage_column: str,
    checkpoint: dict,
) -> Tuple[List[str], np.ndarray, np.ndarray, np.ndarray]:
    label_order = list(checkpoint["label_order"])
    marker_order = list(checkpoint["marker_order"])
    column_infos = list(checkpoint["column_infos"])
    contexts = list(checkpoint["contexts"])
    n_columns = len(checkpoint["selected_cols"])
    region_ids = df[id_col].astype(str).drop_duplicates().tolist()
    region_to_row = {region_id: idx for idx, region_id in enumerate(region_ids)}

    context_to_id = {info["context"]: int(info["context_id"]) for info in contexts}
    stage_to_context_ids: Dict[str, List[int]] = {}
    for info in contexts:
        stage_to_context_ids.setdefault(info["stage"], []).append(int(info["context_id"]))

    marker_to_column_indices: Dict[Tuple[int, str], List[int]] = {}
    for column_idx, info in enumerate(column_infos):
        marker_to_column_indices.setdefault(
            (int(info["context_id"]), info["marker"]),
            [],
        ).append(column_idx)

    marker_lookup = normalize_marker_lookup(df, marker_order)
    counts = np.zeros((len(region_ids), n_columns), dtype=np.float32)
    observed = np.zeros((len(region_ids), n_columns), dtype=bool)
    labels = np.zeros(len(region_ids), dtype=np.int64)
    label_to_id = {label: idx for idx, label in enumerate(label_order)}

    if label_col is not None:
        if label_col not in df.columns:
            raise ValueError(f"Requested label column not found: {label_col}")
        first_labels = df.groupby(id_col, sort=False)[label_col].first()
        unknown = sorted(set(first_labels) - set(label_order))
        if unknown:
            raise ValueError(f"Unexpected labels in {label_col}: {unknown}")
        for region_id, label in first_labels.items():
            labels[region_to_row[str(region_id)]] = label_to_id[label]

    for _, row in df.iterrows():
        region_id = str(row[id_col])
        row_idx = region_to_row[region_id]
        context_id = context_id_for_record(
            row,
            context_column,
            stage_column,
            context_to_id,
            stage_to_context_ids,
        )
        for marker in marker_order:
            source_column = marker_lookup.get(marker) or marker_lookup.get(marker.upper()) or marker_lookup.get(marker.lower())
            if source_column is None or pd.isna(row[source_column]):
                continue
            target_indices = marker_to_column_indices.get((context_id, marker), [])
            if not target_indices:
                continue
            value = float(row[source_column])
            counts[row_idx, target_indices] = value
            observed[row_idx, target_indices] = True

    empty = np.flatnonzero(observed.sum(axis=1) == 0)
    if len(empty):
        examples = ", ".join(region_ids[idx] for idx in empty[:5])
        raise ValueError(f"Some regions have no observed model features. First examples: {examples}")
    return region_ids, counts, observed, labels


def main() -> None:
    args = parse_args()
    device = torch.device(args.device)
    if args.amp and device.type != "cuda":
        raise ValueError("--amp requires a CUDA device")

    checkpoint = load_checkpoint(args.model, device)
    if checkpoint.get("model_type") != "missing_aware_context_attention_mil":
        raise ValueError("Checkpoint is not a missing-aware context attention MIL model")
    model_config = dict(checkpoint["model_config"])
    model = MissingAwareContextGatedAttentionMIL(
        n_contexts=int(model_config["n_contexts"]),
        n_markers=int(model_config["n_markers"]),
        n_classes=int(model_config["n_classes"]),
        context_embedding_dim=int(model_config["context_embedding_dim"]),
        hidden_dim=int(model_config["hidden_dim"]),
        attention_dim=int(model_config["attention_dim"]),
        dropout=float(model_config["dropout"]),
        instance_dropout=float(model_config["instance_dropout"]),
    ).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])

    df = pd.read_csv(args.features, sep="\t")
    id_col = choose_id_column(df, args.id_column)
    context_column = args.context_column
    if context_column is not None and context_column not in df.columns:
        raise ValueError(f"Requested context column not found: {context_column}")
    if context_column is None and args.stage_column not in df.columns:
        raise ValueError(f"Stage column not found: {args.stage_column}")
    label_col = args.label_column or ("Gold" if "Gold" in df.columns else None)
    region_ids, counts, observed, labels = build_feature_arrays(
        df,
        id_col,
        label_col,
        context_column,
        args.stage_column,
        checkpoint,
    )

    dataset = MissingAwareInferenceDataset(
        counts,
        observed,
        labels,
        np.asarray(checkpoint["library_sizes"], dtype=np.float32),
        np.asarray(checkpoint["mean"], dtype=np.float32),
        np.asarray(checkpoint["std"], dtype=np.float32),
        np.asarray(checkpoint["flat_context_marker_ids"], dtype=np.int64),
        int(model_config["n_contexts"]),
        int(model_config["n_markers"]),
    )
    loader = DataLoader(
        dataset,
        batch_size=args.eval_batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        pin_memory=device.type == "cuda",
    )
    rows, y_true, probs, _attn = predict_missing_aware(
        model,
        loader,
        device,
        args.amp,
        int(model_config["n_contexts"]),
        collect_attention=False,
    )
    label_order = list(checkpoint["label_order"])
    y_pred = probs.argmax(axis=1)
    out = pd.DataFrame(
        {
            "region_id": [region_ids[int(idx)] for idx in rows],
            "predicted_label": [label_order[int(idx)] for idx in y_pred],
            "max_probability": probs.max(axis=1),
        }
    )
    if label_col is not None:
        out.insert(1, "true_label", [label_order[int(idx)] for idx in y_true])
    for class_idx, label in enumerate(label_order):
        out[f"prob_{label}"] = probs[:, class_idx]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote predictions to {args.out}", flush=True)

    if label_col is not None:
        metrics_path = args.out.with_suffix(".metrics.tsv")
        pd.DataFrame([metric_row(y_true, y_pred)]).to_csv(metrics_path, sep="\t", index=False)
        print(f"Wrote metrics to {metrics_path}", flush=True)


if __name__ == "__main__":
    main()
