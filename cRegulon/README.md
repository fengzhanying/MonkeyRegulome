# cRegulon

Run the cRegulon workflow in this order:

**CSI → TFS → CSIS → TFModule**

The final analysis uses **K = 12 TF modules**.

## Requirements and input preparation

- Python 3 with `numpy` and `scikit-learn` installed.
- Sample identifiers in [`CSI/SampleNames.txt`](CSI/SampleNames.txt), one per line.
- PECA regulation-score matrices and gene-expression inputs for every listed sample.

The current scripts resolve paths relative to their working directory. Prepare
the following layout at the repository root before running the workflow:

```text
MonkeyRegulome/
├── PECA2/
│   ├── Input/
│   │   └── <sample>.txt
│   └── Results/
│       └── <sample>/
│           ├── TFName.txt
│           ├── TGName.txt
│           └── TFTG_regulationScore.txt
└── cRegulon/
    ├── CSI/
    │   ├── CSI.py
    │   └── SampleNames.txt
    ├── TFS/TFS.py
    ├── CSIS/CSIS.py
    └── TFModule/cRegulon.py
```

`TFTG_regulationScore.txt` is a numeric matrix without headers: rows are TFs
and columns are target genes. `TFName.txt` and `TGName.txt` contain one identifier
per line in the corresponding matrix order. All samples must share the same
TF and target-gene ordering; `CSI.py` reads the names from the first sample.

`PECA2/Input/<sample>.txt` must be a headerless, tab-separated expression table
with gene identifiers in the first column and numeric expression values in the
second. It must include every TF retained by CSI. TFS compares each sample with
the other samples, so at least two samples are required.

These upstream PECA inputs must be supplied separately. If they are stored
elsewhere, arrange the layout above or update the `../../PECA2/` paths in
`CSI/CSI.py` and `TFS/TFS.py` before running.

From the repository root, enter this directory and select your Python environment:

```bash
cd cRegulon
PYTHON_BIN="${PYTHON_BIN:-python3}"
```

The commands below are run from `cRegulon/`. Each parenthesized command changes
directory only for that step, so the next step starts from `cRegulon/` again.

## 1. CSI

[`CSI/CSI.py`](CSI/CSI.py) normalizes the TF–target regulation scores and
calculates a TF-by-TF CSI matrix for each sample. It creates three successive
sets of results: `Clean0`, `Clean1`, and `Clean2`. The later passes filter
TFs/target genes using regulation-score variability and TFs using CSI
variability across the sample collection. Downstream steps use **Clean2**.

Create the output directories before running:

```bash
mkdir -p CSI/Clean0 CSI/Clean1 CSI/Clean2
(cd CSI && "$PYTHON_BIN" CSI.py)
```

Key outputs:

- `CSI/Clean2/TFName.txt` and `CSI/Clean2/TGName.txt`
- `CSI/Clean2/<sample>_NTRS.txt`: normalized TF–target regulation scores
- `CSI/Clean2/<sample>_CSI.txt`: TF-by-TF CSI matrix

## 2. TFS

[`TFS/TFS.py`](TFS/TFS.py) reads the retained TF identifiers and PECA expression
inputs, then calculates TF expression-specificity scores. For each TF, the
sample expression is compared with its mean expression in the other samples;
the resulting score is clipped to the range [0, 1].

```bash
(cd TFS && "$PYTHON_BIN" TFS.py)
```

Output: `TFS/TFES.txt`, with TF identifiers in the first column and sample
specificity scores in the remaining columns. Running this step regenerates
the supplied `TFES.txt` for the current inputs.

## 3. CSIS

[`CSIS/CSIS.py`](CSIS/CSIS.py) combines each sample's CSI matrix with the outer
product of its TF expression-specificity scores.

```bash
(cd CSIS && "$PYTHON_BIN" CSIS.py)
```

Output: `CSIS/<sample>_CSIS.txt`, a TF-by-TF matrix in the order specified by
`CSI/Clean2/TFName.txt`.

Implementation detail: the current expression `(CSI * TFES1)**1/3` evaluates
to `(CSI * TFES1) / 3` in Python, not a cube root. This README describes the
current implementation; no scoring formula is changed by these instructions.

## 4. TFModule — final K = 12

[`TFModule/cRegulon.py`](TFModule/cRegulon.py) jointly models the sample CSIS
matrices to derive shared TF modules and sample-specific module activities.
Run the final configuration with **12 modules**:

```bash
(cd TFModule && "$PYTHON_BIN" cRegulon.py --k 12)
```

The script creates `TFModule/TFModule12/` automatically:

| Output | Description |
| --- | --- |
| `X0.txt`, `X1.txt` | NMF initialization factors |
| `X.txt` | TF-by-module membership weights; columns `C1`–`C12` |
| `XN.txt` | Standardized TF-by-module scores; columns `C1`–`C12` |
| `A.txt` | Sample-by-module activity scores; columns `C1`–`C12` |

Use **`TFModule/TFModule12/`** for the final downstream analysis.

### Optional SLURM execution

[`CSI/CSI.sh`](CSI/CSI.sh) and [`TFModule/cRegulon.sh`](TFModule/cRegulon.sh)
are cluster submission templates. Before using them, adapt the hard-coded
Python executable, partition, and resource requests to your cluster. The
TFModule template scans K = 5–30 by default. To submit only the final K = 12
configuration after completing CSI, TFS, and CSIS, override the array range:

```bash
mkdir -p TFModule/Log
(cd TFModule && sbatch --array=12 cRegulon.sh)
```

SLURM is optional; the Python commands above run the same analysis directly.
