# Synthetic core demo

All coordinates and signals here are synthetic. None are biological findings.

There are 20 chr1 anchors. The GTF has two plus-strand transcripts sharing TSS
1000 and one minus-strand transcript also giving BED position 1000. Deduplication
retains one plus and one minus record. The nearest-TSS fixture uses the first
plus-strand record for all anchors, consistent with the historical first-tie rule.

Four six-column bigWigAverageOverBed-style tables supply mean0=100 at:
- ATAC: first anchor.
- H3K4me3: first and fourth anchors.
- H3K27ac: second and third anchors.
- CTCF: fifth anchor.

All other mean0 values are zero. The tables and nearest-TSS fixture are supplied
directly; no synthetic bigWig files are claimed to exist. SIMULATED_PRECOMPUTED
in the demo manifest is a provenance label, not a usable path.

For N=20 and k high anchors, the sample-standardized log signals have:
- high z = sqrt((N-1)*(N-k)/(N*k))
- low z = -sqrt((N-1)*k/(N*(N-k)))

These closed-form values independently define expected/max_zscores.tsv.
The runner compares z-scores with rtol=1e-6 and atol=1e-7 to accommodate float32
rounding, and compares labels, counts and TSS rows exactly (ignoring line endings).

Expected classes are PLS, pELS, dELS, CA-H3K4me3 and CA-CTCF for the first five
anchors respectively; the remaining 15 are CA. All 20 anchors are retained,
including those with low ATAC scores.

Run from cCRE_annotation:
~~~bash
python scripts/run_demo.py --outdir demo_output
~~~

Use a fresh output directory. See the top-level README for dependencies, timing
and the exact distinction between the portable core demo and the native pipeline.
