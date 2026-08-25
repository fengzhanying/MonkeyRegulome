# Method summary

Macaque cCREs were annotated using an ENCODE-inspired rule set adapted to the
available macaque developmental epigenomic data.

The global accessible anchor universe was defined from the accumulated macaque
ATAC peak set. Chromosome names were normalized to include the `chr` prefix, and
non-canonical contigs and mitochondrial chromosomes were removed. Each retained
anchor was assigned a stable identifier in the form `chr_start_end`.

Transcription start sites were extracted from the Macaca mulatta GTF annotation.
For positive-strand transcripts, the TSS was defined as the transcript start;
for negative-strand transcripts, the TSS was defined as the transcript end.
Coordinates were converted to BED format and deduplicated by chromosome,
position, and strand.

For each ATAC-seq or CUT&Tag bigWig file, `bigWigAverageOverBed` was used to
extract signal over every anchor. For each experiment, `log1p(mean0)` signal was
z-score normalized across all anchors. For each chromatin mark, the maximum
z-score across experiments was used as the global support score.

All ATAC peak-derived anchors were retained as accessible candidate cCREs. This
differs from the official ENCODE registry pipeline, which applies an additional
accessibility z-score threshold, because the macaque fetal dataset has fewer
conditions and the master list already derives from called ATAC peaks.

The classification threshold for H3K4me3, H3K27ac, and CTCF was z-score > 1.64.
The cCRE classes were assigned using the following priority:

1. `PLS`: within 200 bp center-to-center of a TSS and high H3K4me3.
2. `pELS`: within 2 kb of a TSS and high H3K27ac, unless classified as `PLS`.
3. `dELS`: more than 2 kb from a TSS and high H3K27ac.
4. `CA-H3K4me3`: high H3K4me3 without enhancer classification.
5. `CA-CTCF`: high CTCF with low H3K4me3 and low H3K27ac.
6. `CA`: accessible anchor without high H3K4me3, H3K27ac, or CTCF.
