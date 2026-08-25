# Global macaque cCRE results

This directory stores processed global macaque cCRE annotation outputs.

Recommended release files:

- `monkey_cCREs.bed.gz`: BED5 file with `chrom`, `start`, `end`, `id`, and
  `cCRE_class`.
- `monkey_cCREs_annotated.txt.gz`: full annotation table containing genomic
  coordinates, nearest TSS/gene, mark z-scores, and cCRE class.
- `cCRE_class_stats.txt`: counts per cCRE class.
- `Monkey_cCRE_Composition_Stats.csv`: publication-friendly class labels,
  counts, and percentages.
- `monkey_cCRE_developmental_annotation.tsv`: sample, time-point, and tissue
  support summary for global cCRE anchors.

Current local note:

The available local workspace contained `Monkey_cCRE_Composition_Stats.csv` and
`monkey_cCRE_developmental_annotation.tsv`. The coordinate-level global BED and
full annotation table should be copied from the completed workflow output
directory, typically `$CCRE_WORKDIR/04.classification/`.
