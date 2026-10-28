# Global macaque cCRE method and pseudocode

This implementation uses an ENCODE-inspired decision hierarchy on a supplied
ATAC peak-derived global anchor universe. All retained ATAC anchors are
candidates, irrespective of their ATAC maximum z-score. This is not the official
ENCODE cCRE registry pipeline and the class labels are candidate annotations,
not proof of regulatory function.

## Inputs and reference

The historical GTF name is Macaca_mulatta.Mmul_10.114.chr.gtf: Mmul_10 and Ensembl
release 114. The exact GTF, merged ATAC BED, selected experiment manifest and
bigWigs are needed for a raw-input rerun. They are not bundled with this module.
Assembly compatibility of the input tracks must be confirmed by the data owner.

The master-list stage does not merge peaks. It adds a missing chr prefix and
excludes names containing an underscore and chrM/chrMT, preserving the original
filter. This is not a complete canonical-chromosome whitelist. Duplicate
interval IDs and invalid coordinates are now rejected. IDs are chr_start_end.
Intervals use zero-based, half-open BED coordinates.

For each GTF transcript, TSS BED start is transcript_start-1 on the plus strand
or transcript_end-1 on the minus strand; BED end is start+1. TSS records are
deduplicated by chromosome, position and strand. When duplicates have different
genes, the first encountered label is retained. TSS extraction additionally
excludes chrEBV, as in the historical code.

## Signals and normalization

bigWigAverageOverBed supplies six columns: name, size, covered, sum, mean0, mean.
mean0 averages across the full interval, including uncovered bases as zero.
Every expected anchor must appear exactly once in each experiment.

For each experiment, x = log1p(max(mean0, 0)). Calculate mean and sample standard
deviation across all anchors, including zeros. Use ddof=1. If the standard
deviation is zero, assign z=0 to every anchor; otherwise z=(x-mean)/sd.
Cast the per-experiment vector to float32, as in the original workflow.
For each mark, take the maximum z across experiments independently.
Different marks' maxima can therefore originate from different biosamples;
a global class does not imply simultaneous activity in one sample.

Missing experiments, marks or anchors and non-finite numbers cause an error.
They are not converted to biological zeros. Only constant measured signals
receive zero z-scores. Fewer than two anchors cannot be normalized.

## TSS selection and compatibility

The native command is bedtools closest -a anchors.bed -b TSS.bed -D ref -t first,
on chromosome/start sorted inputs. It selects the nearest interval (with the
first B-file record winning ties), then the script calculates the absolute
difference between the anchor midpoint and the selected one-base TSS midpoint.
It does not search globally for the nearest TSS midpoint. An overlapping wide
anchor can consequently select a different TSS from a nearest-midpoint search.

[bedtools documentation](https://bedtools.readthedocs.io/en/latest/content/tools/closest.html)
distinguishes signed upstream distances from absent B records.

The legacy mode is the default for historical reproduction. The original script
interpreted any signed distance of -1 as a missing TSS, even when it represented
a valid upstream neighbour. Legacy mode retains that behavior and the distance
sentinel 1e12. This is a documented compatibility defect, not the recommended
interpretation of signed distances.

The corrected mode identifies a missing TSS by its absent chromosome/negative
coordinate fields. A valid signed distance of -1 retains its calculated center
distance. Set CCRE_DISTANCE_MODE=corrected or --distance-mode corrected to use it.
It does not change the interval-nearest selection algorithm or thresholds.
A corrected full-data classification has not been generated; any effects on the
published counts require a fresh run with the original inputs. Archived results
are preserved. Both modes and parameters are recorded in provenance JSON.

## Classification priority

A high chromatin mark means maximum z > 1.64 (strict). Core and proximal
boundaries are inclusive: distance <=200 bp and <=2000 bp respectively.

~~~text
for each accessible anchor:
    if core and high H3K4me3:                  PLS
    else if proximal and high H3K27ac:         pELS
    else if high H3K27ac:                      dELS
    else if high H3K4me3:                      CA-H3K4me3
    else if high CTCF:                         CA-CTCF
    else:                                     CA
~~~

ATAC scores are reported but do not gate classification. Missing-TSS anchors
have a large distance and therefore cannot be PLS/pELS; they can still be dELS.
New outputs are coordinate sorted; historical files were distance sorted.
Compare labels by stable ID, not by file row order or whole-file checksum.

## Workflow pseudocode

~~~text
validate inputs, required marks and unique experiment IDs
freeze experiment manifest in a fresh workdir
prepare sorted accessible anchors and sorted transcript TSS positions
for each selected bigWig:
    extract full-interval mean0 signal, validate all anchor IDs, publish atomically
for each of the four marks:
    for each experiment:
        compute validated log1p z-scores across the same anchor universe
        accumulate element-wise maximum
select TSS intervals; calculate distances in the explicit compatibility mode
apply ordered six-class rules
write coordinate-level annotation, coordinate/class table and counts
record input/output SHA256 hashes, code hashes, versions and parameters
~~~

The pipeline starts at merged ATAC peaks and processed signal tracks.
Raw sequencing processing, sample/tissue-specific calls and cross-species
analyses are outside this module. The optional developmental table cannot be
regenerated by this global-only workflow.
