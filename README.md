# Monkey ENCODE
This repository deposit the codes for processing data of Monkey ENCODE

## RNA-seq
For RNA-seq data, run as follows:
```bash
# rna_pipeline fq1 fq2 sample_name
rna_pipeline /lustre/home/zhangfy/MacaqueTtoT-test/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R1.fq.gz /lustre/home/zhangfy/MacaqueTtoT-test/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R2.fq.gz Brain_P0_Rep1
```
## ATAC-seq
For ATAC-seq data, run every replicate as follows:
```bash
# atac_pipeline_single fq1 fq2 replicate_name
atac_pipeline_single /lustre/home/zhangfy/MacaqueTtoT-test/Brain/atac/RM22050501-Bulk-ATAC-1/RM22050501-Bulk-ATAC-1_R1.fq.gz /lustre/home/zhangfy/MacaqueTtoT-test/Brain/atac/RM22050501-Bulk-ATAC-1/RM22050501-Bulk-ATAC-1_R2.fq.gz Brain_P0_Rep1
atac_pipeline_single /lustre/home/zhangfy/MacaqueTtoT-test/Brain/atac/RM22050501-Bulk-ATAC-2/RM22050501-Bulk-ATAC-2_R1.fq.gz /lustre/home/zhangfy/MacaqueTtoT-test/Brain/atac/RM22050501-Bulk-ATAC-2/RM22050501-Bulk-ATAC-2_R2.fq.gz Brain_P0_Rep2
```
Then call peak by merging all replicates
```bash
# atac_pipeline_single ta1 ta2 sample_name
atac_pipeline_merge ta1 ta2 Brain_P0
```

## Cut&Tag
For Cut&Tag data, run as follows:
```bash
# cutag_pipeline_align fq1 fq2 sample_name
cutag_pipeline_align ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R1.fq.gz ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R2.fq.gz Brain_P0_Rep1
```
