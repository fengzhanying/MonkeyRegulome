# Monkey ENCODE
This repository deposit the codes for processing data of Monkey ENCODE

For RNA-seq data, run as follows:

```bash
# rna_pipeline fq1 fq2 sample_name
rna_pipeline ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R1.fq.gz ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R2.fq.gz Brain_P0_Rep1
```
For ATAC-seq data, run every replicate as follows:

```bash
# atac_pipeline_single fq1 fq2 sample_name
atac_pipeline_single ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R1.fq.gz ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R2.fq.gz Brain_P0_Rep1
```
Then call peak by merging all replicates
```bash
# atac_pipeline_single ta1 ta2 sample_name
atac_pipeline_callpeak ta1 ta2 Brain_P0
```

For Cut&Tag data, run as follows:

```bash
# atac_pipeline fq1 fq2 sample_name
cutag_pipeline ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R1.fq.gz ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R2.fq.gz Brain_P0_Rep1
```
