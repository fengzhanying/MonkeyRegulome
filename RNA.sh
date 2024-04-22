### STAR Index
../bin/Linux_x86_64_static/STAR --runThreadN 24 --runMode genomeGenerate --genomeDir ./. --genomeFastaFiles ../../../MacaqueTtoT-test/Macaque.t-to-t/Macaque_Assembly.v3-for-encode.fa --sjdbGTFfile ../../../MacaqueTtoT-test/Macaque.t-to-t/Macaque.t-to-t.sorted.final.withoutMT.gtf

./fastp -i ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R1.fq.gz -I ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R2.fq.gz -o ./Brain_P0_R1.trimmed.fq.gz -O ./Brain_P0_R2.trimmed.fq.gz --detect_adapter_for_pe --thread 24 --json Brain_P0_fastp.json --html Brain_P0_fastp.html 2> Brain_P0_fastp.log

### task align
../STAR-2.7.11b/bin/Linux_x86_64_static/STAR  --genomeDir ../STAR-2.7.11b/index/ --readFilesIn ./Brain_P0_R1.trimmed.fq.gz ./Brain_P0_R2.trimmed.fq.gz --readFilesCommand zcat --runThreadN 24 --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1

### task bam_to_signals
../STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph
../STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph --outWigReferencesPrefix chr

bedSort Signal.UniqueMultiple.str1.out.bg Signal.UniqueMultiple.str1.out.bg
bedGraphToBigWig Signal.UniqueMultiple.str1.out.bg {chrom_sizes} *_minusAll.bw

bedSort Signal.Unique.str1.out.bg Signal.Unique.str1.out.bg
bedGraphToBigWig Signal.Unique.str1.out.bg {chrom_sizes} *_minusUniq.bw

bedSort  Signal.UniqueMultiple.str2.out.bg Signal.UniqueMultiple.str2.out.bg
bedGraphToBigWig Signal.UniqueMultiple.str2.out.bg {chrom_sizes} *_plusAll.bw

bedSort Signal.Unique.str2.out.bg Signal.Unique.str2.out.bg
bedGraphToBigWig Signal.Unique.str2.out.bg {chrom_sizes} *_plusUniq.bw


### task rsem_quant

rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed {rnd_seed} -p {ncpus} --no-bam-output --ci-memory {ramGB}000 --forward-prob {fwd_prob} --paired-end \
{anno_bam} rsem_index/rsem {bam_root}_rsem


### task mad_qc

Rscript {path_to_madR} {quants_1} {quants_2}

### task rna_qc

