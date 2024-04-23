### STAR Index
../bin/Linux_x86_64_static/STAR --runThreadN 24 --runMode genomeGenerate --genomeDir ./. --genomeFastaFiles ../../../MacaqueTtoT-test/Macaque.t-to-t/Macaque_Assembly.v3-for-encode.fa --sjdbGTFfile ../../../MacaqueTtoT-test/Macaque.t-to-t/Macaque.t-to-t.sorted.final.withoutMT.gtf

./fastp -i ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R1.fq.gz -I ../../mm10/Brain/rna/RM22050501-P0-RNA/RM22050501-P0-RNA_R2.fq.gz -o ./Brain_P0_R1.trimmed.fq.gz -O ./Brain_P0_R2.trimmed.fq.gz --detect_adapter_for_pe --thread 24 --json Brain_P0_fastp.json --html Brain_P0_fastp.html 2> Brain_P0_fastp.log

### task align
../STAR-2.7.11b/bin/Linux_x86_64_static/STAR  --genomeDir ../STAR-2.7.11b/index/ --readFilesIn ./Brain_P0_R1.trimmed.fq.gz ./Brain_P0_R2.trimmed.fq.gz --readFilesCommand zcat --runThreadN 24 --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1

### task bam_to_signals
../STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph

./bedSort Signal.UniqueMultiple.str1.out.bg Signal.UniqueMultiple.str1.out.sorted.bg
./bedGraphToBigWig Signal.UniqueMultiple.str1.out.sorted.bg Macaca.chrom.sizes Brain_P0_minusAll.bw
./bedSort Signal.UniqueMultiple.str2.out.bg Signal.UniqueMultiple.str2.out.sorted.bg
./bedGraphToBigWig Signal.UniqueMultiple.str2.out.sorted.bg Macaca.chrom.sizes Brain_P0_plusAll.bw

./bedSort Signal.Unique.str1.out.bg Signal.Unique.str1.out.sorted.bg
./bedGraphToBigWig Signal.Unique.str1.out.sorted.bg Macaca.chrom.sizes Brain_P0_minusUniq.bw
./bedSort Signal.Unique.str2.out.bg Signal.Unique.str2.out.sorted.bg
./bedGraphToBigWig Signal.Unique.str2.out.sorted.bg Macaca.chrom.sizes Brain_P0_plusUniq.bw

### task rsem_quant
./rsem-prepare-reference --gtf ../../MacaqueTtoT-test/Macaque.t-to-t/Macaque.t-to-t.sorted.final.withoutMT.gtf ../../MacaqueTtoT-test/Macaque.t-to-t/Macaque_Assembly.v3-for-encode.fa ./Macaca/

../RSEM/rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed 666 -p 24 --no-bam-output --ci-memory 100000 --paired-end Aligned.toTranscriptome.out.bam ../RSEM/Macaca/ Brain_P0_RSEM


### MAD QC for replicates: Rscript {path_to_madR} {quants_1} {quants_2}
