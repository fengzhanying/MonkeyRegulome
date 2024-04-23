fq1=$1
fq2=$2
sample=$3

fastp -i ${fq1} -I ${fq2} -o ./${sample}_R1.trimmed.fq.gz -O ./${sample}_R2.trimmed.fq.gz --detect_adapter_for_pe --thread 24 --json ${sample}_fastp.json --html ${sample}_fastp.html 2> ${sample}_fastp.log

### task align
STAR  --genomeDir ../STAR-2.7.11b/index/ --readFilesIn ./${sample}_R1.trimmed.fq.gz ./${sample}_R2.trimmed.fq.gz --readFilesCommand zcat --runThreadN 24 --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1

### task bam_to_signals
STAR --runMode inputAlignmentsFromBAM --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph

bedSort Signal.UniqueMultiple.str1.out.bg Signal.UniqueMultiple.str1.out.sorted.bg
bedGraphToBigWig Signal.UniqueMultiple.str1.out.sorted.bg /lustre/home/zhangfy/Pipeline/bin/Macaca.chrom.sizes ${sample}_minusAll.bw
bedSort Signal.UniqueMultiple.str2.out.bg Signal.UniqueMultiple.str2.out.sorted.bg
bedGraphToBigWig Signal.UniqueMultiple.str2.out.sorted.bg /lustre/home/zhangfy/Pipeline/bin/Macaca.chrom.sizes ${sample}_plusAll.bw

bedSort Signal.Unique.str1.out.bg Signal.Unique.str1.out.sorted.bg
bedGraphToBigWig Signal.Unique.str1.out.sorted.bg /lustre/home/zhangfy/Pipeline/bin/Macaca.chrom.sizes ${sample}_minusUniq.bw
bedSort Signal.Unique.str2.out.bg Signal.Unique.str2.out.sorted.bg
bedGraphToBigWig Signal.Unique.str2.out.sorted.bg /lustre/home/zhangfy/Pipeline/bin/Macaca.chrom.sizes ${sample}_plusUniq.bw

### task rsem_quant
rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed 666 -p 24 --no-bam-output --ci-memory 100000 --paired-end Aligned.toTranscriptome.out.bam ../RSEM/Macaca/ ${sample}_RSEM


### MAD QC for replicates: Rscript {path_to_madR} {quants_1} {quants_2}
