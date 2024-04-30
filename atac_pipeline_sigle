## QC raw
mkdir ./QC_result
fastqc -t 20 ../../MacaqueTtoT-test/Brain/atac/RM22050501-Bulk-ATAC-1/*.fq.gz -o ./QC_result

## Trim Adapter
fastp -i ../../MacaqueTtoT-test/Brain/atac/RM22050501-Bulk-ATAC-1/RM22050501-Bulk-ATAC-1_R1.fq.gz -I ../../MacaqueTtoT-test/Brain/atac/RM22050501-Bulk-ATAC-1/RM22050501-Bulk-ATAC-1_R2.fq.gz -o Brain_P0_Rep1_R1.trimmed.fq.gz -O Brain_P0_Rep1_R2.trimmed.fq.gz --detect_adapter_for_pe --thread 16 --json Brain_P0_Rep1_fastp.json --html Brain_P0_Rep1_fastp.html 2> Brain_P0_Rep1_fastp.log

## QC trim
fastqc -t 20 *_*.fq.gz -o ./QC_result

## Bowtie2
bowtie2 -p 30 -X 2000 --mm -x ../Bowtie2/Index/Macaca/Macaca -1 Brain_P0_Rep1_R1.trimmed.fq.gz  -2 Brain_P0_Rep1_R2.trimmed.fq.gz -S Brain_P0_Rep1_align.sam --met-file Brain_P0_Rep1_align_result.log &> Brain_P0_Rep1_align.log
samtools sort Brain_P0_Rep1_align.sam -@ 24 -o Brain_P0_Rep1_align_sorted.bam
samtools index Brain_P0_Rep1_align_sorted.bam
samtools flagstat Brain_P0_Rep1_align_sorted.bam > Brain_P0_Rep1_align_sorted_stat.log

## rm mito
samtools view -b -L ../bin/Macaca.chrom_nomito.bed Brain_P0_Rep1_align_sorted.bam > Brain_P0_Rep1_align_nomito.bam
samtools sort Brain_P0_Rep1_align_nomito.bam -@ 24 -o Brain_P0_Rep1_align_nomito_sorted.bam
samtools index Brain_P0_Rep1_align_nomito_sorted.bam
samtools flagstat Brain_P0_Rep1_align_nomito_sorted.bam > Brain_P0_Rep1_align_nomito_sorted_stat.log

### rm unmapped lowq reads ${mapq_thresh}=30?
samtools view -F 1804 -f 2 -q 30 -u Brain_P0_Rep1_align_nomito_sorted.bam | samtools sort -n /dev/stdin -o Brain_P0_Rep1_align_nomito_filt.bam
samtools fixmate -r Brain_P0_Rep1_align_nomito_filt.bam Brain_P0_Rep1_align_nomito_filt_fixmate.bam
samtools sort Brain_P0_Rep1_align_nomito_filt_fixmate.bam -@ 20 -o Brain_P0_Rep1_align_nomito_filt_fixmate_sorted.bam
samtools index Brain_P0_Rep1_align_nomito_filt_fixmate_sorted.bam
samtools flagstat Brain_P0_Rep1_align_nomito_filt_fixmate_sorted.bam > Brain_P0_Rep1_align_nomito_filt_fixmate_sorted_stat.log

## mark dup
java -jar /lustre/home/zhangfy/data0428/picard.jar MarkDuplicates -I Brain_P0_Rep1_align_nomito_filt_fixmate_sorted.bam -O Brain_P0_Rep1_align_nomito_filt_fixmate_markDup.bam -M Brain_P0_Rep1_align_nomito_filt_fixmate_markDup.log
samtools sort Brain_P0_Rep1_align_nomito_filt_fixmate_markDup.bam -@ 20 -o Brain_P0_Rep1_align_nomito_filt_fixmate_markDup_sorted.bam
samtools index Brain_P0_Rep1_align_nomito_filt_fixmate_markDup_sorted.bam
samtools flagstat Brain_P0_Rep1_align_nomito_filt_fixmate_markDup_sorted.bam > Brain_P0_Rep1_align_nomito_filt_fixmate_markDup_sorted_stat.log

## rm duplicates
java -jar /lustre/home/zhangfy/data0428/picard.jar MarkDuplicates -I Brain_P0_Rep1_align_nomito_filt_fixmate_sorted.bam -O Brain_P0_Rep1_align_nomito_filt_fixmate_rmDup.bam --REMOVE_DUPLICATES true -M Brain_P0_Rep1_align_nomito_filt_fixmate_rmDup.log
samtools sort Brain_P0_Rep1_align_nomito_filt_fixmate_rmDup.bam -@ 20 -o Brain_P0_Rep1_align_nomito_filt_fixmate_rmDup_sorted.bam
samtools index Brain_P0_Rep1_align_nomito_filt_fixmate_rmDup_sorted.bam
samtools flagstat Brain_P0_Rep1_align_nomito_filt_fixmate_rmDup_sorted.bam > Brain_P0_Rep1_align_nomito_filt_fixmate_rmDup_sorted_stat.log

## shift +4/-5
alignmentSieve --numberOfProcessors 20 --ATACshift -b Brain_P0_Rep1_align_nomito_filt_fixmate_rmDup_sorted.bam -o Brain_P0_Rep1_align_nomito_filt_fixmate_rmDup_shifted.bam
samtools sort Brain_P0_Rep1_align_nomito_filt_fixmate_rmDup_shifted.bam -@ 20 -o Brain_P0_Rep1_final.bam
samtools index Brain_P0_Rep1_final.bam
samtools flagstat Brain_P0_Rep1_final.bam > Brain_P0_Rep1_final.log

### call peak
samtools sort Brain_P0_Rep1_final.bam -n -@ 24 -o Brain_P0_Rep1.bam
bedtools bamtobed -bedpe -mate1 -i Brain_P0_Rep1.bam | gzip -nc > Brain_P0_Rep1.bedpe.gz
zcat -f Brain_P0_Rep1.bedpe.gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > Brain_P0_Rep1.tagAlign.gz
macs2 callpeak -t Brain_P0_Rep1.tagAlign.gz -g 3077605270 -B --keep-dup all --nomodel --shift -100 --extsize 200 -n Brain_P0_Rep1 --SPMR

### create fc bw file
macs2 bdgcmp -t Brain_P0_Rep1_treat_pileup.bdg -c Brain_P0_Rep1_control_lambda.bdg --o-prefix Brain_P0_Rep1 -m FE
bedtools slop -i Brain_P0_Rep1_FE.bdg -g /lustre/home/zhangfy/Pipeline/bin/Macaca.chrom.sizes -b 0 | bedClip stdin /lustre/home/zhangfy/Pipeline/bin/Macaca.chrom.sizes Brain_P0_Rep1.fc.signal.bedgraph
sort -k1,1 -k2,2n Brain_P0_Rep1.fc.signal.bedgraph | awk 'BEGIN{OFS="\\t"}{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {print $0}; prev_chr=$1; prev_chr_e=$3;}' > Brain_P0_Rep1.fc.signal.srt.bedgraph
bedGraphToBigWig Brain_P0_Rep1.fc.signal.srt.bedgraph /lustre/home/zhangfy/Pipeline/bin/Macaca.chrom.sizes Brain_P0_Rep1.fc.signal.bigwig

### create pvalue bw file
sval=`zcat Brain_P0_Rep1.tagAlign.gz | wc -l`
sval=`expr $sval / 1000000`
macs2 bdgcmp -t Brain_P0_Rep1_treat_pileup.bdg -c Brain_P0_Rep1_control_lambda.bdg --o-prefix Brain_P0_Rep1 -m ppois -S ${sval}
bedtools slop -i Brain_P0_Rep1_ppois.bdg -g /lustre/home/zhangfy/Pipeline/bin/Macaca.chrom.sizes -b 0 | bedClip stdin /lustre/home/zhangfy/Pipeline/bin/Macaca.chrom.sizes Brain_P0_Rep1.pval.signal.bedgraph
sort -k1,1 -k2,2n Brain_P0_Rep1.pval.signal.bedgraph | awk 'BEGIN{{OFS="\\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {{print $0}}; prev_chr=$1; prev_chr_e=$3;}}' > Brain_P0_Rep1.pval.signal.srt.bedgraph
bedGraphToBigWig Brain_P0_Rep1.pval.signal.srt.bedgraph /lustre/home/zhangfy/Pipeline/bin/Macaca.chrom.sizes Brain_P0_Rep1.pval.signal.bigwig

## Fragment distribution
conda activate R3.6
java -jar /lustre/home/zhangfy/data0428/picard.jar CollectInsertSizeMetrics -H $path3/trim_data/align/${file1}_InsertSize.pdf -I $path3/trim_data/align/${file1}_final_sort.bam -O $path3/trim_data/align/${file1}_InsertSize.txt
conda deactivate
