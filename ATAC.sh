### Task Trim Adapter
fastp -i ../../mm10/Brain/atac/RM22050501-Bulk-ATAC-1/RM22050501-Bulk-ATAC-1_R1.fq.gz -I ../../mm10/Brain/atac/RM22050501-Bulk-ATAC-1/RM22050501-Bulk-ATAC-1_R2.fq.gz -o Brain_P0_Rep1_R1.trimmed.fq.gz -O Brain_P0_Rep1_R2.trimmed.fq.gz --detect_adapter_for_pe --thread 16 --json Brain_P0_Rep1_fastp.json --html Brain_P0_Rep1_fastp.html 2 > Brain_P0_Rep1_fastp.log

fastqc -t 20 *_*.fq.gz -o ./QC_result

### Bowtie2
bowtie2 -p 30 -X 2000 --mm -x ../Bowtie2/Index/Macaca/Macaca -1 Brain_P0_Rep1_R1.trimmed.fq.gz  -2 Brain_P0_Rep1_R2.trimmed.fq.gz -S Brain_P0_Rep1_align.sam --met-file Brain_P0_Rep1_align_result.log &> Brain_P0_Rep1_align.log
samtools sort Brain_P0_Rep1_align.sam -@ 20 -o Brain_P0_Rep1_align_sorted.bam
samtools index Brain_P0_Rep1_align_sorted.bam
### Need non-mito BAM？
samtools sort -n Brain_P0_Rep1_align_sorted.bam -O sam | SAMstats --sorted_sam_file - --outf Brain_P0_Rep1_align_sorted.samstat_qc
samtools flagstat Brain_P0_Rep1_align_sorted.bam > Brain_P0_Rep1_align_sorted_stat.log

### rm unmapped lowq reads ${mapq_thresh}=30?
samtools view -F 1804 -f 2 -q 30 -u Brain_P0_Rep1_align_sorted.bam | samtools sort -n /dev/stdin -o Brain_P0_Rep1_align_sorted_filt.bam
samtools fixmate -r Brain_P0_Rep1_align_sorted_filt.bam Brain_P0_Rep1_align_sorted_filt_fixmate.bam
samtools sort Brain_P0_Rep1_align_sorted_filt_fixmate.bam -@ 20 -o Brain_P0_Rep1_align_sorted_filt_fixmate_sorted.bam

##mark dup
java -jar /lustre/home/zhangfy/data0428/picard.jar MarkDuplicates -I Brain_P0_Rep1_align_sorted_filt_fixmate_sorted.bam -O Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_markDup.bam -M Brain_P0_Rep1_align_sorted_markDup.log
samtools sort Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_markDup.bam -@ 20 -o Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_markDup_sorted.bam
samtools index Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_markDup_sorted.bam
samtools sort -n Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_markDup_sorted.bam -O sam | SAMstats --sorted_sam_file - --outf Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_markDup_sorted.samstat_qc
samtools flagstat Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_markDup_sorted.bam > Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_markDup_sorted_stat.log

##rm duplicates
java -jar /lustre/home/zhangfy/data0428/picard.jar MarkDuplicates -I Brain_P0_Rep1_align_sorted_filt_fixmate_sorted.bam -O Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_rmDup.bam --REMOVE_DUPLICATES true -M Brain_P0_Rep1_align_sorted_rmDup.log
samtools sort Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_rmDup.bam -@ 20 -o Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_rmDup_sorted.bam
samtools index Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_rmDup_sorted.bam
samtools sort -n Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_rmDup_sorted.bam -O sam | SAMstats --sorted_sam_file - --outf Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_rmDup_sorted.samstat_qc
samtools flagstat Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_rmDup_sorted.bam > Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_rmDup_sorted_stat.log

##shift +4/-5 移动
alignmentSieve --numberOfProcessors 20 --ATACshift -b Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_rmDup_sorted.bam -o Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_rmDup_sorted_final.bam
samtools sort Brain_P0_Rep1_align_sorted_filt_fixmate_sorted_rmDup_sorted_final.bam -@ 20 -o Brain_P0_Rep1_final.bam
samtools index Brain_P0_Rep1_final.bam
samtools sort -n Brain_P0_Rep1_final.bam -O sam | SAMstats --sorted_sam_file - --outf Brain_P0_Rep1_final.samstat_qc
samtools flagstat Brain_P0_Rep1_final.bam > Brain_P0_Rep1_final.log

###call peak
LC_COLLATE=C bedtools bamtobed -bedpe -mate1 -i Brain_P0_Rep1_final.bam | gzip -nc > Brain_P0_Rep1_final.bedpe.gz
zcat -f Brain_P0_Rep1_final.bedpe.gz | awk \'BEGIN{OFS="\\t"}{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",$1,$2,$3,$9,$4,$5,$6,$10}\' | gzip -nc > Brain_P0_Rep1_final.tagAlign.gz
macs2 callpeak -t Brain_P0_Rep1_final.bam -g 3077605270 -B --keep-dup all --nomodel --shift -100 --extsize 200 -n Brain_P0_Rep1

##Fragment distribution
echo "conda activate R3.6" >> $path1/${file1}.pbs
echo "java -jar /lustre/home/zhangfy/data0428/picard.jar CollectInsertSizeMetrics -H $path3/trim_data/align/${file1}_InsertSize.pdf -I $path3/trim_data/align/${file1}_final_sort.bam -O $path3/trim_data/align/${file1}_InsertSize.txt" >> $path1/${file1}.pbs
echo "conda deactivate" >> $path1/${file1}.pbs
##bam to bw

echo "bamCoverage  --normalizeUsing RPKM --extendReads --binSize 500 -p 4 --bam $path3/trim_data/align/${file1}_final_sort.bam -o $path3/trim_data/align/${file1}_final.bw" >> $path1/${file1}.pbs
##TSS
echo "computeMatrix reference-point --referencePoint TSS -R /lustre/home/zhangfy/data0428/reference/Mmul10_ensembl_TSS.bed -S $path3/trim_data/align/${file1}_final.bw -p 20 -b 3000 -a 3000 --skipZeros -o $path3/trim_data/align/${file1}_final_TSS.gz --outFileSortedRegions $path3/trim_data/align/${file1}_final_TSS.bed --outFileNameMatrix $path3/trim_data/align/${file1}_final_TSS.matirx.txt" >> $path1/${file1}.pbs
echo "plotHeatmap -m $path3/trim_data/align/${file1}_final_TSS.gz -out $path3/trim_data/align/${file1}_final_TSS_heatmap.pdf" >> $path1/${file1}.pbs
##peak center
echo "computeMatrix reference-point --referencePoint center -R $path3/trim_data/peaks/${file1}_summits.bed -S $path3/trim_data/align/${file1}_final.bw -p 20 -b 3000 -a 3000 --skipZeros -o $path3/trim_data/align/${file1}_final_center.gz --outFileSortedRegions $path3/trim_data/align/${file1}_final_center.bed --outFileNameMatrix $path3/trim_data/align/${file1}_final_center.matrix.txt" >> $path1/${file1}.pbs
echo "plotHeatmap -m $path3/trim_data/align/${file1}_final_center.gz -out $path3/trim_data/align/${file1}_final_center_heatmap.pdf" >> $path1/${file1}.pbs
##Fingerprint
echo "plotFingerprint -b $path3/trim_data/align/${file1}_final_sort.bam --numberOfProcessors 20 --labels ${file1} --minMappingQuality 30 --skipZeros --numberOfSamples 500000 -T "Fingerprint" --plotFile $path3/trim_data/align/${file1}_fingerprints.png --outRawCounts $path3/trim_data/align/${file1}.fingerprints.tab" >> $path1/${file1}.pbs
##FRiP
echo "echo \"$file1\" \"total reads\"  > $path3/trim_data/peaks/FRiP_calculation_${file1}.txt" >> $path1/${file1}.pbs
echo "T=\$(cat $path3/trim_data/align/${file1}_final.bed|wc -l)" >> $path1/${file1}.pbs
echo "echo \$T >> $path3/trim_data/peaks/FRiP_calculation_${file1}.txt" >> $path1/${file1}.pbs

echo "echo \"$file1\" \"reads on peak\"  >> $path3/trim_data/peaks/FRiP_calculation_${file1}.txt" >> $path1/${file1}.pbs
echo "P=\$(bedtools intersect -nonamecheck -a $path3/trim_data/align/${file1}_final.bed -b $path3/trim_data/peaks/${file1}_peaks.narrowPeak |wc -l)" >> $path1/${file1}.pbs
echo "echo \$P >> $path3/trim_data/peaks/FRiP_calculation_${file1}.txt" >> $path1/${file1}.pbs

echo "echo \"$file1\" \"FRiP\" >> $path3/trim_data/peaks/FRiP_calculation_${file1}.txt" >> $path1/${file1}.pbs
echo "FRiP=\$(awk \"BEGIN {print \"\$P\"/\"\$T\"}\" )" >> $path1/${file1}.pbs
echo "echo \$FRiP >> $path3/trim_data/peaks/FRiP_calculation_${file1}.txt 2>$path1/log_out" >> $path1/${file1}.pbs
echo "rm $path3/trim_data/align/${file1}_align.sam" >> $path1/${file1}.pbs
qsub $path1/${file1}.pbs
done
