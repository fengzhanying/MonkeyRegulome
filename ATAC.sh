### Task Trim Adapter
fastp -i *_1.fq.gz -I *_2.fq.gz -o *_1.trimmed.fq.gz -O *_2.trimmed.fq.gz --detect_adapter_for_pe --thread 16 --json *_fastp.json --html *_fastp.html 2 > *_fastp.log

fastqc -t 20 *_*.fq.gz -o */QC_result

### Bowtie2
bowtie2 -p 30 -X 2000 --mm -x */Macaca_index -1 *_1.trimmed.fq.gz  -2 *_2.trimmed.fq.gz -S *_align.sam --met-file *_align_result.log &> *_align.log
samtools sort *_align.sam -@ 20 -o *_align_sorted.bam
samtools index *_align_sorted.bam
### Need non-mito BAM？
samtools sort -n *_align_sorted.bam -O sam | SAMstats --sorted_sam_file - --outf *_align_sorted.samstat_qc
samtools flagstat *_align_sorted.bam > *_align_sorted_stat.log

### rm unmapped lowq reads ${mapq_thresh}=30?
samtools view -F 1804 -f 2 -q ${mapq_thresh} -u *_align_sorted.bam | samtools sort -n /dev/stdin -o *_align_sorted_filt.bam
samtools fixmate -r *_align_sorted_filt.bam *_align_sorted_filt_fixmate.bam

##mark dup
java -jar /lustre/home/zhangfy/data0428/picard.jar MarkDuplicates -I *_align_sorted_filt_fixmate.bam -O *_align_sorted_filt_fixmate_markDup.bam -M *_align_sorted_markDup.log
samtools sort *_align_sorted_filt_fixmate_markDup.bam -@ 20 -o *_align_sorted_filt_fixmate_markDup_sorted.bam
samtools index *_align_sorted_filt_fixmate_markDup_sorted.bam
samtools sort -n *_align_sorted_filt_fixmate_markDup_sorted.bam -O sam | SAMstats --sorted_sam_file - --outf *_align_sorted_filt_fixmate_markDup_sorted.samstat_qc
samtools flagstat *_align_sorted_filt_fixmate_markDup_sorted.bam > *_align_sorted_filt_fixmate_markDup_sorted_stat.log

##rm duplicates
java -jar /lustre/home/zhangfy/data0428/picard.jar MarkDuplicates -I *_align_sorted_filt_fixmate.bam -O *_align_sorted_filt_fixmate_rmDup.bam --REMOVE_DUPLICATES true -M *_align_sorted_rmDup.log
samtools sort *_align_sorted_filt_fixmate_rmDup.bam -@ 20 -o *_align_sorted_filt_fixmate_rmDup_sorted.bam
samtools index *_align_sorted_filt_fixmate_rmDup_sorted.bam
samtools sort -n *_align_sorted_filt_fixmate_rmDup_sorted.bam -O sam | SAMstats --sorted_sam_file - --outf *_align_sorted_filt_fixmate_rmDup_sorted.samstat_qc
samtools flagstat *_align_sorted_filt_fixmate_rmDup_sorted.bam > *_align_sorted_filt_fixmate_rmDup_sorted_stat.log

##shift +4/-5 移动
alignmentSieve --numberOfProcessors 20 --ATACshift -b *_align_sorted_filt_fixmate_rmDup_sorted.bam -o *_align_sorted_filt_fixmate_rmDup_sorted_final.bam
samtools sort *_align_sorted_filt_fixmate_rmDup_sorted_final.bam -@ 20 -o *_final.bam
amtools index *_final.bam
samtools sort -n *_final.bam -O sam | SAMstats --sorted_sam_file - --outf *_final.samstat_qc
samtools flagstat *_final.bam > *_final.log

###call peak
echo "mkdir $path3/trim_data/peaks" >> $path1/${file1}.pbs
echo "macs2 callpeak -t $path3/trim_data/align/${file1}_final.bed -g 3077605270 -B --keep-dup all --nomodel --shift -100 --extsize 200 -n ${file1} --outdir $path3/trim_data/peaks/ " >> $path1/${file1}.pbs
##stat peak number
echo "echo  \"${file1}\" >  $path3/trim_data/peaks/samples_peak_number_${file1}.log" >> $path1/${file1}.pbs
echo "wc -l $path3/trim_data/peaks/${file1}_peaks.narrowPeak >> $path3/trim_data/peaks/samples_peak_number_${file1}.log" >> $path1/${file1}.pbs
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
