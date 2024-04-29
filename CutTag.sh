#!/bin/bash

##QC_raw_data
mkdir QC_result
fastqc -t 20 ../../mm10/Brain/cut-tag/P0/CUT-Tag-RM22050502_P0-Brain-H3K27ac-1/*.gz -o ./QC_result/

##trimming fastp
fastp -i ../../mm10/Brain/cut-tag/P0/CUT-Tag-RM22050502_P0-Brain-H3K27ac-1/CUT-Tag-RM22050502_P0-Brain-H3K27ac-1_1.fq.gz -I ../../mm10/Brain/cut-tag/P0/CUT-Tag-RM22050502_P0-Brain-H3K27ac-1/CUT-Tag-RM22050502_P0-Brain-H3K27ac-1_2.fq.gz -o Brain_H3K27ac_P0_Rep1_R1.trimmed.fq.gz -O Brain_H3K27ac_P0_Rep1_R2.trimmed.fq.gz --detect_adapter_for_pe --thread 16 --json Brain_H3K27ac_P0_Rep1_fastp.json --html Brain_H3K27ac_P0_Rep1_fastp.html 2 > Brain_H3K27ac_P0_Rep1_fastp.log
fastqc -t 20 Brain_H3K27ac_P0_Rep1_*.fq.gz -o QC_result
##alignment
bowtie2 -p 30 -X 2000 --mm -x ../Bowtie2/Index/Macaca/Macaca -1 Brain_H3K27ac_P0_Rep1_R1.trimmed.fq.gz  -2 Brain_H3K27ac_P0_Rep1_R2.trimmed.fq.gz -S Brain_H3K27ac_P0_Rep1_align.sam --met-file Brain_H3K27ac_P0_Rep1_align_result.log &> Brain_H3K27ac_P0_Rep1_align.log

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



##marker
#narrow peak
marker1=CTCF
marker3=H3K4me2
marker4=H3K4me3
marker7=H3K27ac
marker5=H3K9ac
##broad peak
marker2=H3K4me1
marker6=H3K9me3
marker8=H3K27me3
marker9=H3K36me3

##IgG
marker10=IgG

###常规窄峰
for monkey in $(ls CUT-Tag*|awk -F "-" '{print $3}'|sort|uniq)
do
        for tissue in $(ls CUT-Tag*|awk -F "-" '{print $4}'|sort|uniq) 
        do
##常规窄峰call peak
                for markerN in $marker1 $marker3 $marker4 $marker5 $marker7
                do
			for rep in 1 2	
			do       
				sample1=CUT-Tag-${monkey}-${tissue}-${markerN}-${rep}  
				if [ ! -d $path1/${sample1} ]
				then echo "$sample1 folder not exist"
				else        
##call peak
           
                	echo "#!/bin/bash" > $path1/$sample1/${sample1}_cuttag_narrow.pbs
					echo "#PBS -N CUT_Tag_narrow " >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
					echo "#PBS -o /lustre/home/zhangfy/Brain/cut-tag/E110/pbs_log/log_pbs.out" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
					echo "#PBS -e /lustre/home/zhangfy/Brain/cut-tag/E110/pbs_log/$log_pbs.err" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
					echo "#PBS -l nodes=cu09:ppn=2" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
					echo "source activate CRE" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs


					mkdir $path2/${sample1}/trim_data/peaks
			

#		if [ ! -d $path1/${sample2} ]
#		then echo "$sample2 folder not exist"
#		else mkdir $path1/${sample2}/trim_data/peaks
#		fi
		
			echo "macs2 callpeak -t $path2/$sample1/trim_data/align/${sample1}_final.sort.bam -c $path2/CUT-Tag-${monkey}-${tissue}-${marker10}-${rep}/trim_data/align/CUT-Tag-${monkey}-${tissue}-${marker10}-${rep}_final.sort.bam -p 1e-3 -g 3077605270 -f BAMPE -B --keep-dup all -n ${sample1} --outdir $path2/$sample1/trim_data/peaks/" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
#                	echo "macs2 callpeak -t $path1/$sample2/trim_data/align/${sample2}_final.sort.bam -c $path1/CUT-Tag-${monkey}-${tissue}-${marker10}-2/trim_data/align/CUT-Tag-${monkey}-${tissue}-${marker10}-2_final.sort.bam -p 1e-3 -g 3077605270 -f BAMPE -B --keep-dup all -n ${sample2} --outdir $path1/$sample2/trim_data/peaks/" >> $path1/${file1}.pbs

##stat peak number
	                echo "echo  \"${sample1}\" >  $path2/$sample1/trim_data/peaks/samples_peak_number_${sample1}.log" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
	                echo "wc -l $path2/$sample1/trim_data/peaks/${sample1}_peaks.narrowPeak >>  $path2/$sample1/trim_data/peaks/samples_peak_number_${sample1}.log" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
##Fragment distribution
					echo "conda activate R3.6" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
	                echo "java -jar /lustre/home/zhangfy/data0428/picard.jar CollectInsertSizeMetrics -H $path2/$sample1/trim_data/align/${sample1}_InsertSize.pdf -I $path2/$sample1/trim_data/align/${sample1}_final.sort.bam -O $path2/$sample1/trim_data/align/${sample1}_InsertSize.txt" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
	                echo "conda deactivate" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
#                echo "java -jar ~/software/picard.jar CollectInsertSizeMetrics -H $path1/$sample2/trim_data/align/${sample2}_InsertSize.pdf -I $path1/$sample2/trim_data/align/${sample2}_final.sort.bam -O $path1/$sample2/trim_data/align/${sample2}_InsertSize.txt" >> $path1/${file1}.pbs
##TSS
			echo "bamCoverage  --normalizeUsing RPKM --extendReads --binSize 500 -p 4 --bam $path2/$sample1/trim_data/align/${sample1}_final.sort.bam -o $path2/$sample1/trim_data/align/${sample1}_final.bw" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
#		echo "bamCoverage  --normalizeUsing RPKM --extendReads --binSize 500 -p 4 --bam $path1/$sample2/trim_data/align/${sample2}_final.sort.bam -o $path1/$sample2/trim_data/align/${sample2}_final.bw" >> $path1/${file1}.pbs
			echo "computeMatrix reference-point --referencePoint TSS -R /lustre/home/zhangfy/data0428/reference/Mmul10_ensembl_TSS.bed -S $path2/$sample1/trim_data/align/${sample1}_final.bw -p 20 -b 3000 -a 3000 --skipZeros -o $path2/$sample1/trim_data/align/${sample1}_final_TSS.gz --outFileSortedRegions $path2/$sample1/trim_data/align/${sample1}_final_TSS.bed --outFileNameMatrix $path2/$sample1/trim_data/align/${sample1}_final_TSS.matrix.txt" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
#		echo "computeMatrix reference-point --referencePoint TSS -R /lustre/home/encode/reference/Macaque_T2T+ensembl_108_Y/Macaca_T2T_Y_TSS.bed -S $path1/$sample2/trim_data/align/${sample2}_final.bw -p 20 -b 3000 -a 3000 --skipZeros -o $path1/$sample2/trim_data/align/${sample2}_final_TSS.gz --outFileSortedRegions $path1/$sample2/trim_data/align/${sample2}_final_TSS.bed" >> $path1/${file1}.pbs
			echo "plotHeatmap -m $path2/$sample1/trim_data/align/${sample1}_final_TSS.gz -out $path2/$sample1/trim_data/align/${sample1}_final_TSS_heatmap.pdf" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
#		echo "plotHeatmap -m $path1/$sample2/trim_data/align/${sample2}_final_TSS.gz -out $path1/$sample2/trim_data/align/${sample2}_final_TSS_heatmap.pdf" >> $path1/${file1}.pbs
##peak center
			echo "computeMatrix reference-point --referencePoint center -R $path2/$sample1/trim_data/peaks/${sample1}_summits.bed -S $path2/$sample1/trim_data/align/${sample1}_final.bw -p 20 -b 3000 -a 3000 --skipZeros -o $path2/$sample1/trim_data/align/${sample1}_final_center.gz --outFileSortedRegions $path2/$sample1/trim_data/align/${sample1}_final_center.bed --outFileNameMatrix $path2/$sample1/trim_data/align/${sample1}_final_center.matrix.txt" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
#		echo "computeMatrix reference-point --referencePoint center -R $path1/$sample2/trim_data/peaks/${sample2}_summits.bed -S $path1/$sample2/trim_data/align/${sample2}_final.bw -p 20 -b 3000 -a 3000 --skipZeros -o $path1/$sample2/trim_data/align/${sample2}_final_center.gz --outFileSortedRegions $path1/$sample2/trim_data/align/${sample2}_final_center.bed" >> $path1/${file1}.pbs

			echo "plotHeatmap -m $path2/$sample1/trim_data/align/${sample1}_final_center.gz -out $path2/$sample1/trim_data/align/${sample1}_final_center_heatmap.pdf" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
#		echo "plotHeatmap -m $path1/$sample2/trim_data/align/${sample2}_final_center.gz -out $path1/$sample2/trim_data/align/${sample2}_final_center_heatmap.pdf" >> $path1/${file1}.pbs
##bam to bed
			echo "bamToBed -i $path2/$sample1/trim_data/align/${sample1}_final.sort.bam > $path2/$sample1/trim_data/align/${sample1}_final.bed" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
#		echo "bamToBed -i $path1/$sample2/trim_data/align/${sample2}_final.sort.bam > $path1/$sample2/trim_data/align/${sample2}_final.bed" >> $path1/${file1}.pbs
##Fingerprint
			echo "plotFingerprint -b $path2/$sample1/trim_data/align/${sample1}_final.sort.bam $path2/CUT-Tag-${monkey}-${tissue}-${marker10}-${rep}/trim_data/align/CUT-Tag-${monkey}-${tissue}-${marker10}-${rep}_final.sort.bam --labels ${sample1} CUT-Tag-${monkey}-${tissue}-${marker10}-${rep} --minMappingQuality 30 --skipZeros --numberOfSamples 500000 -T \"Fingerprint\" --plotFile $path2/$sample1/trim_data/align/${sample1}_fingerprints.png --outRawCounts $path2/$sample1/trim_data/align/${sample1}.fingerprints.tab" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
#		echo "plotFingerprint -b $path1/$sample2/trim_data/align/${sample2}_final.sort.bam $path1/CUT-Tag-${monkey}-${tissue}-${marker10}-2/trim_data/align/CUT-Tag-${monkey}-${tissue}-${marker10}-2_final.sort.bam --labels ${sample2} CUT-Tag-${monkey}-${tissue}-${marker10}-2 --minMappingQuality 30 --skipZeros --numberOfSamples 500000 -T \"Fingerprint\" --plotFile $path1/$sample2/trim_data/align/${sample2}_fingerprints.png --outRawCounts $path1/$sample2/trim_data/align/${sample2}.fingerprints.tab" >> $path1/${file1}.pbs

##FRiP
			echo "echo \"$sample1\" \"total reads\"  > $path2/$sample1/trim_data/peaks/FRiP_calculation_${sample1}.txt" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
			echo "T=\$(cat $path2/$sample1/trim_data/align/${sample1}_final.bed|wc -l)" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
			echo "echo \$T >> $path2/$sample1/trim_data/peaks/FRiP_calculation_${sample1}.txt" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs

			echo "echo \"$sample1\" \"reads on peak\"  >> $path2/$sample1/trim_data/peaks/FRiP_calculation_${sample1}.txt" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
			echo "P=\$(bedtools intersect -nonamecheck -a $path2/$sample1/trim_data/align/${sample1}_final.bed -b $path2/$sample1/trim_data/peaks/${sample1}_peaks.narrowPeak |wc -l)" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
			echo "echo \$P >> $path2/$sample1/trim_data/peaks/FRiP_calculation_${sample1}.txt" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
#			echo "P=\$(bedtools intersect -nonamecheck -a $path1/$sample2/trim_data/align/${sample2}_final.bed -b $path1/$sample2/trim_data/peaks/${sample2}_peaks.broadPeak |wc -l)" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
			echo "echo \"$sample1\" \"FRiP\" >> $path2/$sample1/trim_data/peaks/FRiP_calculation_${sample1}.txt" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
			echo "FRiP=\$(awk \"BEGIN {print \"\$P\"/\"\$T\"}\" )" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
			echo "echo \$FRiP >> $path2/$sample1/trim_data/peaks/FRiP_calculation_${sample1}.txt" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
			echo "mv *.pbs.* ./pbs_log" >> $path1/$sample1/${sample1}_cuttag_narrow.pbs
#######################
		fi
			qsub  $path1/$sample1/${sample1}_cuttag_narrow.pbs
			done
		done
##宽峰call peak
                for markerB in $marker2 $marker6 $marker8 $marker9
                do
                	for rep1 in 1 2
			do
##call peak
                sample2=CUT-Tag-${monkey}-${tissue}-${markerB}-${rep1}
                	if [ ! -d $path1/${sample2} ]
                    then echo "$sample2 folder not exist"
                    else
                	echo "#!/bin/bash" > $path1/$sample2/${sample2}_cuttag_broad.pbs
                	echo "#PBS -N CUT_Tag_broad " >> $path1/$sample2/${sample2}_cuttag_broad.pbs
					echo "#PBS -o /lustre/home/zhangfy/Brain/cut-tag/E110/pbs_log/log_pbs.out" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
					echo "#PBS -e /lustre/home/zhangfy/Brain/cut-tag/E110/pbs_log/log_pbs.err" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
					echo "#PBS -l nodes=cu10:ppn=2" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
                	echo "source activate CRE" >> $path1/$sample2/${sample2}_cuttag_broad.pbs

                   
                    mkdir $path2/${sample2}/trim_data/peaks

                #sample4=CUT-Tag-${monkey}-${tissue}-${markerB}-2
#                echo "mkdir $path1/$sample4/trim_data/peaks" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
                	echo "macs2 callpeak -t $path2/$sample2/trim_data/align/${sample2}_final.sort.bam -c $path2/CUT-Tag-${monkey}-${tissue}-${marker10}-${rep1}/trim_data/align/CUT-Tag-${monkey}-${tissue}-${marker10}-${rep1}_final.sort.bam --broad -p 1e-3 -g 3077605270 -f BAMPE -B --keep-dup all -n ${sample2} --outdir $path2/$sample2/trim_data/peaks/" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#                echo "macs2 callpeak -t $path1/$sample4/trim_data/align/${sample4}_final.sort.bam -c $path1/CUT-Tag-${monkey}-${tissue}-${marker10}-2/trim_data/align/CUT-Tag-${monkey}-${tissue}-${marker10}-2_final.sort.bam --broad -p 1e-3 -g 3077605270 -f BAMPE -B --keep-dup all -n ${sample4} --outdir $path1/$sample4/trim_data/peaks/" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
##stat peak number
                	echo "echo  ${sample2} >  $path2/$sample2/trim_data/peaks/samples_peak_number_${sample2}.log" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
                	echo "wc -l $path2/$sample2/trim_data/peaks/${sample2}_peaks.broadPeak >>  $path2/$sample2/trim_data/peaks/samples_peak_number_${sample2}.log" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#                echo "echo  ${sample4} >  $path1/$sample4/trim_data/peaks/samples_peak_number_${sample4}.log" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#                echo "wc -l $path1/$sample4/trim_data/peaks/${sample4}_peaks.broadPeak >>  $path1/$sample4/trim_data/peaks/samples_peak_number_${sample4}.log" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
##Fragment distribution
					echo "conda activate R3.6" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
                	echo "java -jar /lustre/home/zhangfy/data0428/picard.jar CollectInsertSizeMetrics -H $path2/$sample2/trim_data/align/${sample2}_InsertSize.pdf -I $path2/$sample2/trim_data/align/${sample2}_final.sort.bam -O $path2/$sample2/trim_data/align/${sample2}_InsertSize.txt" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
                	echo "conda deactivate" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#                echo "java -jar ~/software/picard.jar CollectInsertSizeMetrics -H $path1/$sample4/trim_data/align/${sample4}_InsertSize.pdf -I $path1/$sample4/trim_data/align/${sample4}_final.sort.bam -O $path1/$sample4/trim_data/align/${sample4}_InsertSize.txt" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
##TSS
			echo "bamCoverage  --normalizeUsing RPKM --extendReads --binSize 500 -p 4 --bam $path2/$sample2/trim_data/align/${sample2}_final.sort.bam -o $path2/$sample2/trim_data/align/${sample2}_final.bw" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#		echo "bamCoverage  --normalizeUsing RPKM --extendReads --binSize 500 -p 4 --bam $path1/$sample4/trim_data/align/${sample4}_final.sort.bam -o $path1/$sample4/trim_data/align/${sample4}_final.bw" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
			echo "computeMatrix reference-point --referencePoint TSS -R /lustre/home/zhangfy/data0428/reference/Mmul10_ensembl_TSS.bed -S $path2/$sample2/trim_data/align/${sample2}_final.bw -p 20 -b 3000 -a 3000 --skipZeros -o $path2/$sample2/trim_data/align/${sample2}_final_TSS.gz --outFileSortedRegions $path2/$sample2/trim_data/align/${sample2}_final_TSS.bed --outFileNameMatrix $path2/$sample2/trim_data/align/${sample2}_final_TSS.matrix.txt" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#		echo "computeMatrix reference-point --referencePoint TSS -R /lustre/home/encode/reference/Macaque_T2T+ensembl_108_Y/Macaca_T2T_Y_TSS.bed -S $path1/$sample4/trim_data/align/${sample4}_final.bw -p 20 -b 3000 -a 3000 --skipZeros -o $path1/$sample4/trim_data/align/${sample4}_final_TSS.gz --outFileSortedRegions $path1/$sample4/trim_data/align/${sample4}_final_TSS.bed" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
			echo "plotHeatmap -m $path2/$sample2/trim_data/align/${sample2}_final_TSS.gz -out $path2/$sample2/trim_data/align/${sample2}_final_TSS_heatmap.pdf" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#		echo "plotHeatmap -m $path1/$sample4/trim_data/align/${sample4}_final_TSS.gz -out $path1/$sample4/trim_data/align/${sample4}_final_TSS_heatmap.pdf" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
##peak region
			echo "computeMatrix scale-regions -R $path2/$sample2/trim_data/peaks/${sample2}_peaks.broadPeak -S $path2/$sample2/trim_data/align/${sample2}_final.bw --regionBodyLength 2000 --startLabel \"peak start\" --endLabel \"peak end\" -p 20 -b 3000 -a 3000 --skipZeros -o $path2/$sample2/trim_data/align/${sample2}_final_peak_region.gz --outFileSortedRegions $path2/$sample2/trim_data/align/${sample2}_final_peak_region.bed --outFileNameMatrix $path2/$sample2/trim_data/align/${sample2}_final_peak_region.matrix.txt" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#		echo "computeMatrix scale-regions -R $path1/$sample4/trim_data/peaks/${sample4}_peaks.broadPeak -S $path1/$sample4/trim_data/align/${sample4}_final.bw --regionBodyLength 2000 --startLabel \"peak start\" --endLabel \"peak end\" -p 20 -b 3000 -a 3000 --skipZeros -o $path1/$sample4/trim_data/align/${sample4}_final_peak_region.gz --outFileSortedRegions $path1/$sample4/trim_data/align/${sample4}_final_peak_region.bed" >> $path1/$sample2/${sample2}_cuttag_broad.pbs

			echo "plotHeatmap --startLabel \"start\" --endLabel \"end\" -m $path2/$sample2/trim_data/align/${sample2}_final_peak_region.gz -out $path2/$sample2/trim_data/align/${sample2}_final_peak_region_heatmap.pdf" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#		echo "plotHeatmap -m $path1/$sample4/trim_data/align/${sample4}_final_peak_region.gz -out $path1/$sample4/trim_data/align/${sample4}_final_peak_region_heatmap.pdf" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
##bam to bed
			echo "bamToBed -i $path2/$sample2/trim_data/align/${sample2}_final.sort.bam > $path2/$sample2/trim_data/align/${sample2}_final.bed" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#		echo "bamToBed -i $path1/$sample4/trim_data/align/${sample4}_final.sort.bam > $path1/$sample4/trim_data/align/${sample4}_final.bed" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
##Fingerprint
			echo "plotFingerprint -b $path2/$sample2/trim_data/align/${sample2}_final.sort.bam $path2/CUT-Tag-${monkey}-${tissue}-${marker10}-${rep1}/trim_data/align/CUT-Tag-${monkey}-${tissue}-${marker10}-${rep1}_final.sort.bam --labels ${sample2} CUT-Tag-${monkey}-${tissue}-${marker10}-${rep1} --minMappingQuality 30 --skipZeros --numberOfSamples 500000 -T \"Fingerprint\" --plotFile $path2/$sample2/trim_data/align/${sample2}_fingerprints.png --outRawCounts $path2/$sample2/trim_data/align/${sample2}.fingerprints.tab" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#		echo "plotFingerprint -b $path1/$sample4/trim_data/align/${sample4}_final.sort.bam $path1/CUT-Tag-${monkey}-${tissue}-${marker10}-2/trim_data/align/CUT-Tag-${monkey}-${tissue}-${marker10}-2_final.sort.bam --labels ${sample4} CUT-Tag-${monkey}-${tissue}-${marker10}-2 --minMappingQuality 30 --skipZeros --numberOfSamples 500000 -T \"Fingerprint\" --plotFile $path1/$sample4/trim_data/align/${sample4}_fingerprints.png --outRawCounts $path1/$sample4/trim_data/align/${sample4}.fingerprints.tab" >> $path1/$sample2/${sample2}_cuttag_broad.pbs

##FRiP
			echo "echo \"$sample2\" \"total reads\"  > $path2/$sample2/trim_data/peaks/FRiP_calculation_${sample2}.txt" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
			echo "T=\$(cat $path2/$sample2/trim_data/align/${sample2}_final.bed|wc -l)" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
			echo "echo \$T >> $path2/$sample2/trim_data/peaks/FRiP_calculation_${sample2}.txt" >> $path1/$sample2/${sample2}_cuttag_broad.pbs

			echo "echo \"$sample2\" \"reads on peak\"  >> $path2/$sample2/trim_data/peaks/FRiP_calculation_${sample2}.txt" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
			echo "P=\$(bedtools intersect -nonamecheck -a $path2/$sample2/trim_data/align/${sample2}_final.bed -b $path2/$sample2/trim_data/peaks/${sample2}_peaks.broadPeak |wc -l)" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
			echo "echo \$P >> $path2/$sample2/trim_data/peaks/FRiP_calculation_${sample2}.txt" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#			echo "P=\$(bedtools intersect -nonamecheck -a $path1/$sample2/trim_data/align/${sample2}_final.bed -b $path1/$sample2/trim_data/peaks/${sample2}_peaks.broadPeak |wc -l)" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
			echo "echo \"$sample2\" \"FRiP\" >> $path2/$sample2/trim_data/peaks/FRiP_calculation_${sample2}.txt" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
			echo "FRiP=\$(awk \"BEGIN {print \"\$P\"/\"\$T\"}\" )" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
			echo "echo \$FRiP >> $path2/$sample2/trim_data/peaks/FRiP_calculation_${sample2}.txt" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
			echo "mv *.pbs.* ./pbs_log" >> $path1/$sample2/${sample2}_cuttag_broad.pbs
#######################
			fi
			qsub  $path1/$sample2/${sample2}_cuttag_broad.pbs
			done
		done
	done
done


