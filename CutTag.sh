#!/bin/bash
#change the path
cd /lustre/home/zhangfy/Brain/cut-tag/E110
mkdir pbs_log
mkdir result
path2=/lustre/home/zhangfy/Brain/cut-tag/E110/result
for file1 in CUT-Tag-*
do
path1=$PWD/${file1}
mkdir $path2/${file1}
path3=$path2/${file1}
echo "#!/bin/bash" > $path1/${file1}.pbs
echo "#PBS -N CUT_Tag " >> $path1/${file1}.pbs
echo "#PBS -o /lustre/home/zhangfy/Brain/cut-tag/E110/pbs_log/log_pbs.out" >> $path1/${file1}.pbs
echo "#PBS -e /lustre/home/zhangfy/Brain/cut-tag/E110/pbs_log/log_pbs.err" >> $path1/${file1}.pbs
echo "#PBS -l nodes=cu14:ppn=2" >> $path1/${file1}.pbs
echo "source activate CRE" >> $path1/${file1}.pbs
##QC_raw_data
echo "mkdir $path3/QC_result" >> $path1/${file1}.pbs
echo "fastqc -t 20 $PWD/$file1/${file1}_* -o $path1/QC_result " >> $path1/${file1}.pbs
##trimming
echo "mkdir $path3/trim_data" >> $path1/${file1}.pbs
##fastp
echo "fastp -i $PWD/$file1/${file1}_1.fq.gz -I $PWD/$file1/${file1}_2.fq.gz -o $path3/trim_data/${file1}_1.trimmed.fq.gz -O $path3/trim_data/${file1}_2.trimmed.fq.gz --detect_adapter_for_pe --thread 16 --json $path3/trim_data/${file1}_fastp.json --html $path3/trim_data/${file1}_fastp.html 2>$path3/trim_data/${file1}_fastp.log " >> $path1/${file1}.pbs 
#echo "trim_galore -j 30 --quality 25 --stringency 4 --length 30 -e 0.1 --paired --phred33 --gzip --output_dir $path1/trim_data/ $PWD/$file1/${file1}_1.fq.gz $PWD/$file1/${file1}_2.fq.gz" >> $path1/${file1}.pbs
##QC_clean_data
echo "fastqc -t 20 $path3/trim_data/${file1}_*.fq.gz -o $path3/QC_result" >> $path1/${file1}.pbs
##alignment
echo "mkdir $path3/trim_data/align" >> $path1/${file1}.pbs
echo "bowtie2 -p 30 -I 10 -X 700 --phred33 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -x /lustre/home/zhangfy/data0428/reference/Macaca_index  -1 $path3/trim_data/${file1}_1.trimmed.fq.gz  -2 $path3/trim_data/${file1}_2.trimmed.fq.gz -S $path3/trim_data/align/${file1}_align.sam --met-file $path3/trim_data/align/${file1}_align_result.log &> $path3/trim_data/align/${file1}_align.log" >> $path1/${file1}.pbs
echo "samtools sort -@ 20 -o $path3/trim_data/align/${file1}_align_sorted.bam $path3/trim_data/align/${file1}_align.sam" >> $path1/${file1}.pbs
echo "samtools index $path3/trim_data/align/${file1}_align_sorted.bam" >> $path1/${file1}.pbs
echo "samtools flagstat $path3/trim_data/align/${file1}_align_sorted.bam > $path3/trim_data/align/${file1}_align_sorted_stat.log" >> $path1/${file1}.pbs
##mark_dup
echo "java -jar /lustre/home/zhangfy/data0428/picard.jar MarkDuplicates -I $path3/trim_data/align/${file1}_align_sorted.bam -O $path3/trim_data/align/${file1}_align_sorted_markDup.bam -M $path3/trim_data/align/${file1}_align_sorted_markDup.log" >> $path1/${file1}.pbs
echo "samtools sort -@ 20 -o $path3/trim_data/align/${file1}_align_sorted_markDup.sort.bam $path3/trim_data/align/${file1}_align_sorted_markDup.bam" >> $path1/${file1}.pbs
echo "samtools index $path3/trim_data/align/${file1}_align_sorted_markDup.sort.bam" >> $path1/${file1}.pbs
echo "samtools flagstat $path3/trim_data/align/${file1}_align_sorted_markDup.sort.bam > $path3/trim_data/align/${file1}_align_sorted_markDup_stat.log" >> $path1/${file1}.pbs
##rm duplicates
echo "java -jar /lustre/home/zhangfy/data0428/picard.jar MarkDuplicates -I $path3/trim_data/align/${file1}_align_sorted.bam -O $path3/trim_data/align/${file1}_align_sorted_rmDup.bam --REMOVE_DUPLICATES true -M $path3/trim_data/align/${file1}_align_sorted_rmDup.log" >> $path1/${file1}.pbs
echo "samtools sort -@ 20 -o $path3/trim_data/align/${file1}_align_sorted_rmDup.sort.bam $path3/trim_data/align/${file1}_align_sorted_rmDup.bam" >> $path1/${file1}.pbs
echo "samtools index $path3/trim_data/align/${file1}_align_sorted_rmDup.sort.bam" >> $path1/${file1}.pbs
echo "samtools flagstat $path3/trim_data/align/${file1}_align_sorted_rmDup.sort.bam > $path3/trim_data/align/${file1}_align_sorted_rmDup_stat.log" >> $path1/${file1}.pbs
##rm low quality reads and the reads that not map to the same chromosome
echo "samtools view -h -q 30 -f 2 -b $path3/trim_data/align/${file1}_align_sorted_rmDup.sort.bam > $path3/trim_data/align/${file1}_final.bam" >> $path1/${file1}.pbs
echo "samtools sort -@ 20 -o $path3/trim_data/align/${file1}_final.sort.bam $path3/trim_data/align/${file1}_final.bam" >> $path1/${file1}.pbs
echo "samtools index $path3/trim_data/align/${file1}_final.sort.bam" >> $path1/${file1}.pbs
echo "samtools flagstat $path3/trim_data/align/${file1}_final.sort.bam > $path3/trim_data/align/${file1}_final_stat.log" >> $path1/${file1}.pbs
echo "rm $path3/trim_data/align/${file1}_align.sam" >> $path1/${file1}.pbs
echo "mv *.pbs.* ./pbs_log" >> $path1/${file1}.pbs
qsub  $path1/${file1}.pbs
done



