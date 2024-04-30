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

## TSS enrichment: TSSEscore


## IDR values
