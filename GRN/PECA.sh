#/bin/bash

# PECA2 v3.0.1 updated May 27th 2019
# step 1: call peak from bam file
# step 2: motif binding
# step 3: calculate opn
# step 4: corr+dist
# step 5: score(Exp,binding,Opn,weight)

input=$1
genome="/lustre/home/zhangfy/fetalmonkey/Macaque.t-to-t/Macaque_Assembly.v3-for-encode.fa"

cat ../RNA/${input}.txt > ./Input/${input}.txt
mkdir ./Results/${input}/
cat ../ATAC/Openness/${input}.bed | awk '$2>2' | tr '_' '\t' | sortBed | awk -v OFS='\t' '{print $1"_"$2"_"$3,$4}' > ./Results/${input}/openness.bed
cat ../ATAC/Openness/${input}.bed | awk '$2>2' | tr '_' '\t' | sortBed | awk -v OFS='\t' '{print $1,$2,$3}' > ./Results/${input}/region.bed
cat ../ATAC/Openness/${input}.bed | awk '$2>2' | tr '_' '\t' | sortBed | awk -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' > ./Results/${input}/region.txt
cd ./Results/${input}/

echo step 1: motif binding....
findMotifsGenome.pl region.txt ${genome} ./. -p 16 -size given -find ../../Data/all_motif_rmdup -preparsedDir ../../Homer/ > MotifTarget.bed
cat MotifTarget.bed|awk 'NR>1'|cut -f 1,4,6 > MotifTarget.txt
rm MotifTarget.bed
rm motifFindingParameters.txt

echo step 2: calculate opn...
cat openness.bed |tr '_' '\t' > openness1.bed
bedtools intersect -a openness1.bed -b ../../Prior/Opn_median.bed -wa -wb -sorted|cut -f 1-4,8|sed 's/\t/_/1'|sed 's/\t/_/1'|sed 's/\t/_/1'|awk 'BEGIN{OFS="\t"}{ if ($2>a[$1] ) a[$1]=$2 }END{for (i in a) print i,a[i]}'|sed 's/_/\t/3' > openness2.bed
mkdir Enrichment
cat openness2.bed|awk 'BEGIN{OFS="\t"}{print $1,($2+0.5)/($3+0.5)}'|sort -k2nr|cut -f 1|tr '_' '\t'|awk 'BEGIN{OFS="\t"}{if ($3-$2 < 2000) print $0}'|head -10000 > ./Enrichment/region.bed
cat ../../scr/mf_collect.m > ./Enrichment/mf_collect.m 
cd ./Enrichment/
findMotifsGenome.pl region.bed ${genome} ./. -p 16 -size given -mask -nomotif -mknown ../../../Data/all_motif_rmdup -preparsedDir ../../../Homer/
/lustre/software/matlab/bin/matlab -nodisplay -nosplash -nodesktop -r "mf_collect; exit"
cd ../

echo step 3: Prior....
bedtools intersect -a region.bed -b ../../Prior/RE_gene_corr.bed -wa -wb -sorted|cut -f 1-3,7-9|sed 's/\t/\_/1'|sed 's/\t/\_/1'>peak_gene_100k_corr.bed

echo step 4: Network....
cp ../../scr/mfbs.m ./.
sed "s/toreplace/${input}/g" ../../scr/PECA_network.m > PECA_network.m
/lustre/software/matlab/bin/matlab -nodisplay -nosplash -nodesktop -r "PECA_network; exit"
echo ${input} PECA done
