OUT_PATHS="./."
wget -O $OUT_PATHS/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
mkdir -p $OUT_PATHS/hg38_split/
./bin/faSplit byname hg38.fa $OUT_PATHS/hg38_split/
for i in `ls $OUT_PATHS/hg38_split/ | grep fa | sed s/.fa//g`
do
	./bin/faSize $OUT_PATHS/hg38_split/${i}.fa -detailed > $OUT_PATHS/hg38_split/${i}.sizes
	./bin/faToTwoBit $OUT_PATHS/hg38_split/${i}.fa $OUT_PATHS/hg38_split/${i}.2bit
	echo $i
done

mkdir -p $OUT_PATHS/MT2T_hg38_axt_out/
for CHR in {1..22} X Y
do
    ./binlastz-1.04.00 $OUT_PATHS/hg38_split/chr${CHR}.fa MT2T.fa \
        --scores=./Q/H38Q.txt E=150 M=254 O=600 T=2 Y=15000 K=4500 L=3000 H=2000 \
        --format=axt --ambiguous=n --ambiguous=iupac --allocate:traceback=23107200 \
        --output=$OUT_PATHS/MT2T_hg38_axt_out/chr${CHR}.axt
    ./bin/axtChain -faT -faQ -linearGap=medium -minScore=3000 $OUT_PATHS/MT2T_hg38_axt_out/chr${CHR}.axt \
        $OUT_PATHS/hg38_split/chr${CHR}.fa \
        MT2T.fa \
        $OUT_PATHS/MT2T_hg38_axt_out/${CHR}.chain
done

./bin/chainMergeSort $OUT_PATHS/MT2T_hg38_axt_out/chr*.chain > $OUT_PATHS/MT2T_hg38_axt_out/all.chain

./bin/chainSort $OUT_PATHS/MT2T_hg38_axt_out/all.chain $OUT_PATHS/MT2T_hg38_axt_out/all.sorted.chain

./bin/faSize -detailed hg38.fa | grep -E "^chr([1-9]|1[0-9]|2[0-2]|X|Y)\s" > $OUT_PATHS/hg38.chrom.sizes
./bin/faSize -detailed MT2T.fa > $OUT_PATHS/MT2T.chrom.sizes

./bin/chainPreNet $OUT_PATHS/MT2T_hg38_axt_out/all.sorted.chain $OUT_PATHS/hg38.chrom.sizes $OUT_PATHS/MT2T.chrom.sizes $OUT_PATHS/MT2T_hg38_axt_out/all.pre.chain

./bin/chainNet $OUT_PATHS/MT2T_hg38_axt_out/all.pre.chain $OUT_PATHS/hg38.chrom.sizes $OUT_PATHS/MT2T.chrom.sizes $OUT_PATHS/MT2T_hg38_axt_out/hg38.net $OUT_PATHS/MT2T_hg38_axt_out/MT2T.net

./bin/netChainSubset $OUT_PATHS/MT2T_hg38_axt_out/hg38.net $OUT_PATHS/MT2T_hg38_axt_out/all.pre.chain $OUT_PATHS/MT2T_hg38_axt_out/lift.chain

./bin/chainSort $OUT_PATHS/MT2T_hg38_axt_out/lift.chain $OUT_PATHS/MT2T_hg38_axt_out/hg38ToMT2T.over.chain
gzip $OUT_PATHS/MT2T_hg38_axt_out/hg38ToMT2T.over.chain
cp $OUT_PATHS/hg38ToMT2T.over.chain.gz hg38ToMT2T.over.chain.gz

./bin/chainSort $OUT_PATHS/MT2T_hg38_axt_out/lift.chain $OUT_PATHS/MT2T_hg38_axt_out/hg38ToMT2T.over.chain

./bin/chainSwap $OUT_PATHS/MT2T_hg38_axt_out/hg38ToMT2T.over.chain $OUT_PATHS/MT2T_hg38_axt_out/MT2TToHg38.over.chain
./bin/chainSort $OUT_PATHS/MT2T_hg38_axt_out/MT2TToHg38.over.chain $OUT_PATHS/MT2T_hg38_axt_out/MT2TToHg38.sorted.chain
gzip $OUT_PATHS/MT2T_hg38_axt_out/MT2TToHg38.sorted.chain
cp $OUT_PATHS/MT2T_hg38_axt_out/MT2TToHg38.sorted.chain.gz MT2TToHg38.over.chain.gz
