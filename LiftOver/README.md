This is the script to create UCSC chain file from one specie to another.
**Download these softwares before running**:
```bash
mkdir -p ./bin/
for sf in {axtChain,chainMergeSort,chainNet,chainPreNet,chainSort,chainSwap,faSize,faSplit,faToTwoBit,lastz-1.04.00,liftOver,netChainSubset}
do
    wget -O ./bin/${sf} http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v369/${sf}
    chmod +x ./bin/${sf}
done
```
