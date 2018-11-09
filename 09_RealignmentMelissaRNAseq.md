### Need to do some differential expression analyses and get some bigwig files for scnbase uploads


### BIGWIG conversion
```
#/pylon5/mc48o5p/remkv6/Baum/01_GlandRNAseq/08_MelissaRNAAlignment2738Genome
samtools sort  -@ 16 pJ2_Race_3.bam >sortedpJ2_Race_3.bam
bedtools genomecov -ibam sortedpJ2_Race_3.bam -bga -g genome738sl.polished.mitoFixed.fa >sortedpJ2_Race_3.bdg
~/bioawk-master/bioawk -c fastx '{print $name,length($seq)}' genome738sl.polished.mitoFixed.fa >Chr.sizes
~/bedGraphToBigWig sortedpJ2_Race_3.bdg Chr.sizes sortedpJ2_Race_3.bw
#this worked, so all

for f in *bam; do printf "samtools sort  -@ 16 "$f" >sorted"$f"\n bedtools genomecov -ibam sorted"$f" -bga -g genome738sl.polished.mitoFixed.fa
>sorted"$f".bdg\n ~/bedGraphToBigWig sorted"$f".bdg Chr.sizes sorted"$f".bw\n" ; done >bamToBigwig.sh

sh bamToBigwig.sh

#uploaded these to Jbrowse
```
