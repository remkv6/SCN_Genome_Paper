# Helitron prediction, an attempt to see how tandem duplications were made
```
There are only two known transposon types that can create duplications as large as 31kb, Mavericks and Helitrons. We do not have the characteristics of mavericks, so helitrons was the best bet.
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/54_helitronScanner
Run Helitron Scanner


#scan with 3' end homology for helitron ends
java -jar  HelitronScanner/HelitronScanner.jar scanTail -lf TrainingSet/training.set.tail90.N.fa -g genome738.mitofixed.noquiver.fa -th 16 -o tail.helitronscanner.738.out
#scan with 5' end homology for helitron ends
java -jar  HelitronScanner/HelitronScanner.jar scanHead -lf TrainingSet/training.set.head90.N.fa -g genome738.mitofixed.noquiver.fa -th 16 -o head.helitronscanner.738.out
#pairs the two ends
java -jar  HelitronScanner/HelitronScanner.jar pairends -hs head.helitronscanner.738.out -ts tail.helitronscanner.738.out -hlr 200:33000 -o paired.helitrons
#creates a fasta file with all helitrons identified
java -jar  HelitronScanner/HelitronScanner.jar draw -p paired.helitrons -g ../../CamTechGenomeComparison/18_mitochondria/swapMitoScaffolds/genome738sl.polished.mitoFixed.fa -o draw_helitrons_pure --pure
```

### Abundance of helitrons
```
#abundance of helitrons


grep -v ">" paired.helitrons |sed 's/ /\n/g' |grep -v "]" |sort -nr |sed 's/:/\t/g' |awk '{if($2>$1) print $2-$1}' |summary.sh

Total:  2,260,901
Count:  136
Mean:   16,624
Median: 16,530
Min:    261
Max:    32,082

#Create Gff and compare


less draw_helitrons_pure.hel.fa |grep ">" |awk '{print $1}' |sed 's/_/\t/g'|sed 's/>//g'|awk '{print $1,"HScanner","pure_helitron",$3,".","+","."}' |sed 's/-/\t/g'|tr " " "\t" |sort -k1,1 -k4,5n >helitrons.gff

bedtools intersect -wo -a helitrons.gff -b ../50_tandemMummer/same.scaf.deduplicated.gff |awk '{print $9,$12,$13}' |sort|uniq|wc
They overlap with 28 mummer duplicates
bedtools intersect -wo -a helitrons.gff -b ../50_tandemMummer/same.scaf.deduplicated.gff |awk '{print $1,$4,$5,$9,$12,$13}' |less
```
```
000001 727997 756088 000001 740000 745000
000001 727997 756088 000001 734000 739000
000063 291068 298412 000063 294000 329000
000063 328066 357428 000063 294000 329000
000063 328066 357428 000063 337000 371000
000138 821030 835773 000138 805000 837000
000154 196021 227876 000154 209000 214000
000154 196021 227876 000154 217000 233000
000220 208195 228738 000220 186000 215000
000250K 74032 100215 000250K 88000 95000
000250K 74032 100215 000250K 76000 83000
000275 260354 268922 000275 246000 284000
000275 274160 274796 000275 246000 284000
000275 298593 307163 000275 284000 322000
000275 312405 313042 000275 284000 322000
000286 46984 78728 000286 24000 47000
000286 46984 78728 000286 47000 53000
000324 58187 79218 000324 67000 76000
000324 58187 79218 000324 49000 59000
000324 58187 79218 000324 76000 84000
000324 394311 414550 000324 412000 419000
000324 394311 414550 000324 402000 412000
000324 394311 414550 000324 385000 395000
000331 27905 43441 000331 11000 36000
000331 27905 43441 000331 36000 57000
000396 96999 118549 000396 108000 114000
000585K 1788 33460 000585K 28000 39000
000585K 35800 41552 000585K 39000 47000
000585K 35800 41552 000585K 28000 39000
000688 121960 130863 000688 127000 154000
000688 161532 162794 000688 159000 186000
000929K 7111 14878 000929K 13000 25000
000929K 19343 24346 000929K 13000 25000
```
### Other overlaps?
```
bedtools intersect -wo -a helitrons.gff -b ../29_effectorMapping/effector.gmapped.gff3 |less
0
bedtools intersect -wo -a helitrons.gff -b ../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '{print $18}'|sort |uniq -c |sort -k1,1nr|less

164 "Motif:A-rich"
 100 "Motif:GA-rich"
  69 "Motif:rnd-3_family-228"
  53 "Motif:(C)n"
  51 "Motif:rnd-4_family-268"
  50 "Motif:rnd-4_family-352"
  49 "Motif:rnd-4_family-265"
  48 "Motif:rnd-5_family-1435"
  47 "Motif:rnd-3_family-980"
  44 "Motif:rnd-4_family-1041"


  #RepeatExplorer overlap?
  bedtools intersect -wo -a helitrons.gff -b ../43_RepeatExpClusters/repeatexplorer.gff3 |awk '{print $17}'|grep "ID=" |sed 's/Name=/\t/g'| sed 's/Contig.*//g' |awk '{print $2}' |sort|uniq -c |sort -k1,1nr |less
  Nothing significant

  #how much gene overlap
  bedtools intersect -wo -a helitrons.gff -b ../50_tandemMummer/augustus_renamed_with_putative_function.gff3 |awk '$11=="mRNA"' |wc
  564



  #how many helitron and duplication relationships are there?
  bedtools intersect -wo -a helitrons.gff -b ../40_tandemdups/split/sorted.final.renamed.gff |awk '{print $1,$3,$4,$9,$12,$13}' |less
  bedtools intersect -wo -a helitrons.gff -b ../40_tandemdups/split/sorted.final.renamed.gff |awk '{print $1,$3,$4,$9,$12,$13}' |wc
  97
  #how many are uniq helitrons?
  bedtools intersect -wo -a helitrons.gff -b ../40_tandemdups/split/sorted.final.renamed.gff |awk '{print $9,$12,$13}' |sort|uniq|wc

  #what are the annnotations?
  bedtools intersect -wo -a helitrons.gff -b ../50_tandemMummer/augustus_renamed_with_putative_function.gff3 |awk '$11=="mRNA"' |cut -f 17 |cut -d " " -f 3- |tr " " "\n" |sort|uniq -c |sort -k1,1nr|less
  bedtools intersect -wo -a helitrons.gff -b ../50_tandemMummer/augustus_renamed_with_putative_function.gff3 |awk '$11=="mRNA"' |cut -f 17 |cut -d " " -f 3- |grep "helicase" |wc
  9
  bedtools intersect -wo -a helitrons.gff -b ../50_tandemMummer/augustus_renamed_with_putative_function.gff3 |awk '$11=="mRNA"' |cut -f 17 |cut -d " " -f 3- |grep "rep"bedtools intersect -wo -a helitrons.gff -b ../50_tandemMummer/augustus_renamed_with_putative_function.gff3 |awk '$11=="mRNA"' |cut -f 17 |cut -d " " -f 3- |grep "helicase" |wc
  9
  bedtools intersect -wo -a helitrons.gff -b ../50_tandemMummer/augustus_renamed_with_putative_function.gff3 |awk '$11=="mRNA"' |cut -f 17 |cut -d " " -f 3- |grep "replication" |wc
  3
  bedtools intersect -wo -a helitrons.gff -b ../50_tandemMummer/augustus_renamed_with_putative_function.gff3 |awk '$11=="mRNA"' |cut -f 17 |cut -d " " -f 3- |grep "Zinc" |wc
  2
  bedtools intersect -wo -a helitrons.gff -b ../50_tandemMummer/augustus_renamed_with_putative_function.gff3 |awk '$11=="mRNA"' |cut -f 17 |cut -d " " -f 3- |grep "finger"|wc
  3
  bedtools intersect -wo -a helitrons.gff -b ../50_tandemMummer/augustus_renamed_with_putative_function.gff3 |awk '$11=="mRNA"' |cut -f 17 |cut -d " " -f 3- |grep "unknown" |wc
  199

  #synteny overlap?
  #how many helitrons
  bedtools intersect -wo -a helitrons.gff -b ../46_SCN_syntenic_tracks/all.synteny.gff |awk '{print $1,$4,$5}' |sort |uniq |wc
  63
  #How many syntenic regions in the genome?
  bedtools intersect -wo -a helitrons.gff -b ../46_SCN_syntenic_tracks/all.synteny.gff |awk '{print $9,$12,$13}' |sort |uniq |wc
  135
  #how many are shared among species
  bedtools intersect -wo -a helitrons.gff -b ../46_SCN_syntenic_tracks/all.synteny.gff |awk '{print $17}' |sort |uniq |wc
  142

  #So 7 genomic regions with helitrons are syntenic to multiple species.   
```

### Retesting helitronscanner with default settings
```
I reran this under default settings to be sure I was not biasing.  
java -jar  HelitronScanner/HelitronScanner.jar pairends -hs head.helitronscanner.738.out -ts tail.helitronscanner.738.out  -o paired.helitrons.defaultsize
java -jar  HelitronScanner/HelitronScanner.jar draw -p paired.helitrons.defaultsize -g ../../CamTechGenomeComparison/18_mitochondria/swapMitoScaffolds/genome738sl.polished.mitoFixed.fa -o draw_helitrons_pure_default --pure
bioawk -c fastx '{print length($seq)}' draw_helitrons_pure_default.hel.fa |summary.sh
Total:  1,359,209
Count:  131
Mean:   10,375
Median: 10,668
Min:    262
Max:    19,815
This really only gets rid of 5 large helitrons. but knocks out 1MB worth of helitrons.
```
