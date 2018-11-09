# Fix redtandem duplications for scn genome 738

Some scaffolds are missing from the tandem duplication gff, found in circos plot.
```
/work/GIF/remkv6/Baum/CamTechGenomeComparison/40_tandemdups/split/01_rehashMissingScaffolds


#what should be the scaffold count?
awk '{print $1}' ../redtandem.out* |sort|uniq|wc
    669     669    2572
#so the count in the final renamed version should be 668 without the hash "#"

cat ../redtandem.out* |grep -v "#" |sort -k1,1nr |awk '{print "#01010"$1"0101#",$2,$3,$4,$5,$6,$7,$8}'  >redtandemAll.out
cp ../necc_Files/renamer.sh .
less renamer.sh |awk '{print $6,$8}' |less |sed 's/=/\t/g' |sed 's/"//g' |sed 's/>//g' |sed 's/,//g' |awk '{print "#01010"$2"0101#",$3 }' |sed "s/^/sed -i 's\t/g" |sed "s/$/\tg' redtandemAll.out/g" |awk '{print $1,$2,$3"/"$4"/"$5"/"$6,$7}' >renamer2.sh
sh renamer2.sh
cp ../../../58_Renamatorium/33_fixTandem/sorted.final.renamed.bed.renamer.sh  .

sed -i 's/sorted.final.renamed.bed/redtandemAll.out/g' sorted.final.renamed.bed.renamer.sh
#had to change the sorting order of renaming, as 1-5syntmers were collapsing the 11,12,22,55, etc syntmers.
sh sorted.final.renamed.bed.renamer.sh

 less redtandemAll.out |awk '{print $1,$6,$7}' | tr ".." " " |tr "," "\n" |sed '/^$/d' |awk -v I=0 '{if(NF>2) {I=$1" "$2;print $0} else {print I,$1,$2}}' |awk '{print $1,$3,$4,$2}' >redTandemAll.bed
less redTandemAll.bed |awk '{print $1}' |sort|uniq|sed 's/1scaff/scaff/g' |sed 's/2scaff/scaff/g' |sed 's/3scaff/scaff/g' |sed 's/4scaff/scaff/g' |sed 's/5scaff/scaff/g' >redTandemAllFixed.bed

 #are they all here?
 awk '{print $1}' redTandemAll.bed |sort|uniq|wc
    668     668    5905
#yep
```
# survey the duplications
```
How large are the duplications?
[remkv6@condo039 33_fixTandem]$ less redTandemAll.bed |awk '{print $3-$2}' |summary.sh
Total:  18,720,169
Count:  20,577
Mean:   909
Median: 134
Min:    30
Max:    31,662
How many are larger than 500bp?
[remkv6@condo039 33_fixTandem]$ less redTandemAll.bed |awk '{print $3-$2}' |awk '$1>500' |summary.sh
Total:  16,779,467
Count:  5,936
Mean:   2,826
Median: 1,614
Min:    502
Max:    31,662

What are the frequency and size of the consensus TD?
less redTandemAll.bed |awk '{print $3-$2}' |uniq|summary.sh
Total:  4,948,803
Count:  2,284
Mean:   2,166
Median: 979
Min:    30
Max:    31,662

What is the subread overlap stats on the TD larger than the mean (909 bp)?
bedtools intersect -wo -f .9 -a <(awk '$3-$2>909' redTandemAllOldNamesSort.bed |awk '{print $1,$2,$3}' |tr " " "\t") -b ../../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |awk '{print $1,$2,$3}' |sort|uniq|wc
    717    2151   13438

How many are there total?
less redTandemAll.bed |awk '{print $3-$2}' |awk '$1>909' |wc
   4410    4410   22097

What percent is that?
16.258%

How many tandem duplications are spanned by at least 90% of a pread?
bedtools intersect -wo -f .9 -a  <(awk '($3-$2)>909' ../../33_fixTandem/redTandemAllOldNamesSort.bed|tr " " "\t" ) -b preads2genome.gff|awk '{print $1,$2,$3}' |sort|uniq|wc
  4241   12723   87096

bedtools intersect -wo -f .9 -a  <(awk '($3-$2)>909' ../../33_fixTandem/redTandemAllOldNamesSort.bed|tr " " "\t" ) -b ../ccs/ccs2genome.gff|awk '{print $1,$2,$3}' |sort|uniq|wc
      362    1086    6473
```
### Gene overlap counts
```
#genes that tandem duplications overlap
bedtools intersect -wo -a <(tr " " "\t"< redTandemAll.bed) -b <(awk '$3=="gene"' ../1_genomeNgff/fixed.augustus.gff3 ) |awk '$7=="gene"' |cut -f 5- |cut -f 9 |sort|uniq|wc
   6730    6730  273377

#genes that tandem duplications do not overlap TD
bedtools intersect -v -wo -b <(tr " " "\t"< redTandemAll.bed) -a <(awk '$3=="gene"' ../1_genomeNgff/fixed.augustus.gff3 ) |awk '$3=="gene"' |cut -f 9 |sort|uniq|wc
  23039   23039  936046

  bedtools intersect -wo -b <(tr " " "\t"< redTandemAll.bed) -a <(less ../30_BenNewHGT/HighConfHGT.gene.gff|cut -f 1-9 ) |cut -f 9 |sort|uniq|wc
       39      39     819
```

# Gene Overlaps
```
bedtools intersect -wo -a <(tr " " "\t"< redTandemAll.bed) -b <(awk '$3=="gene"' ../1_genomeNgff/fixed.augustus.gff3 ) |awk '$7=="gene"' |cut -f 5- |cut -f 9 |sort|uniq|sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 >TandemDuplicationGeneList

cat SupportedIRFMergeClassified.list ../33_fixTandem/TandemDuplicationGeneList |sort|uniq -c |awk '$1==2' |wc
    969    1938   25194
cat NonRedundantLtrRetroelementMerge.list ../33_fixTandem/TandemDuplicationGeneList |sort|uniq -c |awk '$1==2' |wc
    656    1312   17056
cat AllPredictedEffectors.list ../33_fixTandem/TandemDuplicationGeneList |sort|uniq -c |awk '$1==2' |wc
    136     272    3536
cat <(sed 's/\./\t/2' ../30_BenNewHGT/HighConfHGT.list|cut -f 1) ../33_fixTandem/TandemDuplicationGeneList |sort|uniq -c |awk '$1==2' |wc
         38      76     988


```
