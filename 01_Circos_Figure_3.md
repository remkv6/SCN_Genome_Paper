# Creation of circos figure showing duplication overlap with effectors/HGT

```
Circos plot to represent TE TD HGT and effector predictions

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/31_Fig3Circos
#Copied template karyotype.conf, ideogram.conf, bands.conf, and ticks.conf from Globodera rostochiensis synteny study
cp ../../26_GloboderaSynteny/otherGloboderaGenomes/rostochiensis/circos/* .

#remade the karyttype file so it is only scn
bioawk -c fastx '{print $name,length($seq)}' ../1_genomeNgff/genome738sl.polished.mitoFixed.fa |sort -k2,2nr |head -n 20 |awk '{print "chr","-",$1,$1,"0",$2}' >Top20Scaf.kary


#need to make a histogram or plot points.  I am going with histogram first to get a feel for the presentation
awk '{print $2}'  Top20Scaf.kary |while read line; do echo "awk '\$1==\""$line"\" && \$3==\"gene\"' ../1_genomeNgff/fixed.augustus.gff3 >>Genesoftop20scaffolds.gff" ;done >GetGenesoftop20scaffolds.sh
Genesoftop20scaffolds.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >Genes.histo

awk '{print $2}'  Top20Scaf.kary |while read line; do echo "awk '\$1==\""$line"\" ' ../27_TandemRedo/ >>Genesoftop20scaffolds.gff" ;done >GetGenesoftop20scaffolds.sh

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/30_BenNewHGT
less HighConfHGT.list |sed 's/\.t/\t/g'|cut -f 1 |sort|uniq|grep  -f - <(awk '$3=="gene"' ../1_genomeNgff/fixed.augustus.gff3|sed 's/;/\t/g') >HighConfHGT.gene.gff


#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/31_Fig3Circos
 awk '{print $2}'  Top20Scaf.kary |while read line; do echo "awk '\$1==\""$line"\"' HighConfHGT.gene.gff >>HGTGenesoftop20scaffolds.gff" ;done >HGTGetGenesoftop20scaffolds.sh
 sh HGTGetGenesoftop20scaffolds.sh
less HGTGenesoftop20scaffolds.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >HGTGenes.histo


ln -s ../26_ExpressionSets/AllPredictedEffectors.list
less AllPredictedEffectors.list |grep  -f - <(awk '$3=="gene"' Genesoftop20scaffolds.gff|sed 's/;/\t/g') >AllPredictedEffectorsoftop20scaffolds.gff
less AllPredictedEffectorsoftop20scaffolds.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >AllPredictedEffectors.histo

ln -s ../23_LTR_finder/LtrRetroelementClassified.gff
awk '{print $2}'  Top20Scaf.kary |while read line; do echo "awk '\$1==\""$line"\"' LtrRetroelementClassified.gff >>LTRsoftop20scaffolds.gff" ;done >LTRsoftop20scaffolds.sh
less LTRsoftop20scaffolds.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >LTRsoftop20scaffolds.histo

ln -s ../24_IRF_DNATrans/SupportedIRFMergeClassified.gff
awk '{print $2}'  Top20Scaf.kary |while read line; do echo "awk '\$1==\""$line"\"' SupportedIRFMergeClassified.gff >>TIRsoftop20scaffolds.gff" ;done >TIRsoftop20scaffolds.sh
less TIRsoftop20scaffolds.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >TIRsoftop20scaffolds.histo

ln -s ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff
awk '{print $2}'  Top20Scaf.kary |while read line; do echo "awk '\$1==\""$line"\"' genome738sl.polished.mitoFixed.fa.out.gff >>Repeatsoftop20scaffolds.gff" ;done >Repeatsoftop20scaffolds.sh
less Repeatsoftop20scaffolds.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >Repeatsoftop20scaffolds.histo

ln -s ../27_TandemRedo/TotalSections.gff
awk '{print $2}'  Top20Scaf.kary |while read line; do echo "awk '\$1==\""$line"\"' TotalSections.gff >>TDRoftop20scaffolds.gff" ;done >TDRoftop20scaffolds.sh
less TDRoftop20scaffolds.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >TDRoftop20scaffolds.histo
```

### these needed changed to heatmap, and were not very interesting for the top 20 largest scaffolds.  decided to look at the scaffolds with the greatest numbers of duplications

```
awk '($5-$4)>2000 {print $1}' TotalSections.gff |sort|uniq -c |sort -k1,1nr |awk '{print $2}' |head |cdbyank ../1_genomeNgff/genome738sl.polished.mitoFixed.fa.cidx |bioawk -c fastx '{print "chr","-",$name,$name,"0",length($seq)}' >HighTDRTop10Scaffs.kary

less HighTDRTop10Scafs.kary|awk '{print $2}' |while read line; do echo "awk '\$1==\""$line"\" && \$3==\"gene\"' ../1_genomeNgff/fixed.augustus.gff3 >>HighTDRTopScaffs.gff";done >GetGenesofHighTDRTopScaffs.gff.sh
sh GetGenesofHighTDRTopScaffs.gff.sh

#reiteration of the first run to generate appropriate histo files
less HighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >Genes.histo
less HighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' |awk '{print $1,$2/100000,$3/100000}' |awk '{printf $1; printf "%16f", $2; printf "%17f\n", $3 }' |sed 's/\./\t/g' |awk '{print $1,$2*100000,($4*100000)+100000}' |sort|uniq -c |awk '{print $2,$3,$4,$1}'  >GenesEst.histo

awk '($5-$4)>2000 {print $1}' TotalSections.gff |sort|uniq -c |sort -k1,1nr |awk '{print $2}' |head |grep  -f - <(awk '$3=="gene"' HighTDRTopScaffs.gff|sed 's/;/\t/g') >AllPredictedEffectorsofHighTDRTopScaffs.gff
less HighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >AllPredictedEffectors.histo

 awk '($5-$4)>2000 {print $1}' TotalSections.gff |sort|uniq -c |sort -k1,1nr |awk '{print $2}' |head|while read line; do echo "awk '\$1==\""$line"\"' HighConfHGT.gene.gff >>HGTGenesofHighTDRTopScaffs.gff" ;done >HGTGetGenesofHighTDRTopScaffs.sh
 sh HGTGetGenesofHighTDRTopScaffs.sh
less HGTGenesofHighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >HGTGenesofHighTDRTopScaffs.histo

less AllPredictedEffectors.list |grep  -f - <(awk '$3=="gene"' HighTDRTopScaffs.gff|sed 's/;/\t/g') >AllPredictedEffectorsofHighTDRTopScaffs.gff
less AllPredictedEffectorsofHighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >AllPredictedEffectorsHighTDRTopScaffs.histo


#ln -s ../23_LTR_finder/LtrRetroelementClassified.gff
awk '($5-$4)>2000 {print $1}' TotalSections.gff |sort|uniq -c |sort -k1,1nr |awk '{print $2}' |head|while read line; do echo "awk '\$1==\""$line"\"' LtrRetroelementClassified.gff >>LTRsofHighTDRTopScaffs.gff" ;done >LTRsofHighTDRTopScaffs.sh
sh LTRsofHighTDRTopScaffs.sh
less LTRsofHighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >LTRsofHighTDRTopScaffs.histo

ln -s ../24_IRF_DNATrans/SupportedIRFMergeClassified.gff
awk '($5-$4)>2000 {print $1}' TotalSections.gff |sort|uniq -c |sort -k1,1nr |awk '{print $2}' |head |while read line; do echo "awk '\$1==\""$line"\"' SupportedIRFMergeClassified.gff >>TIRsofHighTDRTopScaffs.gff" ;done >TIRsofHighTDRTopScaffs.gff.sh
sh TIRsofHighTDRTopScaffs.gff.sh
less TIRsofHighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >TIRsHighTDRTopScaffs.histo

ln -s ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff
awk '($5-$4)>2000 {print $1}' TotalSections.gff |sort|uniq -c |sort -k1,1nr |awk '{print $2}' |head |while read line; do echo "awk '\$1==\""$line"\"' genome738sl.polished.mitoFixed.fa.out.gff >>RepeatsofHighTDRTopScaffs.gff" ;done >RepeatsHighTDRTopScaffs.sh
sh Repeatsoftop20scaffolds.sh
less RepeatsofHighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >RepeatsofHighTDRTopScaffs.histo
less RepeatsofHighTDRTopScaffs.histo |awk '{print $1,$2/100000,$3/100000}' |awk '{printf $1; printf "%16f", $2; printf "%17f\n", $3 }'|sed 's/\./\t/g' |awk '{print $1,$2*100000,($4*100000)+100000}' |sort |uniq -c |sort -k1,1nr |awk '{print $2,$3,$4,$1}' >RepeatEst.histo
ln -s ../27_TandemRedo/TotalSections.gff
awk '($5-$4)>2000 {print $1}' TotalSections.gff |sort|uniq -c |sort -k1,1nr |awk '{print $2}' |head |while read line; do echo "awk '\$1==\""$line"\"' TotalSections.gff >>TDRofHighTDRTopScaffs.gff" ;done >TDRofHighTDRTopScaffs.sh
sh TDRofHighTDRTopScaffs.sh
less TDRofHighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >TDRofHighTDRTopScaffs.histo
```


### Figured out what needed displayed here, effecotr/hgt overlap with duplicative tracks.  So grabbign top ten scaffolds from this list and repeating above analysis

```
cat ../26_ExpressionSets/AllPredictedEffectors.list ../26_ExpressionSets/HGTRevised.list |sort|uniq|cat  - <(cat ../26_ExpressionSets/SupportedIRFMergeClassified.list ../26_ExpressionSets/NonRedundantLtrRetroelementMerge.list ../26_ExpressionSets/tandem.gene.list |sort|uniq) |sort|uniq -c |awk '$1==2' |awk '{print $2}' |grep -w -f - <(sed 's/;/\t/g' ../1_genomeNgff/fixed.augustus.gff3 ) |awk '$3=="gene"' |awk '{print $1}' |sort|uniq -c |sort -k1,1nr |head |awk '{print $2}' >TopScaffoldsofEffHGToverlappingTIRLTRTDR.list


module load bioawk
module load cdbfasta
less TopScaffoldsofEffHGToverlappingTIRLTRTDR.list |cdbyank ../1_genomeNgff/genome738sl.polished.mitoFixed.fa.cidx |bioawk -c fastx '{print "chr","-",$name,$name,"0",length($seq),"green"}' >HighTDRTop10Scaffs.kary

less HighTDRTop10Scaffs.kary|awk '{print $3}' |while read line; do echo "awk '\$1==\""$line"\" && \$3==\"gene\"' ../1_genomeNgff/fixed.augustus.gff3 >>HighTDRTopScaffs.gff";done >GetGenesofHighTDRTopScaffs.gff.sh
sh GetGenesofHighTDRTopScaffs.gff.sh
```
### reiteration of the first run to generate appropriate histo files
```
less HighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >Genes.histo
less HighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' |awk '{print $1,$2/100000,$3/100000}' |awk '{printf $1; printf "%16f", $2; printf "%17f\n", $3 }' |sed 's/\./\t/g' |awk '{print $1,$2*100000,($4*100000)+100000}' |sort|uniq -c |awk '{print $2,$3,$4,$1}'  >GenesEst.histo

less TopScaffoldsofEffHGToverlappingTIRLTRTDR.list |grep  -f - <(awk '$3=="gene"' HighTDRTopScaffs.gff|sed 's/;/\t/g') >AllPredictedEffectorsofHighTDRTopScaffs.gff
less HighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >AllPredictedEffectors.histo

less TopScaffoldsofEffHGToverlappingTIRLTRTDR.list|while read line; do echo "awk '\$1==\""$line"\"' HighConfHGT.gene.gff >>HGTGenesofHighTDRTopScaffs.gff" ;done >HGTGetGenesofHighTDRTopScaffs.sh
 sh HGTGetGenesofHighTDRTopScaffs.sh
less HGTGenesofHighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >HGTGenesofHighTDRTopScaffs.histo

less AllPredictedEffectors.list |grep  -f - <(awk '$3=="gene"' HighTDRTopScaffs.gff|sed 's/;/\t/g') >AllPredictedEffectorsofHighTDRTopScaffs.gff
less AllPredictedEffectorsofHighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >AllPredictedEffectorsHighTDRTopScaffs.histo


#ln -s ../23_LTR_finder/LtrRetroelementClassified.gff
less TopScaffoldsofEffHGToverlappingTIRLTRTDR.list|while read line; do echo "awk '\$1==\""$line"\"' LtrRetroelementClassified.gff >>LTRsofHighTDRTopScaffs.gff" ;done >LTRsofHighTDRTopScaffs.sh
sh LTRsofHighTDRTopScaffs.sh
less LTRsofHighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >LTRsofHighTDRTopScaffs.histo
q
ln -s ../24_IRF_DNATrans/SupportedIRFMergeClassified.gff
less TopScaffoldsofEffHGToverlappingTIRLTRTDR.list |while read line; do echo "awk '\$1==\""$line"\"' SupportedIRFMergeClassified.gff >>TIRsofHighTDRTopScaffs.gff" ;done >TIRsofHighTDRTopScaffs.gff.sh
sh TIRsofHighTDRTopScaffs.gff.sh
less TIRsofHighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >TIRsHighTDRTopScaffs.histo

ln -s ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff
less TopScaffoldsofEffHGToverlappingTIRLTRTDR.list |while read line; do echo "awk '\$1==\""$line"\"' genome738sl.polished.mitoFixed.fa.out.gff >>RepeatsofHighTDRTopScaffs.gff" ;done >RepeatsHighTDRTopScaffs.sh
sh RepeatsHighTDRTopScaffs.sh
less RepeatsofHighTDRTopScaffs.gff|awk '{print $1"\t"$4"\t"$5"\t100"}' >RepeatsofHighTDRTopScaffs.histo
less RepeatsofHighTDRTopScaffs.histo |awk '{print $1,$2/100000,$3/100000}' |awk '{printf $1; printf "%16f", $2; printf "%17f\n", $3 }'|sed 's/\./\t/g' |awk '{print $1,$2*100000,($4*100000)+100000}' |sort |uniq -c |sort -k1,1nr |awk '{print $2,$3,$4,$1}' >RepeatEst.histo
ln -s ../27_TandemRedo/TotalSections.gff
less TopScaffoldsofEffHGToverlappingTIRLTRTDR.list |while read line; do echo "awk '\$1==\""$line"\"' redTandemAll.bed >>TDRofHighTDRTopScaffs.bed" ;done >TDRofHighTDRTopScaffs.sh
sh TDRofHighTDRTopScaffs.sh
cp TDRofHighTDRTopScaffs.bed TDRofHighTDRTopScaffs.histo
```

![circosPlot](assets/Figure3.jpg)
