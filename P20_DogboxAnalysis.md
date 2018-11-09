# Dogbox containing genes comparisons to other genomic strata
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/25_DogBox
vi Dorsal-likeGenes.list
#copy paste Sebastians list

#how many known effectors have a dog box
cat ../20_Expression/121EffectorOnlyNewGeneNames.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
     60     120    1560

#how many genes overlap the 1265 line genes
cat ../20_Expression/Family1265.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
      0       0       0

#how many overlap family976 genes?
cat ../20_Expression/Family976.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
     21      42     546

#how many are in synteny?
cat ../20_Expression/GenesInSynteny.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
     36      72     936

#How many low confidence HGT have dog boxes?
cat ../20_Expression/HGT.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
     73     146    1898

#How many high confidence HGT have dog boxes?
 cat ../20_Expression/HGTNewNames.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
      2       4      52

#how many genes with dog boxes are found in ltr retroelements?
cat ../20_Expression/NonRedundantLtrRetroelementMerge.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
      3       6      78

#How many are found in DNA transposons
cat ../20_Expression/SupportedIRFMergeClassified.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
     15      30     390

#how many repeat affected genes have dog boxes?
cat ../../72_RepeatNTandemFinalGFF/RepeatAffectedGenesNoSimple.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
     41      82    1066


#how many dog box genes are in tandem duplications
cat ../20_Expression/tandem.gene.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
     24      48     624

#How many of the predicted effectors have a dog box?
cat ../20_Expression/PredictedEffector.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
     50     100    1300

#How many are secreted
cat ../20_Expression/Secretome.gene.list Dorsal-likeGenes.list |sort|uniq -c |awk '$1>1' |wc
    128     256    3328
```
