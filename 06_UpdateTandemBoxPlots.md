# 400 extra genes were found in tandem duplicates, need to modify figures

```
/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/38_BoxPlotsRedoExp
for f in ../26_ExpressionSets/*Fold.out; do ln -s $f; done
ln -s ../33_fixTandem/TandemDuplicationGeneList
ln -s ../26_ExpressionSets/LogFoldChange


#All of these need to be recalculated
unlink tandem.gene.listFold.out
unlink AllPredictedEffectorTandem.listFold.out
unlink DorsalEffectorTandem.listFold.out
unlink KnownEffectorTandem.listFold.out
unlink PredictedEffectorTandem.listFold.out
```   

### Recalculate overlap and add fold
```
grep -w -f TandemDuplicationGeneList  LogFoldChange |cut -f 2 >TandemDuplicationGeneListFold.out &

cat ../26_ExpressionSets/AllPredictedEffectors.list TandemDuplicationGeneList |sort|uniq -c |awk '$1==2 {print $2}' >AllPredictedEffectorTandem.list
grep -w -f AllPredictedEffectorTandem.list  LogFoldChange |cut -f 2 >AllPredictedEffectorTandem.listFold.out &

cat ../26_ExpressionSets/Dorsal-likeGenesFixed.list TandemDuplicationGeneList |sort|uniq -c |awk '$1==2 {print $2}' >DorsalEffectorTandem.list
grep -w -f DorsalEffectorTandem.list  LogFoldChange |cut -f 2 >DorsalEffectorTandem.listFold.out &

cat ../26_ExpressionSets/121EffectorOnlyNewGeneNames.list TandemDuplicationGeneList |sort|uniq -c |awk '$1==2 {print $2}' >KnownEffectorTandem.list
grep -w -f KnownEffectorTandem.list  LogFoldChange |cut -f 2 > KnownEffectorTandem.listFold.out &

cat ../26_ExpressionSets/AllPredictedEffectors.list TandemDuplicationGeneList |sort|uniq -c |awk '$1==2 {print $2}' >PredictedEffectorTandem.list
grep -w -f PredictedEffectorTandem.list  LogFoldChange |cut -f 2 > PredictedEffectorTandem.listFold.out &

```

### Get these guys in a file for boxplots
```
awk '{print "Known Effectors\t" $0}' 121EffectorOnlyNewGeneNames.listFold.out         >>FoldChange4BoxPlot                
awk '{print "All Predicted Effectors\t" $0}' AllPredictedEffectorLTR.listFold.out             >>FoldChange4BoxPlot
awk '{print "All Predicted Effectors in DNA Transposons\t" $0}' AllPredictedEffectorTIR.listFold.out             >>FoldChange4BoxPlot
awk '{print "All Predicted Effectors in TDR\t" $0}' AllPredictedEffectorTandem.listFold.out          >>FoldChange4BoxPlot
awk '{print "All Predicted Effectors in LTR Retroelements\t" $0}' AllPredictedEffectors.listFold.out               >>FoldChange4BoxPlot
awk '{print "DOG-box Genes in LTR Retroelements\t" $0}' DorsalEffectorLTR.listFold.out                   >>FoldChange4BoxPlot
awk '{print "DOG-box Genes in DNA Transposons\t" $0}' DorsalEffectorTIR.listFold.out                   >>FoldChange4BoxPlot
awk '{print "DOG-box Genes in TDR\t" $0}' DorsalEffectorTandem.listFold.out                >>FoldChange4BoxPlot
awk '{print "DOG-box Genes\t" $0}' Dorsal-likeGenesFixed.listFold.out               >>FoldChange4BoxPlot
awk '{print "Genes with a Rnd_4-family976 Repeat\t" $0}' Family976.listFold.out                           >>FoldChange4BoxPlot
awk '{print "Genes with a Rnd_4-family1265 Repeat\t" $0}' Family1265.listFold.out                          >>FoldChange4BoxPlot
awk '{print "High Confidence HGT\t" $0}' HGTRevised.listFold.out                       >>FoldChange4BoxPlot
awk '{print "Known Effectors in LTR Retroelements\t" $0}' KnownEffectorLTR.listFold.out                    >>FoldChange4BoxPlot
awk '{print "Known Effectors in DNA Transposons\t" $0}' KnownEffectorTIR.listFold.out                    >>FoldChange4BoxPlot
awk '{print "Known Effectors in TDR\t" $0}' KnownEffectorTandem.listFold.out                 >>FoldChange4BoxPlot
awk '{print "Genes in LTR Retroelements\t" $0}' NonRedundantLtrRetroelementMerge.listFold.out    >>FoldChange4BoxPlot
awk '{print "Predicted Effectors\t" $0}' PredictedEffector.listFold.out                   >>FoldChange4BoxPlot
awk '{print "Predicted Effectors in LTR Retroelements\t" $0}' PredictedEffectorLTR.listFold.out                >>FoldChange4BoxPlot
awk '{print "Predicted Effectors in DNA Transposons\t" $0}' PredictedEffectorTIR.listFold.out                >>FoldChange4BoxPlot
awk '{print "Predicted Effectors in TDR\t" $0}' PredictedEffectorTandem.listFold.out             >>FoldChange4BoxPlot
awk '{print "Genes with Repeat-overlapping Exons\t" $0}' RepeatAffectedGenesNoSimple.listFold.out         >>FoldChange4BoxPlot
awk '{print "Secreted Genes\t" $0}' Secretome.gene.listFold.out                      >>FoldChange4BoxPlot
awk '{print "Genes in DNA Transposons\t" $0}' SupportedIRFMergeClassified.listFold.out         >>FoldChange4BoxPlot
awk '{print "Genes in TDR\t" $0}' TandemDuplicationGeneListFold.out                        >>FoldChange4BoxPlot
awk '{print "All Genes\t" $2}' LogFoldChange >>FoldChange4BoxPlot
```
### Make the plot in R
```
library(ggplot2)
FoldChange4BoxPlot <- read.table("FoldChange4BoxPlot" ,header=FALSE, sep="\t")

ggplot(data=FoldChange4BoxPlot, aes(y=V2,x =V1))+geom_violin(size=.5,fill = "#FF6666") +labs(title="Gene Expression of Genomic Strata",y="log(Fold Change)",x="") +coord_flip(ylim=c(-10,10)) + xlim( "Genes with Repeat-overlapping Exons", "All Predicted Effectors in DNA Transposons", "DOG-box Genes in DNA Transposons", "Predicted Effectors in DNA Transposons", "Known Effectors in DNA Transposons", "Genes in DNA Transposons", "All Predicted Effectors in LTR Retroelements", "DOG-box Genes in LTR Retroelements", "Predicted Effectors in LTR Retroelements", "Known Effectors in LTR Retroelements", "Genes in LTR Retroelements", "All Predicted Effectors in TDR", "DOG-box Genes in TDR", "Predicted Effectors in TDR", "Known Effectors in TDR", "Genes in TDR", "High Confidence HGT", "All Predicted Effectors", "DOG-box Genes",  "Predicted Effectors", "Known Effectors", "All Genes") + scale_colour_manual(values = cols) +  theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size = 10),title = element_text(size = 10),axis.ticks = element_line(colour = "black", size = 1),text = element_text( size = 1))
```

#  Do the same for the snp density plots
```
/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/39_BoxPlotsRedoSNPDensity
for f in ../28_SNPDensity/*SNP.out; do ln -s $f; done

unlink AllPredictedEffectorTandem.listSNP.out
unlink DorsalEffectorTandem.listSNP.out
unlink KnownEffectorTandem.listSNP.out
unlink PredictedEffectorTandem.listSNP.out
unlink tandem.gene.listSNP.out

ln -s ../38_BoxPlotsRedoExp/AllPredictedEffectorTandem.list
ln -s ../38_BoxPlotsRedoExp/DorsalEffectorTandem.list
ln -s ../38_BoxPlotsRedoExp/KnownEffectorTandem.list
ln -s ../38_BoxPlotsRedoExp/PredictedEffectorTandem.list
ln -s ../38_BoxPlotsRedoExp/TandemDuplicationGeneList

ln -s ../28_SNPDensity/SNPDensity
```

### get snp density for each of these GENES
```
grep -w -f TandemDuplicationGeneList  SNPDensity |cut -f 2 >TandemDuplicationGeneListSNP.out &
grep -w -f AllPredictedEffectorTandem.list  SNPDensity |cut -f 2 >AllPredictedEffectorTandem.listSNP.out &
grep -w -f DorsalEffectorTandem.list  SNPDensity |cut -f 2 >DorsalEffectorTandem.listSNP.out &
grep -w -f KnownEffectorTandem.list  SNPDensity |cut -f 2 >KnownEffectorTandem.listSNP.out &
grep -w -f PredictedEffectorTandem.list  SNPDensity |cut -f 2 >PredictedEffectorTandem.listSNP.out &
```

### get them into a box plot format for R
```
awk '{print "Known Effectors\t" log($4log(2))}' 121EffectorOnlyNewGeneNames.listSNP.out         >>SNPDensity4BoxPlot                
awk '{print "All Predicted Effectors\t" log($4log(2))}' AllPredictedEffectorLTR.listSNP.out             >>SNPDensity4BoxPlot
awk '{print "All Predicted Effectors in DNA Transposons\t" log($4log(2))}' AllPredictedEffectorTIR.listSNP.out             >>SNPDensity4BoxPlot
awk '{print "All Predicted Effectors in TDR\t" log($4log(2))}' AllPredictedEffectorTandem.listSNP.out          >>SNPDensity4BoxPlot
awk '{print "All Predicted Effectors in LTR Retroelements\t" log($4log(2))}' AllPredictedEffectors.listSNP.out               >>SNPDensity4BoxPlot
awk '{print "DOG-box Genes in LTR Retroelements\t" log($4log(2))}' DorsalEffectorLTR.listSNP.out                   >>SNPDensity4BoxPlot
awk '{print "DOG-box Genes in DNA Transposons\t" log($4log(2))}' DorsalEffectorTIR.listSNP.out                   >>SNPDensity4BoxPlot
awk '{print "DOG-box Genes in TDR\t" log($4log(2))}' DorsalEffectorTandem.listSNP.out                >>SNPDensity4BoxPlot
awk '{print "DOG-box Genes\t" log($4log(2))}' Dorsal-likeGenesFixed.listSNP.out               >>SNPDensity4BoxPlot
awk '{print "Genes with a Rnd_4-family976 Repeat\t" log($4log(2))}' Family976.listSNP.out                           >>SNPDensity4BoxPlot
awk '{print "Genes with a Rnd_4-family1265 Repeat\t" log($4log(2))}' Family1265.listSNP.out                          >>SNPDensity4BoxPlot
awk '{print "High Confidence HGT\t" log($4log(2))}' HGTRevised.listSNP.out                         >>SNPDensity4BoxPlot
awk '{print "Known Effectors in LTR Retroelements\t" log($4log(2))}' KnownEffectorLTR.listSNP.out                    >>SNPDensity4BoxPlot
awk '{print "Known Effectors in DNA Transposons\t" log($4log(2))}' KnownEffectorTIR.listSNP.out                    >>SNPDensity4BoxPlot
awk '{print "Known Effectors in TDR\t" log($4log(2))}' KnownEffectorTandem.listSNP.out                 >>SNPDensity4BoxPlot
awk '{print "Genes in LTR Retroelements\t" log($4log(2))}' NonRedundantLtrRetroelementMerge.listSNP.out    >>SNPDensity4BoxPlot
awk '{print "Predicted Effectors\t" log($4log(2))}' PredictedEffector.listSNP.out                   >>SNPDensity4BoxPlot
awk '{print "Predicted Effectors in LTR Retroelements\t" log($4log(2))}' PredictedEffectorLTR.listSNP.out                >>SNPDensity4BoxPlot
awk '{print "Predicted Effectors in DNA Transposons\t" log($4log(2))}' PredictedEffectorTIR.listSNP.out                >>SNPDensity4BoxPlot
awk '{print "Predicted Effectors in TDR\t" log($4log(2))}' PredictedEffectorTandem.listSNP.out             >>SNPDensity4BoxPlot
awk '{print "Genes with Repeat-overlapping Exons\t" log($4log(2))}' RepeatAffectedGenesNoSimple.listSNP.out         >>SNPDensity4BoxPlot
awk '{print "Secreted Genes\t" log($4log(2))}' Secretome.gene.listSNP.out                      >>SNPDensity4BoxPlot
awk '{print "Genes in DNA Transposons\t" log($4log(2))}' SupportedIRFMergeClassified.listSNP.out         >>SNPDensity4BoxPlot
awk '{print "Genes in TDR\t" log($4log(2))}' TandemDuplicationGeneListSNP.out                         >>SNPDensity4BoxPlot
awk '{print "All Genes\t" log($4log(2))}' SNPDensity >>SNPDensity4BoxPlot
```

#### The GG plot
```
ggplot(data=SNPDensity4BoxPlot, aes(y=V2,x =V1))+geom_violin(size=.5,fill = "#FF6666") +labs(title="SNP Density of Genomic Strata",y="log(SNP/CDS Length)",x="") +coord_flip(ylim=c(-10,0)) + xlim("Genes with Repeat-overlapping Exons", "All Predicted Effectors in DNA Transposons", "DOG-box Genes in DNA Transposons", "Predicted Effectors in DNA Transposons", "Known Effectors in DNA Transposons", "Genes in DNA Transposons","All Predicted Effectors in LTR Retroelements", "DOG-box Genes in LTR Retroelements", "Predicted Effectors in LTR Retroelements", "Known Effectors in LTR Retroelements", "Genes in LTR Retroelements", "All Predicted Effectors in TDR", "DOG-box Genes in TDR", "Predicted Effectors in TDR", "Known Effectors in TDR", "Genes in TDR", "High Confidence HGT", "All Predicted Effectors", "DOG-box Genes",  "Predicted Effectors", "Known Effectors", "All Genes") + scale_colour_manual(values = cols) +  theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size = 10),title = element_text(size = 10),axis.ticks = element_line(colour = "black", size = 1),text = element_text( size = 1))
```

# Calculate significance with gene overlap package
```
#grabbing important gene lists
for f in ../39_BoxPlotsRedoSNPDensity/*list; do ln -s $f; done
for f in ../39_BoxPlotsRedoSNPDensity/*List; do ln -s $f; done
for f in ../26_ExpressionSets/*list; do ln -s $f; done
 for f in ../28_SNPDensity/*list; do ln -s $f; done
ln -s ../38_BoxPlotsRedoExp/TandemDuplicationGeneList TandemDuplicationGene.list

#these are no longer of relevance
unlink tandem.gene.list
unlink SNPDensityAll10Percentile.list
unlink HighSNPDensity.list
unlink HighSNPDensity2Dev.list
unlink HighSNPDensity1Dev.list
unlink HighSNPDensity10Perc.list
unlink HGT.list
unlink HGTNewNames.list

```

### import all of the tables
```
for f in *; do echo $f" <- read.table(\""$f"\")" ;done |less
#make sure to fix this dirty bastard, variables cant start with numbers in R


x121EffectorOnlyNewGeneNames.list <- read.table("121EffectorOnlyNewGeneNames.list")
AllGenes.list <- read.table("AllGenes.list")
AllHGTRevised.list <- read.table("AllHGTRevised.list")
AllPredictedEffectorLTR.list <- read.table("AllPredictedEffectorLTR.list")
AllPredictedEffectors.list <- read.table("AllPredictedEffectors.list")
AllPredictedEffectorTandem.list <- read.table("AllPredictedEffectorTandem.list")
AllPredictedEffectorTIR.list <- read.table("AllPredictedEffectorTIR.list")
BuscoDNA.list <- read.table("BuscoDNA.list")
BuscoLTR.list <- read.table("BuscoLTR.list")
BuscoTandem.list <- read.table("BuscoTandem.list")
DorsalEffectorLTR.list <- read.table("DorsalEffectorLTR.list")
DorsalEffectorTandem.list <- read.table("DorsalEffectorTandem.list")
DorsalEffectorTIR.list <- read.table("DorsalEffectorTIR.list")
DorsallikeGenesFixed.list <- read.table("Dorsal-likeGenesFixed.list")
DuplicatedBuscosGene.list <- read.table("DuplicatedBuscosGene.list")
Family1265.list <- read.table("Family1265.list")
Family976.list <- read.table("Family976.list")
HGTRevised.list <- read.table("HGTRevised.list")
HighSignificant.list <- read.table("HighSignificant.list")
KnownEffectorLTR.list <- read.table("KnownEffectorLTR.list")
KnownEffectorTandem.list <- read.table("KnownEffectorTandem.list")
KnownEffectorTIR.list <- read.table("KnownEffectorTIR.list")
LowSignificant.list <- read.table("LowSignificant.list")
NonRedundantLtrRetroelementMerge.list <- read.table("NonRedundantLtrRetroelementMerge.list")
NotAllPredictedEffectorLTR.list <- read.table("NotAllPredictedEffectorLTR.list")
NotAllPredictedEffectorTandem.list <- read.table("NotAllPredictedEffectorTandem.list")
NotAllPredictedEffectorTIR.list <- read.table("NotAllPredictedEffectorTIR.list")
NotBuscoDNA.list <- read.table("NotBuscoDNA.list")
NotBuscoLTR.list <- read.table("NotBuscoLTR.list")
NotBuscoTandem.list <- read.table("NotBuscoTandem.list")
NotDorsalEffectorLTR.list <- read.table("NotDorsalEffectorLTR.list")
NotDorsalEffectorTandem.list <- read.table("NotDorsalEffectorTandem.list")
NotDorsalEffectorTIR.list <- read.table("NotDorsalEffectorTIR.list")
NotKnownEffectorLTR.list <- read.table("NotKnownEffectorLTR.list")
NotKnownEffectorTandem.list <- read.table("NotKnownEffectorTandem.list")
NotKnownEffectorTIR.list <- read.table("NotKnownEffectorTIR.list")
NotPredictedEffectorLTR.list <- read.table("NotPredictedEffectorLTR.list")
NotPredictedEffectorTandem.list <- read.table("NotPredictedEffectorTandem.list")
NotPredictedEffectorTIR.list <- read.table("NotPredictedEffectorTIR.list")
PredictedEffector.list <- read.table("PredictedEffector.list")
PredictedEffectorLTR.list <- read.table("PredictedEffectorLTR.list")
PredictedEffectorTandem.list <- read.table("PredictedEffectorTandem.list")
PredictedEffectorTIR.list <- read.table("PredictedEffectorTIR.list")
RepeatAffectedGenesNoSimple.list <- read.table("RepeatAffectedGenesNoSimple.list")
Secretome.gene.list <- read.table("Secretome.gene.list")
SNPDensityAll0SNPS.list <- read.table("SNPDensityAll0SNPS.list")
SNPDensityAll10percAfterRemoveZeros.list <- read.table("SNPDensityAll10percAfterRemoveZeros.list")
SNPDensityAll90Percentile.list <- read.table("SNPDensityAll90Percentile.list")
SupportedIRFMergeClassified.list <- read.table("SupportedIRFMergeClassified.list")
TandemDuplicationGene.list <- read.table("TandemDuplicationGene.list")
```

### Prepare comparison with High expression
```
ls -1 *list |awk '{print "go.obj <- newGeneOverlap("$1"$V1, HighSignificant.list$V1, genome.size=AllGenes.list)"}' |awk '{print $0"\ngo.obj"NR" <- testGeneOverlap(go.obj)"}' |less

go.obj <- newGeneOverlap(x121EffectorOnlyNewGeneNames.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj1 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllGenes.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj2 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllHGTRevised.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj3 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj4 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj5 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj6 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectors.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj7 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoDNA.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj8 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoLTR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj9 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoTandem.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj10 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj11 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj12 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj13 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsallikeGenesFixed.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj14 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DuplicatedBuscosGene.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj15 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family976.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj16 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family1265.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj17 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTRevised.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj18 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSignificant.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj19 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj20 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj21 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj22 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(LowSignificant.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj23 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NonRedundantLtrRetroelementMerge.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj24 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj25 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj26 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj27 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoDNA.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj28 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoLTR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj29 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoTandem.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj30 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj31 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj32 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj33 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj34 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj35 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj36 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj37 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj38 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj39 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffector.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj40 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj41 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj42 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj43 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(RepeatAffectedGenesNoSimple.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj44 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll0SNPS.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj45 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll10percAfterRemoveZeros.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj46 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll90Percentile.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj47 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Secretome.gene.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj48 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SupportedIRFMergeClassified.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj49 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(TandemDuplicationGene.list$V1, HighSignificant.list$V1, genome.size=AllGenes.list)
go.obj50 <- testGeneOverlap(go.obj)


AllCombine <- c(go.obj1, go.obj2, go.obj3, go.obj4, go.obj5, go.obj6, go.obj7, go.obj8, go.obj9, go.obj10, go.obj11, go.obj12, go.obj13, go.obj14, go.obj15, go.obj16, go.obj17, go.obj18, go.obj19, go.obj20, go.obj21, go.obj22, go.obj23, go.obj24, go.obj25, go.obj26, go.obj27, go.obj28, go.obj29, go.obj30, go.obj31, go.obj32, go.obj33, go.obj34, go.obj35, go.obj36, go.obj37, go.obj38, go.obj39, go.obj40, go.obj41, go.obj42, go.obj43, go.obj44, go.obj45, go.obj46, go.obj47, go.obj48, go.obj49, go.obj50)

library(RJSONIO)
 exportJSON <- toJSON(AllCombine)
 write(exportJSON,"HighSignificant.out")
#exit R
 ls -1 *list |paste - <(less HighSignificant.out|sed 's/\[/\t/g' |awk 'NF<15' |grep "\"" | sed '0~5 s/$/\nwait/g' |tr "\n" "\t" |sed 's/wait/\n/g' |sed 's/,/\t/g' |sed 's/"//g' |cut -f 1-7,9- |sed 's/}//g') |less -S

```

Result of above
```
121EffectorOnlyNewGeneNames.list           notA: 28480            inA: 80          notA: 1170             inA: 41          ]    pval: 7.3351e-27        Jaccard: 0.031758       is.teste
AllGenes.list              notA: 28559            inA: 1           notA: 1211             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
AllHGTRevised.list                 notA: 27096            inA: 1464        notA: 1117             inA: 94       odds.ratio: 1.5575      pval: 8.9494e-05        Jaccard: 0.03514
AllPredictedEffectorLTR.list               notA: 28540            inA: 20          notA: 1209             inA: 2        odds.ratio: 2.3605      pval: 0.22477   Jaccard: 0.0016247
AllPredictedEffectorTIR.list               notA: 28510            inA: 50          notA: 1203             inA: 8        odds.ratio:  3.792      pval: 0.0023057 Jaccard: 0.0063442
AllPredictedEffectorTandem.list            notA: 28448            inA: 112         notA: 1187             inA: 24       odds.ratio:  5.135      pval: 1.2833e-09        Jaccard: 0.01814
AllPredictedEffectors.list                 notA: 28232            inA: 328         notA: 1108             inA: 103      odds.ratio:  8.001      pval: 1.2287e-49        Jaccard: 0.06692
BuscoDNA.list              notA: 28536            inA: 24          notA: 1211             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
BuscoLTR.list              notA: 28540            inA: 20          notA: 1211             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
BuscoTandem.list                   notA: 28483            inA: 77          notA: 1209             inA: 2        odds.ratio: 0.61193     pval: 0.83681   Jaccard: 0.0015528      is.teste
DorsalEffectorLTR.list             notA: 28555            inA: 5           notA: 1210             inA: 1        odds.ratio: 4.7193      pval: 0.22057   Jaccard: 0.00082237     is.teste
DorsalEffectorTIR.list             notA: 28552            inA: 8           notA: 1204             inA: 7        odds.ratio: 20.743      pval: 8.746e-07 Jaccard: 0.0057424      is.teste
DorsalEffectorTandem.list                  notA: 28537            inA: 23          notA: 1193             inA: 18       odds.ratio: 18.714      pval: 6.8316e-15        Jaccard: 0.01458
Dorsal-likeGenesFixed.list                 notA: 28471            inA: 89          notA: 1140             inA: 71       odds.ratio: 19.917      pval: 2.6556e-55        Jaccard: 0.05461
DuplicatedBuscosGene.list                  notA: 28232            inA: 328         notA: 1196             inA: 15       odds.ratio: 1.0795      pval: 0.42433   Jaccard: 0.0097466
Family976.list             notA: 28473            inA: 87          notA: 1197             inA: 14       odds.ratio: 3.8278      pval: 5.9472e-05        Jaccard: 0.010786       is.teste
Family1265.list            notA: 28451            inA: 109         notA: 1208             inA: 3        odds.ratio: 0.64823     pval: 0.83887   Jaccard: 0.0022727      is.tested: true
HGTRevised.list            notA: 28424            inA: 136         notA: 1196             inA: 15       odds.ratio: 2.6212      pval: 0.0013058 Jaccard: 0.011136       is.tested: true
HighSignificant.list               notA: 28560            inA: 0           notA: 0                inA: 1211     odds.ratio:   Infinity  pval:      0    Jaccard:      1 is.tested: true
KnownEffectorLTR.list              notA: 28555            inA: 5           notA: 1211             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
KnownEffectorTIR.list              notA: 28548            inA: 12          notA: 1206             inA: 5        odds.ratio: 9.8611      pval: 0.00045379        Jaccard: 0.0040883
KnownEffectorTandem.list                   notA: 28528            inA: 32          notA: 1202             inA: 9        odds.ratio: 6.6741      pval: 3.1968e-05        Jaccard: 0.00724
LowSignificant.list                notA: 27992            inA: 568         notA: 1211             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
NonRedundantLtrRetroelementMerge.list              notA: 27170            inA: 1390        notA: 1200             inA: 11       odds.ratio: 0.17918     pval:      1    Jaccard: 0.00422
NotAllPredictedEffectorLTR.list            notA: 27190            inA: 1370        notA: 1202             inA: 9        odds.ratio: 0.14861     pval:      1    Jaccard: 0.003487
NotAllPredictedEffectorTIR.list            notA: 26732            inA: 1828        notA: 1182             inA: 29       odds.ratio: 0.35879     pval:      1    Jaccard: 0.0095426
NotAllPredictedEffectorTandem.list                 notA: 22487            inA: 6073        notA: 1063             inA: 148      odds.ratio: 0.51556     pval:      1    Jaccard: 0.02031
NotBuscoDNA.list                   notA: 26706            inA: 1854        notA: 1174             inA: 37       odds.ratio: 0.45398     pval:      1    Jaccard: 0.012072       is.teste
NotBuscoLTR.list                   notA: 27190            inA: 1370        notA: 1200             inA: 11       odds.ratio: 0.18193     pval:      1    Jaccard: 0.0042619      is.teste
NotBuscoTandem.list                notA: 22472            inA: 6088        notA: 1041             inA: 170      odds.ratio: 0.60281     pval:      1    Jaccard: 0.023291       is.teste
NotDorsalEffectorLTR.list                  notA: 27175            inA: 1385        notA: 1201             inA: 10       odds.ratio: 0.16338     pval:      1    Jaccard: 0.0038521
NotDorsalEffectorTIR.list                  notA: 26690            inA: 1870        notA: 1181             inA: 30       odds.ratio: 0.36257     pval:      1    Jaccard: 0.0097371
NotDorsalEffectorTandem.list               notA: 22415            inA: 6145        notA: 1057             inA: 154      odds.ratio: 0.53148     pval:      1    Jaccard: 0.020935
NotKnownEffectorLTR.list                   notA: 27175            inA: 1385        notA: 1200             inA: 11       odds.ratio: 0.17986     pval:      1    Jaccard: 0.0042373
NotKnownEffectorTIR.list                   notA: 26694            inA: 1866        notA: 1179             inA: 32       odds.ratio: 0.38828     pval:      1    Jaccard: 0.0104 is.teste
NotKnownEffectorTandem.list                notA: 22424            inA: 6136        notA: 1048             inA: 163      odds.ratio: 0.56843     pval:      1    Jaccard: 0.022186
NotPredictedEffectorLTR.list               notA: 27187            inA: 1373        notA: 1202             inA: 9        odds.ratio: 0.14827     pval:      1    Jaccard: 0.003483
NotPredictedEffectorTIR.list               notA: 26728            inA: 1832        notA: 1180             inA: 31       odds.ratio: 0.38329     pval:      1    Jaccard: 0.010187
NotPredictedEffectorTandem.list            notA: 22471            inA: 6089        notA: 1053             inA: 158      odds.ratio: 0.55377     pval:      1    Jaccard: 0.021644
PredictedEffector.list             notA: 28322            inA: 238         notA: 1157             inA: 54       odds.ratio: 5.5533      pval: 7.2954e-21        Jaccard: 0.037267
PredictedEffectorLTR.list                  notA: 28543            inA: 17          notA: 1209             inA: 2        odds.ratio: 2.7773      pval: 0.17969   Jaccard: 0.0016287
PredictedEffectorTIR.list                  notA: 28514            inA: 46          notA: 1205             inA: 6        odds.ratio: 3.0863      pval: 0.018495  Jaccard: 0.0047733
PredictedEffectorTandem.list               notA: 28448            inA: 112         notA: 1187             inA: 24       odds.ratio:  5.135      pval: 1.2833e-09        Jaccard: 0.01814
RepeatAffectedGenesNoSimple.list                   notA: 16503            inA: 12057       notA: 911              inA: 300      odds.ratio: 0.45076     pval:      1    Jaccard: 0.02261
SNPDensityAll0SNPS.list            notA: 23974            inA: 4586        notA: 1185             inA: 26       odds.ratio: 0.1147      pval:      1    Jaccard: 0.0044851      is.teste
SNPDensityAll10percAfterRemoveZeros.list                   notA: 26133            inA: 2427        notA: 1122             inA: 89       odds.ratio: 0.8541      pval: 0.93013   Jaccard:
SNPDensityAll90Percentile.list             notA: 25680            inA: 2880        notA: 1122             inA: 89       odds.ratio: 0.7073      pval: 0.99949   Jaccard: 0.021755
Secretome.gene.list                notA: 26685            inA: 1875        notA: 864              inA: 347      odds.ratio: 5.7144      pval: 6.9306e-115       Jaccard: 0.11244
SupportedIRFMergeClassified.list                   notA: 26682            inA: 1878        notA: 1174             inA: 37       odds.ratio: 0.44778     pval:      1    Jaccard: 0.01197
TandemDuplicationGene.list                 notA: 22013            inA: 6547        notA: 1028             inA: 183      odds.ratio: 0.59858     pval:      1    Jaccard: 0.023589

```

### Prepare comparison with Low expression

```
ls -1 *list |awk '{print "go.obj <- newGeneOverlap("$1"$V1, LowSignificant.list$V1, genome.size=AllGenes.list)"}' |awk '{print $0"\ngo.obj"NR" <- testGeneOverlap(go.obj)"}' |less

go.obj <- newGeneOverlap(x121EffectorOnlyNewGeneNames.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj1 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllGenes.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj2 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllHGTRevised.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj3 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj4 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj5 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj6 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectors.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj7 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoDNA.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj8 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoLTR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj9 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoTandem.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj10 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj11 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj12 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj13 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsallikeGenesFixed.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj14 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DuplicatedBuscosGene.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj15 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family976.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj16 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family1265.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj17 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTRevised.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj18 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSignificant.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj19 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj20 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj21 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj22 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(LowSignificant.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj23 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NonRedundantLtrRetroelementMerge.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj24 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj25 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj26 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj27 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoDNA.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj28 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoLTR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj29 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoTandem.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj30 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj31 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj32 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj33 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj34 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj35 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj36 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj37 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj38 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj39 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffector.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj40 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj41 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj42 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj43 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(RepeatAffectedGenesNoSimple.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj44 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll0SNPS.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj45 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll10percAfterRemoveZeros.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj46 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll90Percentile.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj47 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Secretome.gene.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj48 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SupportedIRFMergeClassified.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj49 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(TandemDuplicationGene.list$V1, LowSignificant.list$V1, genome.size=AllGenes.list)
go.obj50 <- testGeneOverlap(go.obj)


AllCombine <- c(go.obj1, go.obj2, go.obj3, go.obj4, go.obj5, go.obj6, go.obj7, go.obj8, go.obj9, go.obj10, go.obj11, go.obj12, go.obj13, go.obj14, go.obj15, go.obj16, go.obj17, go.obj18, go.obj19, go.obj20, go.obj21, go.obj22, go.obj23, go.obj24, go.obj25, go.obj26, go.obj27, go.obj28, go.obj29, go.obj30, go.obj31, go.obj32, go.obj33, go.obj34, go.obj35, go.obj36, go.obj37, go.obj38, go.obj39, go.obj40, go.obj41, go.obj42, go.obj43, go.obj44, go.obj45, go.obj46, go.obj47, go.obj48, go.obj49, go.obj50, go.obj51, go.obj52, go.obj53)

library(RJSONIO)
 exportJSON <- toJSON(AllCombine)
 write(exportJSON,"LowSignificant.out")

#format the file so it is a readable table
ls -1 *list |paste - <(less LowSignificant.out|sed 's/\[/\t/g' |awk 'NF<15' |grep "\"" | sed '0~5 s/$/\nwait/g' |tr "\n" "\t" |sed 's/wait/\n/g' |sed 's/,/\t/g' |sed 's/"//g' |cut -f 1-7,9- |sed 's/}//g') |less -S
```
### Results from Low significant expression
```
121EffectorOnlyNewGeneNames.list           notA: 29093            inA: 110         notA: 557              inA: 11          ]    pval: 2.1312e-05        Jaccard: 0.016224       is.teste
AllGenes.list              notA: 29202            inA: 1           notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
AllHGTRevised.list                 notA: 27678            inA: 1525        notA: 535              inA: 33       odds.ratio: 1.1195      pval: 0.29201   Jaccard: 0.015767       is.teste
AllPredictedEffectorLTR.list               notA: 29181            inA: 22          notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.teste
AllPredictedEffectorTIR.list               notA: 29147            inA: 56          notA: 566              inA: 2        odds.ratio: 1.8392      pval: 0.30378   Jaccard: 0.0032051
AllPredictedEffectorTandem.list            notA: 29075            inA: 128         notA: 560              inA: 8        odds.ratio: 3.2447      pval: 0.0046783 Jaccard: 0.011494
AllPredictedEffectors.list                 notA: 28793            inA: 410         notA: 547              inA: 21       odds.ratio: 2.6958      pval: 9.7619e-05        Jaccard: 0.02147
BuscoDNA.list              notA: 29179            inA: 24          notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
BuscoLTR.list              notA: 29183            inA: 20          notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
BuscoTandem.list                   notA: 29126            inA: 77          notA: 566              inA: 2        odds.ratio: 1.3366      pval: 0.44644   Jaccard: 0.0031008      is.teste
DorsalEffectorLTR.list             notA: 29197            inA: 6           notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
DorsalEffectorTIR.list             notA: 29190            inA: 13          notA: 566              inA: 2        odds.ratio: 7.9331      pval: 0.032376  Jaccard: 0.0034423      is.teste
DorsalEffectorTandem.list                  notA: 29167            inA: 36          notA: 563              inA: 5        odds.ratio: 7.1928      pval: 0.0010557 Jaccard: 0.0082781
Dorsal-likeGenesFixed.list                 notA: 29056            inA: 147         notA: 555              inA: 13       odds.ratio: 4.6295      pval: 1.3198e-05        Jaccard: 0.01818
DuplicatedBuscosGene.list                  notA: 28868            inA: 335         notA: 560              inA: 8        odds.ratio:  1.231      pval: 0.33264   Jaccard: 0.0088594
Family976.list             notA: 29105            inA: 98          notA: 565              inA: 3        odds.ratio: 1.5769      pval: 0.30332   Jaccard: 0.0045045      is.tested: true
Family1265.list            notA: 29091            inA: 112         notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
HGTRevised.list            notA: 29061            inA: 142         notA: 559              inA: 9        odds.ratio: 3.2947      pval: 0.002538  Jaccard: 0.012676       is.tested: true
HighSignificant.list               notA: 27992            inA: 1211        notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
KnownEffectorLTR.list              notA: 29198            inA: 5           notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
KnownEffectorTIR.list              notA: 29187            inA: 16          notA: 567              inA: 1        odds.ratio:  3.217      pval: 0.27932   Jaccard: 0.0017123      is.teste
KnownEffectorTandem.list                   notA: 29167            inA: 36          notA: 563              inA: 5        odds.ratio: 7.1928      pval: 0.0010557 Jaccard: 0.0082781
LowSignificant.list                notA: 29203            inA: 0           notA: 0                inA: 568      odds.ratio:   Infinity  pval:      0    Jaccard:      1 is.tested: true
NonRedundantLtrRetroelementMerge.list              notA: 27809            inA: 1394        notA: 561              inA: 7        odds.ratio: 0.24893     pval:      1    Jaccard: 0.00356
NotAllPredictedEffectorLTR.list            notA: 27831            inA: 1372        notA: 561              inA: 7        odds.ratio: 0.25312     pval:      1    Jaccard: 0.0036082
NotAllPredictedEffectorTIR.list            notA: 27362            inA: 1841        notA: 552              inA: 16       odds.ratio: 0.43081     pval: 0.99995   Jaccard: 0.0066418
NotAllPredictedEffectorTandem.list                 notA: 23025            inA: 6178        notA: 525              inA: 43       odds.ratio: 0.30526     pval:      1    Jaccard: 0.00637
NotBuscoDNA.list                   notA: 27330            inA: 1873        notA: 550              inA: 18       odds.ratio: 0.47755     pval: 0.99979   Jaccard: 0.007374       is.teste
NotBuscoLTR.list                   notA: 27829            inA: 1374        notA: 561              inA: 7        odds.ratio: 0.25273     pval:      1    Jaccard: 0.0036045      is.teste
NotBuscoTandem.list                notA: 22994            inA: 6209        notA: 519              inA: 49       odds.ratio: 0.34964     pval:      1    Jaccard: 0.0072303      is.teste
NotDorsalEffectorLTR.list                  notA: 27815            inA: 1388        notA: 561              inA: 7        odds.ratio: 0.25006     pval:      1    Jaccard: 0.0035787
NotDorsalEffectorTIR.list                  notA: 27319            inA: 1884        notA: 552              inA: 16       odds.ratio: 0.42031     pval: 0.99997   Jaccard: 0.0065253
NotDorsalEffectorTandem.list               notA: 22950            inA: 6253        notA: 522              inA: 46       odds.ratio: 0.32343     pval:      1    Jaccard: 0.0067439
NotKnownEffectorLTR.list                   notA: 27814            inA: 1389        notA: 561              inA: 7        odds.ratio: 0.24987     pval:      1    Jaccard: 0.0035769
NotKnownEffectorTIR.list                   notA: 27322            inA: 1881        notA: 551              inA: 17       odds.ratio: 0.44816     pval: 0.99992   Jaccard: 0.0069416
NotKnownEffectorTandem.list                notA: 22950            inA: 6253        notA: 522              inA: 46       odds.ratio: 0.32343     pval:      1    Jaccard: 0.0067439
NotPredictedEffectorLTR.list               notA: 27828            inA: 1375        notA: 561              inA: 7        odds.ratio: 0.25254     pval:      1    Jaccard: 0.0036027
NotPredictedEffectorTIR.list               notA: 27357            inA: 1846        notA: 551              inA: 17       odds.ratio: 0.45724     pval: 0.99988   Jaccard: 0.0070423
NotPredictedEffectorTandem.list            notA: 23001            inA: 6202        notA: 523              inA: 45       odds.ratio: 0.3191      pval:      1    Jaccard: 0.006647
PredictedEffector.list             notA: 28918            inA: 285         notA: 561              inA: 7        odds.ratio: 1.2661      pval: 0.32415   Jaccard: 0.0082063      is.teste
PredictedEffectorLTR.list                  notA: 29184            inA: 19          notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.teste
PredictedEffectorTIR.list                  notA: 29152            inA: 51          notA: 567              inA: 1        odds.ratio: 1.0081      pval: 0.63306   Jaccard: 0.0016155
PredictedEffectorTandem.list               notA: 29075            inA: 128         notA: 560              inA: 8        odds.ratio: 3.2447      pval: 0.0046783 Jaccard: 0.011494
RepeatAffectedGenesNoSimple.list                   notA: 16989            inA: 12214       notA: 425              inA: 143      odds.ratio: 0.46803     pval:      1    Jaccard: 0.01118
SNPDensityAll0SNPS.list            notA: 24614            inA: 4589        notA: 545              inA: 23       odds.ratio: 0.22636     pval:      1    Jaccard: 0.00446        is.teste
SNPDensityAll10percAfterRemoveZeros.list                   notA: 26732            inA: 2471        notA: 523              inA: 45       odds.ratio: 0.93083     pval: 0.69804   Jaccard:
SNPDensityAll90Percentile.list             notA: 26270            inA: 2933        notA: 532              inA: 36       odds.ratio: 0.6061      pval: 0.99923   Jaccard: 0.010283
Secretome.gene.list                notA: 27081            inA: 2122        notA: 468              inA: 100      odds.ratio: 2.7267      pval: 6.4594e-16        Jaccard: 0.037175
SupportedIRFMergeClassified.list                   notA: 27306            inA: 1897        notA: 550              inA: 18       odds.ratio: 0.4711      pval: 0.99984   Jaccard: 0.00730
TandemDuplicationGene.list                 notA: 22536            inA: 6667        notA: 505              inA: 63       odds.ratio: 0.42168     pval:      1    Jaccard: 0.0087077

```

# get the low snp density gene overlap ready
```
ls -1 *list |awk '{print "go.obj <- newGeneOverlap("$1"$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)"}' |awk '{print $0"\ngo.obj"NR" <- testGeneOverlap(go.obj)"}' |less

go.obj <- newGeneOverlap(x121EffectorOnlyNewGeneNames.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj1 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllGenes.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj2 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllHGTRevised.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj3 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj4 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj5 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj6 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectors.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj7 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoDNA.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj8 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj9 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj10 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj11 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj12 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj13 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsallikeGenesFixed.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj14 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DuplicatedBuscosGene.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj15 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family976.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj16 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family1265.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj17 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTRevised.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj18 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSignificant.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj19 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj20 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj21 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj22 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(LowSignificant.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj23 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NonRedundantLtrRetroelementMerge.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj24 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj25 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj26 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj27 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoDNA.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj28 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj29 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj30 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj31 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj32 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj33 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj34 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj35 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj36 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj37 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj38 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj39 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffector.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj40 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj41 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj42 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj43 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(RepeatAffectedGenesNoSimple.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj44 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll0SNPS.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj45 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll10percAfterRemoveZeros.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj46 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll90Percentile.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj47 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Secretome.gene.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj48 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SupportedIRFMergeClassified.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj49 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(TandemDuplicationGene.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=AllGenes.list)
go.obj50 <- testGeneOverlap(go.obj)

AllCombine <- c(go.obj1, go.obj2, go.obj3, go.obj4, go.obj5, go.obj6, go.obj7, go.obj8, go.obj9, go.obj10, go.obj11, go.obj12, go.obj13, go.obj14, go.obj15, go.obj16, go.obj17, go.obj18, go.obj19, go.obj20, go.obj21, go.obj22, go.obj23, go.obj24, go.obj25, go.obj26, go.obj27, go.obj28, go.obj29, go.obj30, go.obj31, go.obj32, go.obj33, go.obj34, go.obj35, go.obj36, go.obj37, go.obj38, go.obj39, go.obj40, go.obj41, go.obj42, go.obj43, go.obj44, go.obj45, go.obj46, go.obj47, go.obj48, go.obj49, go.obj50)
library(RJSONIO)
exportJSON <- toJSON(AllCombine)
write(exportJSON,"SNPDensityAll10percAfterRemoveZeros.out")
ls -1 *list |paste - <(less SNPDensityAll10percAfterRemoveZeros.out|sed 's/\[/\t/g' |awk 'NF<15' |grep "\"" | sed '0~5 s/$/\nwait/g' |tr "\n" "\t" |sed 's/wait/\n/g' |sed 's/,/\t/g' |sed 's/"//g' |cut -f 1-7,9- |sed 's/}//g') |less -S
```
### Results of overlap of low snp density
```

121EffectorOnlyNewGeneNames.list           notA: 27138            inA: 117         notA: 2512             inA: 4           ]    pval: 0.99319   Jaccard: 0.0015192      is.tested: true
AllGenes.list              notA: 27254            inA: 1           notA: 2516             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
AllHGTRevised.list                 notA: 25866            inA: 1389        notA: 2347             inA: 169      odds.ratio: 1.3409      pval: 0.00042846        Jaccard: 0.043278
AllPredictedEffectorLTR.list               notA: 27238            inA: 17          notA: 2511             inA: 5        odds.ratio: 3.1902      pval: 0.033545  Jaccard: 0.0019739
AllPredictedEffectorTIR.list               notA: 27205            inA: 50          notA: 2508             inA: 8        odds.ratio: 1.7355      pval: 0.11367   Jaccard: 0.0031177
AllPredictedEffectorTandem.list            notA: 27132            inA: 123         notA: 2503             inA: 13       odds.ratio: 1.1457      pval: 0.36351   Jaccard: 0.0049261
AllPredictedEffectors.list                 notA: 26859            inA: 396         notA: 2481             inA: 35       odds.ratio: 0.95683     pval: 0.62347   Jaccard: 0.012019
BuscoDNA.list              notA: 27233            inA: 22          notA: 2514             inA: 2        odds.ratio: 0.98478     pval: 0.61381   Jaccard: 0.00078802     is.tested: true
BuscoLTR.list              notA: 27236            inA: 19          notA: 2515             inA: 1        odds.ratio:   0.57      pval: 0.82908   Jaccard: 0.00039448     is.tested: true
BuscoTandem.list                   notA: 27181            inA: 74          notA: 2511             inA: 5        odds.ratio: 0.73138     pval: 0.80798   Jaccard: 0.0019305      is.teste
DorsalEffectorLTR.list             notA: 27251            inA: 4           notA: 2514             inA: 2        odds.ratio: 5.4188      pval: 0.085166  Jaccard: 0.00079365     is.teste
DorsalEffectorTIR.list             notA: 27240            inA: 15          notA: 2516             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
DorsalEffectorTandem.list                  notA: 27218            inA: 37          notA: 2512             inA: 4        odds.ratio: 1.1714      pval: 0.46027   Jaccard: 0.0015668
Dorsal-likeGenesFixed.list                 notA: 27107            inA: 148         notA: 2504             inA: 12       odds.ratio: 0.87774     pval: 0.70861   Jaccard: 0.0045045
DuplicatedBuscosGene.list                  notA: 26943            inA: 312         notA: 2485             inA: 31       odds.ratio: 1.0773      pval: 0.37475   Jaccard: 0.010962
Family976.list             notA: 27160            inA: 95          notA: 2510             inA: 6        odds.ratio: 0.68345     pval: 0.86488   Jaccard: 0.002298       is.tested: true
Family1265.list            notA: 27152            inA: 103         notA: 2507             inA: 9        odds.ratio: 0.94635     pval: 0.6123    Jaccard: 0.0034364      is.tested: true
HGTRevised.list            notA: 27120            inA: 135         notA: 2500             inA: 16       odds.ratio: 1.2857      pval: 0.20657   Jaccard: 0.0060355      is.tested: true
HighSignificant.list               notA: 26133            inA: 1122        notA: 2427             inA: 89       odds.ratio: 0.8541      pval: 0.93013   Jaccard: 0.024464       is.teste
KnownEffectorLTR.list              notA: 27250            inA: 5           notA: 2516             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
KnownEffectorTIR.list              notA: 27240            inA: 15          notA: 2514             inA: 2        odds.ratio: 1.4447      pval: 0.42735   Jaccard: 0.0007902      is.teste
KnownEffectorTandem.list                   notA: 27217            inA: 38          notA: 2513             inA: 3        odds.ratio: 0.85503     pval: 0.68495   Jaccard: 0.0011746
LowSignificant.list                notA: 26732            inA: 523         notA: 2471             inA: 45       odds.ratio: 0.93083     pval: 0.69804   Jaccard: 0.014808       is.teste
NonRedundantLtrRetroelementMerge.list              notA: 25977            inA: 1278        notA: 2393             inA: 123      odds.ratio: 1.0448      pval: 0.33958   Jaccard: 0.03242
NotAllPredictedEffectorLTR.list            notA: 25994            inA: 1261        notA: 2398             inA: 118      odds.ratio: 1.0144      pval: 0.45726   Jaccard: 0.031242
NotAllPredictedEffectorTIR.list            notA: 25568            inA: 1687        notA: 2346             inA: 170      odds.ratio: 1.0982      pval: 0.13991   Jaccard: 0.040447
NotAllPredictedEffectorTandem.list                 notA: 21618            inA: 5637        notA: 1932             inA: 584      odds.ratio: 1.1592      pval: 0.001698  Jaccard: 0.07163
NotBuscoDNA.list                   notA: 25540            inA: 1715        notA: 2340             inA: 176      odds.ratio: 1.1201      pval: 0.091333  Jaccard: 0.041598       is.teste
NotBuscoLTR.list                   notA: 25996            inA: 1259        notA: 2394             inA: 122      odds.ratio: 1.0523      pval: 0.31417   Jaccard: 0.032318       is.teste
NotBuscoTandem.list                notA: 21588            inA: 5667        notA: 1925             inA: 591      odds.ratio: 1.1696      pval: 0.00091418        Jaccard: 0.072223
NotDorsalEffectorLTR.list                  notA: 25981            inA: 1274        notA: 2395             inA: 121      odds.ratio: 1.0303      pval: 0.39417   Jaccard: 0.031926
NotDorsalEffectorTIR.list                  notA: 25533            inA: 1722        notA: 2338             inA: 178      odds.ratio: 1.1288      pval: 0.075997  Jaccard: 0.042001
NotDorsalEffectorTandem.list               notA: 21548            inA: 5707        notA: 1924             inA: 592      odds.ratio: 1.1617      pval: 0.0014094 Jaccard: 0.071993
NotKnownEffectorLTR.list                   notA: 25982            inA: 1273        notA: 2393             inA: 123      odds.ratio: 1.0491      pval: 0.32437   Jaccard: 0.032462
NotKnownEffectorTIR.list                   notA: 25533            inA: 1722        notA: 2340             inA: 176      odds.ratio: 1.1152      pval: 0.10006   Jaccard: 0.041529
NotKnownEffectorTandem.list                notA: 21549            inA: 5706        notA: 1923             inA: 593      odds.ratio: 1.1645      pval: 0.001196  Jaccard: 0.072124
NotPredictedEffectorLTR.list               notA: 25993            inA: 1262        notA: 2396             inA: 120      odds.ratio: 1.0316      pval: 0.38998   Jaccard: 0.031763
NotPredictedEffectorTIR.list               notA: 25562            inA: 1693        notA: 2346             inA: 170      odds.ratio: 1.0941      pval:   0.15    Jaccard: 0.04039
NotPredictedEffectorTandem.list            notA: 21597            inA: 5658        notA: 1927             inA: 589      odds.ratio: 1.1667      pval: 0.0010875 Jaccard: 0.072058
PredictedEffector.list             notA: 26985            inA: 270         notA: 2494             inA: 22       odds.ratio: 0.88162     pval: 0.74404   Jaccard: 0.0078966      is.teste
PredictedEffectorLTR.list                  notA: 27239            inA: 16          notA: 2513             inA: 3        odds.ratio: 2.0323      pval: 0.21327   Jaccard: 0.0011848
PredictedEffectorTIR.list                  notA: 27211            inA: 44          notA: 2508             inA: 8        odds.ratio: 1.9726      pval: 0.069072  Jaccard: 0.003125
PredictedEffectorTandem.list               notA: 27132            inA: 123         notA: 2503             inA: 13       odds.ratio: 1.1457      pval: 0.36351   Jaccard: 0.0049261
RepeatAffectedGenesNoSimple.list                   notA: 15733            inA: 11522       notA: 1681             inA: 835      odds.ratio: 0.67828     pval:      1    Jaccard: 0.05948
SNPDensityAll0SNPS.list            notA: 22643            inA: 4612        notA: 2516             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
SNPDensityAll10percAfterRemoveZeros.list                   notA: 27255            inA: 0           notA: 0                inA: 2516     odds.ratio:   Infinity  pval:      0    Jaccard:
SNPDensityAll90Percentile.list             notA: 24286            inA: 2969        notA: 2516             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.teste
Secretome.gene.list                notA: 25189            inA: 2066        notA: 2360             inA: 156      odds.ratio: 0.80593     pval: 0.99556   Jaccard: 0.034046       is.teste
SupportedIRFMergeClassified.list                   notA: 25518            inA: 1737        notA: 2338             inA: 178      odds.ratio: 1.1184      pval: 0.092974  Jaccard: 0.04185
TandemDuplicationGene.list                 notA: 21165            inA: 6090        notA: 1876             inA: 640      odds.ratio: 1.1857      pval: 0.00024778        Jaccard: 0.07436

```

# Start gene overlap for genes in the 90 percentile of snps
```
ls -1 *list |awk '{print "go.obj <- newGeneOverlap("$1"$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)"}' |awk '{print $0"\ngo.obj"NR" <- testGeneOverlap(go.obj)"}' |less

go.obj <- newGeneOverlap(x121EffectorOnlyNewGeneNames.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj1 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllGenes.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj2 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllHGTRevised.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj3 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj4 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj5 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj6 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectors.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj7 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoDNA.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj8 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj9 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj10 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj11 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj12 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj13 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsallikeGenesFixed.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj14 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DuplicatedBuscosGene.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj15 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family976.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj16 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family1265.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj17 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTRevised.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj18 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSignificant.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj19 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj20 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj21 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj22 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(LowSignificant.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj23 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NonRedundantLtrRetroelementMerge.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj24 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj25 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj26 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj27 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoDNA.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj28 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj29 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj30 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj31 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj32 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj33 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj34 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj35 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj36 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj37 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj38 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj39 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffector.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj40 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj41 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj42 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj43 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(RepeatAffectedGenesNoSimple.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj44 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll0SNPS.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj45 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll10percAfterRemoveZeros.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj46 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll90Percentile.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj47 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Secretome.gene.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj48 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SupportedIRFMergeClassified.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj49 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(TandemDuplicationGene.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=AllGenes.list)
go.obj50 <- testGeneOverlap(go.obj)


AllCombine <- c(go.obj1, go.obj2, go.obj3, go.obj4, go.obj5, go.obj6, go.obj7, go.obj8, go.obj9, go.obj10, go.obj11, go.obj12, go.obj13, go.obj14, go.obj15, go.obj16, go.obj17, go.obj18, go.obj19, go.obj20, go.obj21, go.obj22, go.obj23, go.obj24, go.obj25, go.obj26, go.obj27, go.obj28, go.obj29, go.obj30, go.obj31, go.obj32, go.obj33, go.obj34, go.obj35, go.obj36, go.obj37, go.obj38, go.obj39, go.obj40, go.obj41, go.obj42, go.obj43, go.obj44, go.obj45, go.obj46, go.obj47, go.obj48, go.obj49, go.obj50)
library(RJSONIO)
exportJSON <- toJSON(AllCombine)
write(exportJSON,"SNPDensityAll90Percentile.out")
```

```
ls -1 *list |paste - <(less SNPDensityAll90Percentile.out|sed 's/\[/\t/g' |awk 'NF<15' |grep "\"" | sed '0~5 s/$/\nwait/g' |tr "\n" "\t" |sed 's/wait/\n/g' |sed 's/,/\t/g' |sed 's/"//g' |cut -f 1-7,9- |sed 's/}//g') |less -S
```
### Results from High snp density overlap analysis
```
121EffectorOnlyNewGeneNames.list           notA: 26698            inA: 104         notA: 2952             inA: 17          ]    pval: 0.093069  Jaccard: 0.0055321      is.tested: true
AllGenes.list              notA: 26801            inA: 1           notA: 2969             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
AllHGTRevised.list                 notA: 25389            inA: 1413        notA: 2824             inA: 145      odds.ratio: 0.92258     pval: 0.82732   Jaccard: 0.03309        is.tested: true
AllPredictedEffectorLTR.list               notA: 26781            inA: 21          notA: 2968             inA: 1        odds.ratio: 0.42969     pval: 0.90095   Jaccard: 0.00033445     is.tested: true
AllPredictedEffectorTIR.list               notA: 26751            inA: 51          notA: 2962             inA: 7        odds.ratio: 1.2396      pval: 0.35702   Jaccard: 0.0023179      is.tested: true
AllPredictedEffectorTandem.list            notA: 26680            inA: 122         notA: 2955             inA: 14       odds.ratio: 1.0361      pval: 0.49193   Jaccard: 0.0045293      is.tested: true
AllPredictedEffectors.list                 notA: 26411            inA: 391         notA: 2929             inA: 40       odds.ratio: 0.92246     pval: 0.70865   Jaccard: 0.011905       is.tested: true
BuscoDNA.list              notA: 26781            inA: 21          notA: 2966             inA: 3        odds.ratio: 1.2899      pval: 0.43395   Jaccard: 0.0010033      is.tested: true
BuscoLTR.list              notA: 26786            inA: 16          notA: 2965             inA: 4        odds.ratio: 2.2584      pval: 0.13191   Jaccard: 0.00134        is.tested: true
BuscoTandem.list                   notA: 26738            inA: 64          notA: 2954             inA: 15       odds.ratio: 2.1214      pval: 0.010671  Jaccard: 0.0049456      is.tested: true
DorsalEffectorLTR.list             notA: 26796            inA: 6           notA: 2969             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
DorsalEffectorTIR.list             notA: 26790            inA: 12          notA: 2966             inA: 3        odds.ratio:  2.258      pval: 0.18297   Jaccard: 0.0010064      is.tested: true
DorsalEffectorTandem.list                  notA: 26768            inA: 34          notA: 2962             inA: 7        odds.ratio: 1.8605      pval: 0.10886   Jaccard: 0.002331       is.tested: true
Dorsal-likeGenesFixed.list                 notA: 26660            inA: 142         notA: 2951             inA: 18       odds.ratio: 1.1452      pval: 0.33097   Jaccard: 0.0057859      is.tested: true
DuplicatedBuscosGene.list                  notA: 26495            inA: 307         notA: 2933             inA: 36       odds.ratio: 1.0593      pval: 0.39872   Jaccard: 0.010989       is.tested: true
Family976.list             notA: 26711            inA: 91          notA: 2959             inA: 10       odds.ratio: 0.99198     pval: 0.55838   Jaccard: 0.003268       is.tested: true
Family1265.list            notA: 26698            inA: 104         notA: 2961             inA: 8        odds.ratio: 0.6936      pval: 0.88089   Jaccard: 0.0026033      is.tested: true
HGTRevised.list            notA: 26655            inA: 147         notA: 2965             inA: 4        odds.ratio: 0.24463     pval: 0.99988   Jaccard: 0.0012837      is.tested: true
HighSignificant.list               notA: 25680            inA: 1122        notA: 2880             inA: 89       odds.ratio: 0.7073      pval: 0.99949   Jaccard: 0.021755       is.tested: true
KnownEffectorLTR.list              notA: 26797            inA: 5           notA: 2969             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
KnownEffectorTIR.list              notA: 26788            inA: 14          notA: 2966             inA: 3        odds.ratio: 1.9354      pval: 0.2369    Jaccard: 0.0010057      is.tested: true
KnownEffectorTandem.list                   notA: 26767            inA: 35          notA: 2963             inA: 6        odds.ratio: 1.5486      pval: 0.22079   Jaccard: 0.0019973      is.tested: true
LowSignificant.list                notA: 26270            inA: 532         notA: 2933             inA: 36       odds.ratio: 0.6061      pval: 0.99923   Jaccard: 0.010283       is.tested: true
NonRedundantLtrRetroelementMerge.list              notA: 25661            inA: 1141        notA: 2709             inA: 260      odds.ratio: 2.1584      pval: 1.0265e-23        Jaccard: 0.06326        is.tested:
NotAllPredictedEffectorLTR.list            notA: 25682            inA: 1120        notA: 2710             inA: 259      odds.ratio: 2.1916      pval: 2.1326e-24        Jaccard: 0.063341       is.tested: true
NotAllPredictedEffectorTIR.list            notA: 25209            inA: 1593        notA: 2705             inA: 264      odds.ratio: 1.5444      pval: 1.2563e-09        Jaccard: 0.057869       is.tested: true
NotAllPredictedEffectorTandem.list                 notA: 21396            inA: 5406        notA: 2154             inA: 815      odds.ratio: 1.4975      pval: 1.5573e-19        Jaccard: 0.097313       is.tested:
NotBuscoDNA.list                   notA: 25179            inA: 1623        notA: 2701             inA: 268      odds.ratio: 1.5393      pval: 1.265e-09 Jaccard: 0.058362       is.tested: true
NotBuscoLTR.list                   notA: 25677            inA: 1125        notA: 2713             inA: 256      odds.ratio: 2.1536      pval: 2.7539e-23        Jaccard: 0.062531       is.tested: true
NotBuscoTandem.list                notA: 21358            inA: 5444        notA: 2155             inA: 814      odds.ratio: 1.4819      pval: 1.2002e-18        Jaccard: 0.096755       is.tested: true
NotDorsalEffectorLTR.list                  notA: 25667            inA: 1135        notA: 2709             inA: 260      odds.ratio: 2.1703      pval: 5.4307e-24        Jaccard: 0.063353       is.tested: true
NotDorsalEffectorTIR.list                  notA: 25170            inA: 1632        notA: 2701             inA: 268      odds.ratio: 1.5303      pval: 2.0133e-09        Jaccard: 0.058248       is.tested: true
NotDorsalEffectorTandem.list               notA: 21325            inA: 5477        notA: 2147             inA: 822      odds.ratio: 1.4907      pval: 2.963e-19 Jaccard: 0.097324       is.tested: true
NotKnownEffectorLTR.list                   notA: 25666            inA: 1136        notA: 2709             inA: 260      odds.ratio: 2.1684      pval: 6.0412e-24        Jaccard: 0.063337       is.tested: true
NotKnownEffectorTIR.list                   notA: 25172            inA: 1630        notA: 2701             inA: 268      odds.ratio: 1.5323      pval: 1.8169e-09        Jaccard: 0.058274       is.tested: true
NotKnownEffectorTandem.list                notA: 21326            inA: 5476        notA: 2146             inA: 823      odds.ratio: 1.4935      pval: 1.9777e-19        Jaccard: 0.097454       is.tested: true
NotPredictedEffectorLTR.list               notA: 25679            inA: 1123        notA: 2710             inA: 259      odds.ratio: 2.1854      pval: 2.9499e-24        Jaccard: 0.063294       is.tested: true
NotPredictedEffectorTIR.list               notA: 25203            inA: 1599        notA: 2705             inA: 264      odds.ratio: 1.5383      pval: 1.7189e-09        Jaccard: 0.057793       is.tested: true
NotPredictedEffectorTandem.list            notA: 21371            inA: 5431        notA: 2153             inA: 816      odds.ratio: 1.4914      pval: 3.3172e-19        Jaccard: 0.097143       is.tested: true
PredictedEffector.list             notA: 26544            inA: 258         notA: 2935             inA: 34       odds.ratio: 1.1918      pval: 0.19305   Jaccard: 0.010536       is.tested: true
PredictedEffectorLTR.list                  notA: 26784            inA: 18          notA: 2968             inA: 1        odds.ratio: 0.50136     pval: 0.86422   Jaccard: 0.00033478     is.tested: true
PredictedEffectorTIR.list                  notA: 26757            inA: 45          notA: 2962             inA: 7        odds.ratio: 1.4052      pval: 0.25856   Jaccard: 0.0023225      is.tested: true
PredictedEffectorTandem.list               notA: 26680            inA: 122         notA: 2955             inA: 14       odds.ratio: 1.0361      pval: 0.49193   Jaccard: 0.0045293      is.tested: true
RepeatAffectedGenesNoSimple.list                   notA: 16334            inA: 10468       notA: 1080             inA: 1889     odds.ratio: 2.7291      pval: 2.0423e-144       Jaccard: 0.14058        is.tested:
SNPDensityAll0SNPS.list            notA: 22190            inA: 4612        notA: 2969             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
SNPDensityAll10percAfterRemoveZeros.list                   notA: 24286            inA: 2516        notA: 2969             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
SNPDensityAll90Percentile.list             notA: 26802            inA: 0           notA: 0                inA: 2969     odds.ratio:   Infinity  pval:      0    Jaccard:      1 is.tested: true
Secretome.gene.list                notA: 24756            inA: 2046        notA: 2793             inA: 176      odds.ratio: 0.76247     pval: 0.99976   Jaccard: 0.035095       is.tested: true
SupportedIRFMergeClassified.list                   notA: 25158            inA: 1644        notA: 2698             inA: 271      odds.ratio: 1.5371      pval: 1.1741e-09        Jaccard: 0.058747       is.tested:
TandemDuplicationGene.list                 notA: 20920            inA: 5882        notA: 2121             inA: 848      odds.ratio:  1.422      pval: 8.6615e-16        Jaccard: 0.095808       is.tested: true
```


### Running gene overlap for gene categories with genes with 0 snps
```
ls -1 *list |awk '{print "go.obj <- newGeneOverlap("$1"$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)"}' |awk '{print $0"\ngo.obj"NR" <- testGeneOverlap(go.obj)"}' |less

go.obj <- newGeneOverlap(121EffectorOnlyNewGeneNames.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj1 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllGenes.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj2 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllHGTRevised.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj3 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj4 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj5 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj6 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectors.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj7 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoDNA.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj8 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj9 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj10 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj11 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj12 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj13 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Dorsal-likeGenesFixed.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj14 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DuplicatedBuscosGene.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj15 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family976.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj16 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family1265.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj17 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTRevised.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj18 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSignificant.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj19 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj20 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj21 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj22 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(LowSignificant.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj23 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NonRedundantLtrRetroelementMerge.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj24 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj25 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj26 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj27 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoDNA.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj28 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj29 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj30 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj31 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj32 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj33 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj34 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj35 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj36 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj37 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj38 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj39 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffector.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj40 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj41 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj42 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj43 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(RepeatAffectedGenesNoSimple.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj44 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll0SNPS.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj45 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll10percAfterRemoveZeros.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj46 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll90Percentile.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj47 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Secretome.gene.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj48 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SupportedIRFMergeClassified.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj49 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(TandemDuplicationGene.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=AllGenes.list)
go.obj50 <- testGeneOverlap(go.obj)



AllCombine <- c(go.obj1, go.obj2, go.obj3, go.obj4, go.obj5, go.obj6, go.obj7, go.obj8, go.obj9, go.obj10, go.obj11, go.obj12, go.obj13, go.obj14, go.obj15, go.obj16, go.obj17, go.obj18, go.obj19, go.obj20, go.obj21, go.obj22, go.obj23, go.obj24, go.obj25, go.obj26, go.obj27, go.obj28, go.obj29, go.obj30, go.obj31, go.obj32, go.obj33, go.obj34, go.obj35, go.obj36, go.obj37, go.obj38, go.obj39, go.obj40, go.obj41, go.obj42, go.obj43, go.obj44, go.obj45, go.obj46, go.obj47, go.obj48, go.obj49, go.obj50)
library(RJSONIO)
exportJSON <- toJSON(AllCombine)
write(exportJSON,"SNPDensityAll0SNPs.out")

ls -1 *list |paste - <(less SNPDensityAll0SNPs.out|sed 's/\[/\t/g' |awk 'NF<15' |grep "\"" | sed '0~5 s/$/\nwait/g' |tr "\n" "\t" |sed 's/wait/\n/g' |sed 's/,/\t/g' |sed 's/"//g' |cut -f 1-7,9- |sed 's/}//g') |less -S
```

### Results of zero snp gene overlapping
```
121EffectorOnlyNewGeneNames.list           notA: 20920            inA: 5882        notA: 2121             inA: 848         ]    pval: 8.6615e-16        Jaccard: 0.095808       is.tested: true
AllGenes.list              notA: 25158            inA: 1           notA: 4612             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
AllHGTRevised.list                 notA: 23861            inA: 1298        notA: 4352             inA: 260      odds.ratio: 1.0983      pval: 0.096817  Jaccard: 0.043993       is.tested: true
AllPredictedEffectorLTR.list               notA: 25145            inA: 14          notA: 4604             inA: 8        odds.ratio: 3.1207      pval: 0.01373   Jaccard: 0.0017294      is.tested: true
AllPredictedEffectorTIR.list               notA: 25112            inA: 47          notA: 4601             inA: 11       odds.ratio: 1.2774      pval: 0.28128   Jaccard: 0.002361       is.tested: true
AllPredictedEffectorTandem.list            notA: 25060            inA: 99          notA: 4575             inA: 37       odds.ratio: 2.0472      pval: 0.00032212        Jaccard: 0.007854       is.tested: true
AllPredictedEffectors.list                 notA: 24806            inA: 353         notA: 4534             inA: 78       odds.ratio: 1.2089      pval: 0.077142  Jaccard: 0.01571        is.tested: true
BuscoDNA.list              notA: 25142            inA: 17          notA: 4605             inA: 7        odds.ratio:  2.248      pval: 0.066099  Jaccard: 0.0015122      is.tested: true
BuscoLTR.list              notA: 25148            inA: 11          notA: 4603             inA: 9        odds.ratio: 4.4699      pval: 0.0016742 Jaccard: 0.0019468      is.tested: true
BuscoTandem.list                   notA: 25094            inA: 65          notA: 4598             inA: 14       odds.ratio: 1.1754      pval: 0.33627   Jaccard: 0.0029934      is.tested: true
DorsalEffectorLTR.list             notA: 25153            inA: 6           notA: 4612             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
DorsalEffectorTIR.list             notA: 25144            inA: 15          notA: 4612             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
DorsalEffectorTandem.list                  notA: 25120            inA: 39          notA: 4610             inA: 2        odds.ratio: 0.27944     pval: 0.99146   Jaccard: 0.00043002     is.tested: true
Dorsal-likeGenesFixed.list                 notA: 25120            inA: 39          notA: 4610             inA: 2        odds.ratio: 0.27944     pval: 0.99146   Jaccard: 0.00043002     is.tested: true
DuplicatedBuscosGene.list                  notA: 24878            inA: 281         notA: 4550             inA: 62       odds.ratio: 1.2064      pval: 0.10623   Jaccard: 0.012671       is.tested: true
Family976.list             notA: 25082            inA: 77          notA: 4588             inA: 24       odds.ratio: 1.7039      pval: 0.019185  Jaccard: 0.0051184      is.tested: true
Family1265.list            notA: 25081            inA: 78          notA: 4578             inA: 34       odds.ratio: 2.3881      pval: 5.4716e-05        Jaccard: 0.0072495      is.tested: true
HGTRevised.list            notA: 25034            inA: 125         notA: 4586             inA: 26       odds.ratio: 1.1354      pval: 0.31032   Jaccard: 0.0054887      is.tested: true
HighSignificant.list               notA: 23974            inA: 1185        notA: 4586             inA: 26       odds.ratio: 0.1147      pval:      1    Jaccard: 0.0044851      is.tested: true
KnownEffectorLTR.list              notA: 25156            inA: 3           notA: 4610             inA: 2        odds.ratio: 3.6377      pval: 0.1739    Jaccard: 0.00043337     is.tested: true
KnownEffectorTIR.list              notA: 25144            inA: 15          notA: 4610             inA: 2        odds.ratio: 0.72724     pval: 0.76468   Jaccard: 0.00043225     is.tested: true
KnownEffectorTandem.list                   notA: 25126            inA: 33          notA: 4604             inA: 8        odds.ratio:  1.323      pval:  0.297    Jaccard: 0.0017223      is.tested: true
LowSignificant.list                notA: 24614            inA: 545         notA: 4589             inA: 23       odds.ratio: 0.22636     pval:      1    Jaccard: 0.00446        is.tested: true
NonRedundantLtrRetroelementMerge.list              notA: 24055            inA: 1104        notA: 4315             inA: 297      odds.ratio: 1.4997      pval: 3.8962e-09        Jaccard: 0.051959       is.tested:
NotAllPredictedEffectorLTR.list            notA: 24069            inA: 1090        notA: 4323             inA: 289      odds.ratio: 1.4762      pval: 2.0179e-08        Jaccard: 0.050684       is.tested: true
NotAllPredictedEffectorTIR.list            notA: 23752            inA: 1407        notA: 4162             inA: 450      odds.ratio: 1.8251      pval: 2.4807e-24        Jaccard: 0.074763       is.tested: true
NotAllPredictedEffectorTandem.list                 notA: 20280            inA: 4879        notA: 3270             inA: 1342     odds.ratio: 1.7058      pval: 2.5775e-47        Jaccard: 0.1414 is.tested: true
NotBuscoDNA.list                   notA: 23722            inA: 1437        notA: 4158             inA: 454      odds.ratio: 1.8024      pval: 1.1368e-23        Jaccard: 0.075054       is.tested: true
NotBuscoLTR.list                   notA: 24066            inA: 1093        notA: 4324             inA: 288      odds.ratio: 1.4665      pval: 3.4587e-08        Jaccard: 0.050482       is.tested: true
NotBuscoTandem.list                notA: 20256            inA: 4903        notA: 3257             inA: 1355     odds.ratio: 1.7187      pval: 7.6354e-49        Jaccard: 0.14241        is.tested: true
NotDorsalEffectorLTR.list                  notA: 24061            inA: 1098        notA: 4315             inA: 297      odds.ratio: 1.5083      pval: 2.4763e-09        Jaccard: 0.052014       is.tested: true
NotDorsalEffectorTIR.list                  notA: 23720            inA: 1439        notA: 4151             inA: 461      odds.ratio: 1.8305      pval: 4.7813e-25        Jaccard: 0.076186       is.tested: true
NotDorsalEffectorTandem.list               notA: 20227            inA: 4932        notA: 3245             inA: 1367     odds.ratio: 1.7276      pval: 5.509e-50 Jaccard: 0.14323        is.tested: true
NotKnownEffectorLTR.list                   notA: 24058            inA: 1101        notA: 4317             inA: 295      odds.ratio: 1.4932      pval: 6.1294e-09        Jaccard: 0.051637       is.tested: true
NotKnownEffectorTIR.list                   notA: 23720            inA: 1439        notA: 4153             inA: 459      odds.ratio: 1.8217      pval: 1.2733e-24        Jaccard: 0.075855       is.tested: true
NotKnownEffectorTandem.list                notA: 20221            inA: 4938        notA: 3251             inA: 1361     odds.ratio: 1.7142      pval: 1.4478e-48        Jaccard: 0.14251        is.tested: true
NotPredictedEffectorLTR.list               notA: 24066            inA: 1093        notA: 4323             inA: 289      odds.ratio: 1.4719      pval: 2.5008e-08        Jaccard: 0.050657       is.tested: true
NotPredictedEffectorTIR.list               notA: 23746            inA: 1413        notA: 4162             inA: 450      odds.ratio: 1.8169      pval: 4.9852e-24        Jaccard: 0.074689       is.tested: true
NotPredictedEffectorTandem.list            notA: 20255            inA: 4904        notA: 3269             inA: 1343     odds.ratio: 1.6968      pval: 1.7104e-46        Jaccard: 0.14113        is.tested: true
PredictedEffector.list             notA: 24936            inA: 223         notA: 4543             inA: 69       odds.ratio: 1.6983      pval: 0.00016941        Jaccard: 0.014271       is.tested: true
PredictedEffectorLTR.list                  notA: 25148            inA: 11          notA: 4604             inA: 8        odds.ratio: 3.9723      pval: 0.004993  Jaccard: 0.0017305      is.tested: true
PredictedEffectorTIR.list                  notA: 25118            inA: 41          notA: 4601             inA: 11       odds.ratio: 1.4647      pval: 0.17231   Jaccard: 0.0023641      is.tested: true
PredictedEffectorTandem.list               notA: 25060            inA: 99          notA: 4575             inA: 37       odds.ratio: 2.0472      pval: 0.00032212        Jaccard: 0.007854       is.tested: true
RepeatAffectedGenesNoSimple.list                   notA: 14920            inA: 10239       notA: 2494             inA: 2118     odds.ratio: 1.2375      pval: 2.3104e-11        Jaccard: 0.14262        is.tested:
SNPDensityAll0SNPS.list            notA: 25159            inA: 0           notA: 0                inA: 4612     odds.ratio:   Infinity  pval:      0    Jaccard:      1 is.tested: true
SNPDensityAll10percAfterRemoveZeros.list                   notA: 22643            inA: 2516        notA: 4612             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
SNPDensityAll90Percentile.list             notA: 22190            inA: 2969        notA: 4612             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
Secretome.gene.list                notA: 23166            inA: 1993        notA: 4383             inA: 229      odds.ratio: 0.60732     pval:      1    Jaccard: 0.034671       is.tested: true
SupportedIRFMergeClassified.list                   notA: 23705            inA: 1454        notA: 4151             inA: 461      odds.ratio: 1.8105      pval: 2.7293e-24        Jaccard: 0.075997       is.tested:
TandemDuplicationGene.list                 notA: 19888            inA: 5271        notA: 3153             inA: 1459     odds.ratio: 1.7458      pval: 5.4809e-54        Jaccard: 0.14763        is.tested: true

```
