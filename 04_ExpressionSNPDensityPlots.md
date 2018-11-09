# Expression and SNP density gene overlap analyses
These are primarily sets of genes that I developed from multiple analyses on the genome.

This is a modification of this pipeline to include the updated HGT candidates.
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/26_ExpressionSets
library(GeneOverlap)
sed 's/\./\t/2' ../30_BenNewHGT/AllHGT.list |awk '{print $1}' >AllHGTRevised.list
sed 's/\./\t/2' ../30_BenNewHGT/HighConfHGT.list |awk '{print $1}' >HGTRevised.list

ls -1 *list|awk '{print $1" <- read.table(\""$1"\")"}' |less
```
### High Expression 90 percentile gene list
```
x121EffectorOnlyNewGeneNames.list <- read.table("121EffectorOnlyNewGeneNames.list")
AllGenes.list <- read.table("AllGenes.list")
AllHGTRevised.list <- read.table("AllHGTRevised.list")
AllPredictedEffectorLTR.list <- read.table("AllPredictedEffectorLTR.list")
AllPredictedEffectorTIR.list <- read.table("AllPredictedEffectorTIR.list")
AllPredictedEffectorTandem.list <- read.table("AllPredictedEffectorTandem.list")
AllPredictedEffectors.list <- read.table("AllPredictedEffectors.list")
BuscoDNA.list <- read.table("BuscoDNA.list")
BuscoLTR.list <- read.table("BuscoLTR.list")
BuscoTandem.list <- read.table("BuscoTandem.list")
DorsalEffectorLTR.list <- read.table("DorsalEffectorLTR.list")
DorsalEffectorTIR.list <- read.table("DorsalEffectorTIR.list")
DorsalEffectorTandem.list <- read.table("DorsalEffectorTandem.list")
DorsallikeGenesFixed.list <- read.table("Dorsal-likeGenesFixed.list")
DuplicatedBuscosGene.list <- read.table("DuplicatedBuscosGene.list")
Family976.list <- read.table("Family976.list")
Family1265.list <- read.table("Family1265.list")
HGT.list <- read.table("HGT.list")
HGTNewNames.list <- read.table("HGTNewNames.list")
HGTRevised.list <- read.table("HGTRevised.list")
HighSNPDensity.list <- read.table("HighSNPDensity.list")
HighSNPDensity1Dev.list <- read.table("HighSNPDensity1Dev.list")
HighSNPDensity2Dev.list <- read.table("HighSNPDensity2Dev.list")
HighSNPDensity10Perc.list <- read.table("HighSNPDensity10Perc.list")
HighSignificant.list <- read.table("HighSignificant.list")
KnownEffectorLTR.list <- read.table("KnownEffectorLTR.list")
KnownEffectorTIR.list <- read.table("KnownEffectorTIR.list")
KnownEffectorTandem.list <- read.table("KnownEffectorTandem.list")
LowSignificant.list <- read.table("LowSignificant.list")
NonRedundantLtrRetroelementMerge.list <- read.table("NonRedundantLtrRetroelementMerge.list")
NotAllPredictedEffectorLTR.list <- read.table("NotAllPredictedEffectorLTR.list")
NotAllPredictedEffectorTIR.list <- read.table("NotAllPredictedEffectorTIR.list")
NotAllPredictedEffectorTandem.list <- read.table("NotAllPredictedEffectorTandem.list")
NotBuscoDNA.list <- read.table("NotBuscoDNA.list")
NotBuscoLTR.list <- read.table("NotBuscoLTR.list")
NotBuscoTandem.list <- read.table("NotBuscoTandem.list")
NotDorsalEffectorLTR.list <- read.table("NotDorsalEffectorLTR.list")
NotDorsalEffectorTIR.list <- read.table("NotDorsalEffectorTIR.list")
NotDorsalEffectorTandem.list <- read.table("NotDorsalEffectorTandem.list")
NotKnownEffectorLTR.list <- read.table("NotKnownEffectorLTR.list")
NotKnownEffectorTIR.list <- read.table("NotKnownEffectorTIR.list")
NotKnownEffectorTandem.list <- read.table("NotKnownEffectorTandem.list")
NotPredictedEffectorLTR.list <- read.table("NotPredictedEffectorLTR.list")
NotPredictedEffectorTIR.list <- read.table("NotPredictedEffectorTIR.list")
NotPredictedEffectorTandem.list <- read.table("NotPredictedEffectorTandem.list")
PredictedEffector.list <- read.table("PredictedEffector.list")
PredictedEffectorLTR.list <- read.table("PredictedEffectorLTR.list")
PredictedEffectorTIR.list <- read.table("PredictedEffectorTIR.list")
PredictedEffectorTandem.list <- read.table("PredictedEffectorTandem.list")
RepeatAffectedGenesNoSimple.list <- read.table("RepeatAffectedGenesNoSimple.list")
Secretome.gene.list <- read.table("Secretome.gene.list")
SupportedIRFMergeClassified.list <- read.table("SupportedIRFMergeClassified.list")
tandem.gene.list <- read.table("tandem.gene.list")
```

```
ls -1 *list |awk '{print "go.obj <- newGeneOverlap("$1"$V1, HighSignificant.list$V1, genome.size=All.genes)"}' |awk '{print $0"\ngo.obj"NR" <- testGeneOverlap(go.obj)"}' |less
Output of the above
##################################
go.obj <- newGeneOverlap(x121EffectorOnlyNewGeneNames.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj1 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllGenes.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj2 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllHGTRevised.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj3 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj4 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj5 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj6 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectors.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj7 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoDNA.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj8 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoLTR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj9 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoTandem.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj10 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj11 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj12 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj13 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsallikeGenesFixed.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj14 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DuplicatedBuscosGene.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj15 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family976.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj16 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family1265.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj17 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGT.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj18 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTNewNames.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj19 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTRevised.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj20 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSNPDensity.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj21 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSNPDensity1Dev.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj22 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSNPDensity2Dev.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj23 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSNPDensity10Perc.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj24 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSignificant.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj25 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj26 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj27 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj28 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(LowSignificant.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj29 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NonRedundantLtrRetroelementMerge.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj30 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj31 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj32 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj33 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoDNA.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj34 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoLTR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj35 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoTandem.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj36 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj37 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj38 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj39 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj40 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj41 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj42 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj43 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj44 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj45 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffector.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj46 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorLTR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj47 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTIR.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj48 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTandem.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj49 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(RepeatAffectedGenesNoSimple.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj50 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Secretome.gene.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj51 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SupportedIRFMergeClassified.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj52 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(tandem.gene.list$V1, HighSignificant.list$V1, genome.size=All.genes)
go.obj53 <- testGeneOverlap(go.obj)
##################################
```
Last few commands to get the test Blast2go_finished
```
AllCombine <- c(go.obj1, go.obj2, go.obj3, go.obj4, go.obj5, go.obj6, go.obj7, go.obj8, go.obj9, go.obj10, go.obj11, go.obj12, go.obj13, go.obj14, go.obj15, go.obj16, go.obj17, go.obj18, go.obj19, go.obj20, go.obj21, go.obj22, go.obj23, go.obj24, go.obj25, go.obj26, go.obj27, go.obj28, go.obj29, go.obj30, go.obj31, go.obj32, go.obj33, go.obj34, go.obj35, go.obj36, go.obj37, go.obj38, go.obj39, go.obj40, go.obj41, go.obj42, go.obj43, go.obj44, go.obj45, go.obj46, go.obj47, go.obj48, go.obj49, go.obj50, go.obj51, go.obj52, go.obj53)

library(RJSONIO)
 exportJSON <- toJSON(AllCombine)
 write(exportJSON,"HighSignificant.out")

 ls -1 *list |paste - <(less HighSignificant.out|sed 's/\[/\t/g' |awk 'NF<15' |grep "\"" | sed '0~5 s/$/\nwait/g' |tr "\n" "\t" |sed 's/wait/\n/g' |sed 's/,/\t/g' |sed 's/"//g' |cut -f 1-7,9- |sed 's/}//g') |less -S
```
### output
```
121EffectorOnlyNewGeneNames.list           notA: 28480            inA: 80          notA: 1170             inA: 41          ]    pval: 7.3351e-27        Jaccard: 0.031758       is.tested: true
AllGenes.list              notA: 28559            inA: 1           notA: 1211             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
AllHGTRevised.list                 notA: 27096            inA: 1464        notA: 1117             inA: 94       odds.ratio: 1.5575      pval: 8.9494e-05        Jaccard: 0.03514        is.tested: true
AllPredictedEffectorLTR.list               notA: 28540            inA: 20          notA: 1209             inA: 2        odds.ratio: 2.3605      pval: 0.22477   Jaccard: 0.0016247      is.tested: true
AllPredictedEffectorTIR.list               notA: 28510            inA: 50          notA: 1203             inA: 8        odds.ratio:  3.792      pval: 0.0023057 Jaccard: 0.0063442      is.tested: true
AllPredictedEffectorTandem.list            notA: 28468            inA: 92          notA: 1187             inA: 24       odds.ratio: 6.2545      pval: 4.1118e-11        Jaccard: 0.018419       is.tested: true
AllPredictedEffectors.list                 notA: 28232            inA: 328         notA: 1108             inA: 103      odds.ratio:  8.001      pval: 1.2287e-49        Jaccard: 0.066927       is.tested: true
BuscoDNA.list              notA: 28536            inA: 24          notA: 1211             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
BuscoLTR.list              notA: 28540            inA: 20          notA: 1211             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
BuscoTandem.list                   notA: 28483            inA: 77          notA: 1209             inA: 2        odds.ratio: 0.61193     pval: 0.83681   Jaccard: 0.0015528      is.tested: true
DorsalEffectorLTR.list             notA: 28555            inA: 5           notA: 1210             inA: 1        odds.ratio: 4.7193      pval: 0.22057   Jaccard: 0.00082237     is.tested: true
DorsalEffectorTIR.list             notA: 28552            inA: 8           notA: 1204             inA: 7        odds.ratio: 20.743      pval: 8.746e-07 Jaccard: 0.0057424      is.tested: true
DorsalEffectorTandem.list                  notA: 28540            inA: 20          notA: 1193             inA: 18       odds.ratio: 21.523      pval: 1.2745e-15        Jaccard: 0.014622       is.tested: true
Dorsal-likeGenesFixed.list                 notA: 28471            inA: 89          notA: 1140             inA: 71       odds.ratio: 19.917      pval: 2.6556e-55        Jaccard: 0.054615       is.tested: true
DuplicatedBuscosGene.list                  notA: 28232            inA: 328         notA: 1196             inA: 15       odds.ratio: 1.0795      pval: 0.42433   Jaccard: 0.0097466      is.tested: true
Family976.list             notA: 28473            inA: 87          notA: 1197             inA: 14       odds.ratio: 3.8278      pval: 5.9472e-05        Jaccard: 0.010786       is.tested: true
Family1265.list            notA: 28451            inA: 109         notA: 1208             inA: 3        odds.ratio: 0.64823     pval: 0.83887   Jaccard: 0.0022727      is.tested: true
HGT.list                   notA: 27581            inA: 979         notA: 1177             inA: 34       odds.ratio: 0.81382     pval: 0.89653   Jaccard: 0.015525       is.tested: true
HGTNewNames.list                   notA: 28418            inA: 142         notA: 1206             inA: 5        odds.ratio: 0.82972     pval: 0.71839   Jaccard: 0.0036955      is.tested: true
HGTRevised.list            notA: 28424            inA: 136         notA: 1196             inA: 15       odds.ratio: 2.6212      pval: 0.0013058 Jaccard: 0.011136       is.tested: true
HighSNPDensity.list                notA: 28093            inA: 467         notA: 1184             inA: 27       odds.ratio: 1.3718      pval: 0.075161  Jaccard: 0.016091       is.tested: true
HighSNPDensity1Dev.list            notA: 26682            inA: 1878        notA: 1111             inA: 100      odds.ratio: 1.2788      pval: 0.014393  Jaccard: 0.032373       is.tested: true
HighSNPDensity2Dev.list            notA: 27728            inA: 832         notA: 1171             inA: 40       odds.ratio: 1.1384      pval: 0.23762   Jaccard: 0.019579       is.tested: true
HighSNPDensity10Perc.list                  notA: 25745            inA: 2815        notA: 1050             inA: 161      odds.ratio: 1.4023      pval: 0.00010255        Jaccard: 0.03999        is.tested: true
HighSignificant.list               notA: 28560            inA: 0           notA: 0                inA: 1211     odds.ratio:   Infinity  pval:      0    Jaccard:      1 is.tested: true
KnownEffectorLTR.list              notA: 28555            inA: 5           notA: 1211             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
KnownEffectorTIR.list              notA: 28548            inA: 12          notA: 1206             inA: 5        odds.ratio: 9.8611      pval: 0.00045379        Jaccard: 0.0040883      is.tested: true
KnownEffectorTandem.list                   notA: 28531            inA: 29          notA: 1202             inA: 9        odds.ratio: 7.3639      pval: 1.66e-05  Jaccard: 0.0072581      is.tested: true
LowSignificant.list                notA: 27992            inA: 568         notA: 1211             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
NonRedundantLtrRetroelementMerge.list              notA: 27170            inA: 1390        notA: 1200             inA: 11       odds.ratio: 0.17918     pval:      1    Jaccard: 0.0042291      is.tested: true
NotAllPredictedEffectorLTR.list            notA: 27190            inA: 1370        notA: 1202             inA: 9        odds.ratio: 0.14861     pval:      1    Jaccard: 0.003487       is.tested: true
NotAllPredictedEffectorTIR.list            notA: 26732            inA: 1828        notA: 1182             inA: 29       odds.ratio: 0.35879     pval:      1    Jaccard: 0.0095426      is.tested: true
NotAllPredictedEffectorTandem.list                 notA: 22487            inA: 6073        notA: 1063             inA: 148      odds.ratio: 0.51556     pval:      1    Jaccard: 0.020319       is.tested: true
NotBuscoDNA.list                   notA: 26706            inA: 1854        notA: 1174             inA: 37       odds.ratio: 0.45398     pval:      1    Jaccard: 0.012072       is.tested: true
NotBuscoLTR.list                   notA: 27190            inA: 1370        notA: 1200             inA: 11       odds.ratio: 0.18193     pval:      1    Jaccard: 0.0042619      is.tested: true
NotBuscoTandem.list                notA: 22472            inA: 6088        notA: 1041             inA: 170      odds.ratio: 0.60281     pval:      1    Jaccard: 0.023291       is.tested: true
NotDorsalEffectorLTR.list                  notA: 27175            inA: 1385        notA: 1201             inA: 10       odds.ratio: 0.16338     pval:      1    Jaccard: 0.0038521      is.tested: true
NotDorsalEffectorTIR.list                  notA: 26690            inA: 1870        notA: 1181             inA: 30       odds.ratio: 0.36257     pval:      1    Jaccard: 0.0097371      is.tested: true
NotDorsalEffectorTandem.list               notA: 22415            inA: 6145        notA: 1057             inA: 154      odds.ratio: 0.53148     pval:      1    Jaccard: 0.020935       is.tested: true
NotKnownEffectorLTR.list                   notA: 27175            inA: 1385        notA: 1200             inA: 11       odds.ratio: 0.17986     pval:      1    Jaccard: 0.0042373      is.tested: true
NotKnownEffectorTIR.list                   notA: 26694            inA: 1866        notA: 1179             inA: 32       odds.ratio: 0.38828     pval:      1    Jaccard: 0.0104 is.tested: true
NotKnownEffectorTandem.list                notA: 22424            inA: 6136        notA: 1048             inA: 163      odds.ratio: 0.56843     pval:      1    Jaccard: 0.022186       is.tested: true
NotPredictedEffectorLTR.list               notA: 27187            inA: 1373        notA: 1202             inA: 9        odds.ratio: 0.14827     pval:      1    Jaccard: 0.003483       is.tested: true
NotPredictedEffectorTIR.list               notA: 26728            inA: 1832        notA: 1180             inA: 31       odds.ratio: 0.38329     pval:      1    Jaccard: 0.010187       is.tested: true
NotPredictedEffectorTandem.list            notA: 22471            inA: 6089        notA: 1053             inA: 158      odds.ratio: 0.55377     pval:      1    Jaccard: 0.021644       is.tested: true
PredictedEffector.list             notA: 28322            inA: 238         notA: 1157             inA: 54       odds.ratio: 5.5533      pval: 7.2954e-21        Jaccard: 0.037267       is.tested: true
PredictedEffectorLTR.list                  notA: 28543            inA: 17          notA: 1209             inA: 2        odds.ratio: 2.7773      pval: 0.17969   Jaccard: 0.0016287      is.tested: true
PredictedEffectorTIR.list                  notA: 28514            inA: 46          notA: 1205             inA: 6        odds.ratio: 3.0863      pval: 0.018495  Jaccard: 0.0047733      is.tested: true
PredictedEffectorTandem.list               notA: 28484            inA: 76          notA: 1197             inA: 14       odds.ratio: 4.3831      pval: 1.5865e-05        Jaccard: 0.010878       is.tested: true
RepeatAffectedGenesNoSimple.list                   notA: 16503            inA: 12057       notA: 911              inA: 300      odds.ratio: 0.45076     pval:      1    Jaccard: 0.022611       is.tested: true
Secretome.gene.list                notA: 26685            inA: 1875        notA: 864              inA: 347      odds.ratio: 5.7144      pval: 6.9306e-115       Jaccard: 0.11244        is.tested: true
SupportedIRFMergeClassified.list                   notA: 26682            inA: 1878        notA: 1174             inA: 37       odds.ratio: 0.44778     pval:      1    Jaccard: 0.011978       is.tested: true
tandem.gene.list                   notA: 22395            inA: 6165        notA: 1039             inA: 172      odds.ratio: 0.60138     pval:      1    Jaccard: 0.023319       is.tested: true
```

#  Low significant expression overlap analyses (10 percentile)
```
go.obj <- newGeneOverlap(x121EffectorOnlyNewGeneNames.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj1 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllGenes.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj2 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllHGTRevised.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj3 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj4 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj5 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj6 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectors.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj7 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoDNA.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj8 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoLTR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj9 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(BuscoTandem.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj10 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj11 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj12 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj13 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsallikeGenesFixed.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj14 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DuplicatedBuscosGene.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj15 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family976.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj16 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family1265.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj17 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGT.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj18 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTNewNames.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj19 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTRevised.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj20 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSNPDensity.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj21 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSNPDensity1Dev.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj22 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSNPDensity2Dev.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj23 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HighSNPDensity10Perc.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj24 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(LowSignificant.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj25 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj26 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj27 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj28 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(LowSignificant.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj29 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NonRedundantLtrRetroelementMerge.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj30 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj31 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj32 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj33 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoDNA.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj34 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoLTR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj35 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotBuscoTandem.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj36 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj37 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj38 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj39 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj40 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj41 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj42 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj43 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj44 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj45 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffector.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj46 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorLTR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj47 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTIR.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj48 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTandem.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj49 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(RepeatAffectedGenesNoSimple.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj50 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Secretome.gene.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj51 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SupportedIRFMergeClassified.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj52 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(tandem.gene.list$V1, LowSignificant.list$V1, genome.size=All.genes)
go.obj53 <- testGeneOverlap(go.obj)
```
Last few commands to get the chart
```
AllCombine <- c(go.obj1, go.obj2, go.obj3, go.obj4, go.obj5, go.obj6, go.obj7, go.obj8, go.obj9, go.obj10, go.obj11, go.obj12, go.obj13, go.obj14, go.obj15, go.obj16, go.obj17, go.obj18, go.obj19, go.obj20, go.obj21, go.obj22, go.obj23, go.obj24, go.obj25, go.obj26, go.obj27, go.obj28, go.obj29, go.obj30, go.obj31, go.obj32, go.obj33, go.obj34, go.obj35, go.obj36, go.obj37, go.obj38, go.obj39, go.obj40, go.obj41, go.obj42, go.obj43, go.obj44, go.obj45, go.obj46, go.obj47, go.obj48, go.obj49, go.obj50, go.obj51, go.obj52, go.obj53)

library(RJSONIO)
 exportJSON <- toJSON(AllCombine)
 write(exportJSON,"LowSignificant.out")

#format the file so it is a readable table
ls -1 *list |paste - <(less LowSignificant.out|sed 's/\[/\t/g' |awk 'NF<15' |grep "\"" | sed '0~5 s/$/\nwait/g' |tr "\n" "\t" |sed 's/wait/\n/g' |sed 's/,/\t/g' |sed 's/"//g' |cut -f 1-7,9- |sed 's/}//g') |less -S
```
### Results from Low significant expression
```
121EffectorOnlyNewGeneNames.list           notA: 29093            inA: 110         notA: 557              inA: 11          ]    pval: 2.1312e-05        Jaccard: 0.016224       is.tested: true
AllGenes.list              notA: 29202            inA: 1           notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
AllHGTRevised.list                 notA: 27678            inA: 1525        notA: 535              inA: 33       odds.ratio: 1.1195      pval: 0.29201   Jaccard: 0.015767       is.tested: true
AllPredictedEffectorLTR.list               notA: 29181            inA: 22          notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
AllPredictedEffectorTIR.list               notA: 29147            inA: 56          notA: 566              inA: 2        odds.ratio: 1.8392      pval: 0.30378   Jaccard: 0.0032051      is.tested: true
AllPredictedEffectorTandem.list            notA: 29095            inA: 108         notA: 560              inA: 8        odds.ratio: 3.8483      pval: 0.001754  Jaccard: 0.011834       is.tested: true
AllPredictedEffectors.list                 notA: 28793            inA: 410         notA: 547              inA: 21       odds.ratio: 2.6958      pval: 9.7619e-05        Jaccard: 0.021472       is.tested: true
BuscoDNA.list              notA: 29179            inA: 24          notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
BuscoLTR.list              notA: 29183            inA: 20          notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
BuscoTandem.list                   notA: 29126            inA: 77          notA: 566              inA: 2        odds.ratio: 1.3366      pval: 0.44644   Jaccard: 0.0031008      is.tested: true
DorsalEffectorLTR.list             notA: 29197            inA: 6           notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
DorsalEffectorTIR.list             notA: 29190            inA: 13          notA: 566              inA: 2        odds.ratio: 7.9331      pval: 0.032376  Jaccard: 0.0034423      is.tested: true
DorsalEffectorTandem.list                  notA: 29170            inA: 33          notA: 563              inA: 5        odds.ratio: 7.8493      pval: 0.00074105        Jaccard: 0.0083195      is.tested: true
Dorsal-likeGenesFixed.list                 notA: 29056            inA: 147         notA: 555              inA: 13       odds.ratio: 4.6295      pval: 1.3198e-05        Jaccard: 0.018182       is.tested: true
DuplicatedBuscosGene.list                  notA: 28868            inA: 335         notA: 560              inA: 8        odds.ratio:  1.231      pval: 0.33264   Jaccard: 0.0088594      is.tested: true
Family976.list             notA: 29105            inA: 98          notA: 565              inA: 3        odds.ratio: 1.5769      pval: 0.30332   Jaccard: 0.0045045      is.tested: true
Family1265.list            notA: 29091            inA: 112         notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
HGT.list                   notA: 28207            inA: 996         notA: 551              inA: 17       odds.ratio: 0.87377     pval: 0.73913   Jaccard: 0.01087        is.tested: true
HGTNewNames.list                   notA: 29059            inA: 144         notA: 565              inA: 3        odds.ratio: 1.0715      pval: 0.53396   Jaccard: 0.0042135      is.tested: true
HGTRevised.list            notA: 29061            inA: 142         notA: 559              inA: 9        odds.ratio: 3.2947      pval: 0.002538  Jaccard: 0.012676       is.tested: true
HighSNPDensity.list                notA: 28717            inA: 486         notA: 560              inA: 8        odds.ratio: 0.84412     pval: 0.7281    Jaccard: 0.0075901      is.tested: true
HighSNPDensity1Dev.list            notA: 27267            inA: 1936        notA: 526              inA: 42       odds.ratio: 1.1246      pval: 0.25659   Jaccard: 0.016773       is.tested: true
HighSNPDensity2Dev.list            notA: 28350            inA: 853         notA: 549              inA: 19       odds.ratio: 1.1502      pval: 0.30939   Jaccard: 0.013371       is.tested: true
HighSNPDensity10Perc.list                  notA: 26290            inA: 2913        notA: 505              inA: 63       odds.ratio: 1.1259      pval: 0.20769   Jaccard: 0.018098       is.tested: true
HighSignificant.list               notA: 29203            inA: 0           notA: 0                inA: 568      odds.ratio:   Infinity  pval:      0    Jaccard:      1 is.tested: true
KnownEffectorLTR.list              notA: 29198            inA: 5           notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
KnownEffectorTIR.list              notA: 29187            inA: 16          notA: 567              inA: 1        odds.ratio:  3.217      pval: 0.27932   Jaccard: 0.0017123      is.tested: true
KnownEffectorTandem.list                   notA: 29170            inA: 33          notA: 563              inA: 5        odds.ratio: 7.8493      pval: 0.00074105        Jaccard: 0.0083195      is.tested: true
LowSignificant.list                notA: 29203            inA: 0           notA: 0                inA: 568      odds.ratio:   Infinity  pval:      0    Jaccard:      1 is.tested: true
NonRedundantLtrRetroelementMerge.list              notA: 27809            inA: 1394        notA: 561              inA: 7        odds.ratio: 0.24893     pval:      1    Jaccard: 0.0035678      is.tested: true
NotAllPredictedEffectorLTR.list            notA: 27831            inA: 1372        notA: 561              inA: 7        odds.ratio: 0.25312     pval:      1    Jaccard: 0.0036082      is.tested: true
NotAllPredictedEffectorTIR.list            notA: 27362            inA: 1841        notA: 552              inA: 16       odds.ratio: 0.43081     pval: 0.99995   Jaccard: 0.0066418      is.tested: true
NotAllPredictedEffectorTandem.list                 notA: 23025            inA: 6178        notA: 525              inA: 43       odds.ratio: 0.30526     pval:      1    Jaccard: 0.0063741      is.tested: true
NotBuscoDNA.list                   notA: 27330            inA: 1873        notA: 550              inA: 18       odds.ratio: 0.47755     pval: 0.99979   Jaccard: 0.007374       is.tested: true
NotBuscoLTR.list                   notA: 27829            inA: 1374        notA: 561              inA: 7        odds.ratio: 0.25273     pval:      1    Jaccard: 0.0036045      is.tested: true
NotBuscoTandem.list                notA: 22994            inA: 6209        notA: 519              inA: 49       odds.ratio: 0.34964     pval:      1    Jaccard: 0.0072303      is.tested: true
NotDorsalEffectorLTR.list                  notA: 27815            inA: 1388        notA: 561              inA: 7        odds.ratio: 0.25006     pval:      1    Jaccard: 0.0035787      is.tested: true
NotDorsalEffectorTIR.list                  notA: 27319            inA: 1884        notA: 552              inA: 16       odds.ratio: 0.42031     pval: 0.99997   Jaccard: 0.0065253      is.tested: true
NotDorsalEffectorTandem.list               notA: 22950            inA: 6253        notA: 522              inA: 46       odds.ratio: 0.32343     pval:      1    Jaccard: 0.0067439      is.tested: true
NotKnownEffectorLTR.list                   notA: 27814            inA: 1389        notA: 561              inA: 7        odds.ratio: 0.24987     pval:      1    Jaccard: 0.0035769      is.tested: true
NotKnownEffectorTIR.list                   notA: 27322            inA: 1881        notA: 551              inA: 17       odds.ratio: 0.44816     pval: 0.99992   Jaccard: 0.0069416      is.tested: true
NotKnownEffectorTandem.list                notA: 22950            inA: 6253        notA: 522              inA: 46       odds.ratio: 0.32343     pval:      1    Jaccard: 0.0067439      is.tested: true
NotPredictedEffectorLTR.list               notA: 27828            inA: 1375        notA: 561              inA: 7        odds.ratio: 0.25254     pval:      1    Jaccard: 0.0036027      is.tested: true
NotPredictedEffectorTIR.list               notA: 27357            inA: 1846        notA: 551              inA: 17       odds.ratio: 0.45724     pval: 0.99988   Jaccard: 0.0070423      is.tested: true
NotPredictedEffectorTandem.list            notA: 23001            inA: 6202        notA: 523              inA: 45       odds.ratio: 0.3191      pval:      1    Jaccard: 0.006647       is.tested: true
PredictedEffector.list             notA: 28918            inA: 285         notA: 561              inA: 7        odds.ratio: 1.2661      pval: 0.32415   Jaccard: 0.0082063      is.tested: true
PredictedEffectorLTR.list                  notA: 29184            inA: 19          notA: 568              inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
PredictedEffectorTIR.list                  notA: 29152            inA: 51          notA: 567              inA: 1        odds.ratio: 1.0081      pval: 0.63306   Jaccard: 0.0016155      is.tested: true
PredictedEffectorTandem.list               notA: 29119            inA: 84          notA: 562              inA: 6        odds.ratio: 3.7009      pval: 0.0075792 Jaccard: 0.0092025      is.tested: true
RepeatAffectedGenesNoSimple.list                   notA: 16989            inA: 12214       notA: 425              inA: 143      odds.ratio: 0.46803     pval:      1    Jaccard: 0.011188       is.tested: true
Secretome.gene.list                notA: 27081            inA: 2122        notA: 468              inA: 100      odds.ratio: 2.7267      pval: 6.4594e-16        Jaccard: 0.037175       is.tested: true
SupportedIRFMergeClassified.list                   notA: 27306            inA: 1897        notA: 550              inA: 18       odds.ratio: 0.4711      pval: 0.99984   Jaccard: 0.0073022      is.tested: true
tandem.gene.list                   notA: 22917            inA: 6286        notA: 517              inA: 51       odds.ratio: 0.35964     pval:      1    Jaccard: 0.0074409      is.tested: true
```

# Gene overlap with low SNP density (10 percentile, no 0's)
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/28_SNPDensity
sed 's/\./\t/2' ../30_BenNewHGT/AllHGT.list |awk '{print $1}' >AllHGTRevised.list
sed 's/\./\t/2' ../30_BenNewHGT/HighConfHGT.list |awk '{print $1}' >HGTRevised.list

ls -1 *list|awk '{print $1" <- read.table(\""$1"\")"}' |less

###############################################
x121EffectorOnlyNewGeneNames.list <- read.table("121EffectorOnlyNewGeneNames.list")
AllGenes.list <- read.table("AllGenes.list")
AllHGTRevised.list <- read.table("AllHGTRevised.list")
AllPredictedEffectorLTR.list <- read.table("AllPredictedEffectorLTR.list")
AllPredictedEffectorTIR.list <- read.table("AllPredictedEffectorTIR.list")
AllPredictedEffectorTandem.list <- read.table("AllPredictedEffectorTandem.list")
AllPredictedEffectors.list <- read.table("AllPredictedEffectors.list")
DorsalEffectorLTR.list <- read.table("DorsalEffectorLTR.list")
DorsalEffectorTIR.list <- read.table("DorsalEffectorTIR.list")
DorsalEffectorTandem.list <- read.table("DorsalEffectorTandem.list")
DorsallikeGenesFixed.list <- read.table("Dorsal-likeGenesFixed.list")
Family976.list <- read.table("Family976.list")
Family1265.list <- read.table("Family1265.list")
HGT.list <- read.table("HGT.list")
HGTNewNames.list <- read.table("HGTNewNames.list")
HGTRevised.list <- read.table("HGTRevised.list")
KnownEffectorLTR.list <- read.table("KnownEffectorLTR.list")
KnownEffectorTIR.list <- read.table("KnownEffectorTIR.list")
KnownEffectorTandem.list <- read.table("KnownEffectorTandem.list")
NonRedundantLtrRetroelementMerge.list <- read.table("NonRedundantLtrRetroelementMerge.list")
NotAllPredictedEffectorLTR.list <- read.table("NotAllPredictedEffectorLTR.list")
NotAllPredictedEffectorTIR.list <- read.table("NotAllPredictedEffectorTIR.list")
NotAllPredictedEffectorTandem.list <- read.table("NotAllPredictedEffectorTandem.list")
NotDorsalEffectorLTR.list <- read.table("NotDorsalEffectorLTR.list")
NotDorsalEffectorTIR.list <- read.table("NotDorsalEffectorTIR.list")
NotDorsalEffectorTandem.list <- read.table("NotDorsalEffectorTandem.list")
NotKnownEffectorLTR.list <- read.table("NotKnownEffectorLTR.list")
NotKnownEffectorTIR.list <- read.table("NotKnownEffectorTIR.list")
NotKnownEffectorTandem.list <- read.table("NotKnownEffectorTandem.list")
NotPredictedEffectorLTR.list <- read.table("NotPredictedEffectorLTR.list")
NotPredictedEffectorTIR.list <- read.table("NotPredictedEffectorTIR.list")
NotPredictedEffectorTandem.list <- read.table("NotPredictedEffectorTandem.list")
PredictedEffector.list <- read.table("PredictedEffector.list")
PredictedEffectorLTR.list <- read.table("PredictedEffectorLTR.list")
PredictedEffectorTIR.list <- read.table("PredictedEffectorTIR.list")
PredictedEffectorTandem.list <- read.table("PredictedEffectorTandem.list")
RepeatAffectedGenesNoSimple.list <- read.table("RepeatAffectedGenesNoSimple.list")
SNPDensityAll0SNPS.list <- read.table("SNPDensityAll0SNPS.list")
SNPDensityAll10Percentile.list <- read.table("SNPDensityAll10Percentile.list")
SNPDensityAll10percAfterRemoveZeros.list <- read.table("SNPDensityAll10percAfterRemoveZeros.list")
SNPDensityAll90Percentile.list <- read.table("SNPDensityAll90Percentile.list")
Secretome.gene.list <- read.table("Secretome.gene.list")
SupportedIRFMergeClassified.list <- read.table("SupportedIRFMergeClassified.list")
tandem.gene.list <- read.table("tandem.gene.list")
###############################################

ls -1 *list |awk '{print "go.obj <- newGeneOverlap("$1"$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)"}' |awk '{print $0"\ngo.obj"NR" <- testGeneOverlap(go.obj)"}' |less
#################################################
go.obj <- newGeneOverlap(x121EffectorOnlyNewGeneNames.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj1 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllGenes.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj2 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllHGTRevised.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj3 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj4 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj5 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj6 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectors.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj7 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj8 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj9 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj10 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsallikeGenesFixed.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj11 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family976.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj12 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family1265.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj13 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGT.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj14 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTNewNames.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj15 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTRevised.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj16 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj17 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj18 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj19 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NonRedundantLtrRetroelementMerge.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj20 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj21 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj22 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj23 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj24 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj25 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj26 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj27 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj28 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj29 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj30 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj31 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj32 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffector.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj33 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorLTR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj34 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTIR.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj35 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTandem.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj36 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(RepeatAffectedGenesNoSimple.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj37 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll0SNPS.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj38 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll10Percentile.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj39 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll10percAfterRemoveZeros.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj40 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll90Percentile.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj41 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Secretome.gene.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj42 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SupportedIRFMergeClassified.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj43 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(tandem.gene.list$V1, SNPDensityAll10percAfterRemoveZeros.list$V1, genome.size=All.genes)
go.obj44 <- testGeneOverlap(go.obj)
#################################################

#last few commands to perform the tested
AllCombine <- c(go.obj1, go.obj2, go.obj3, go.obj4, go.obj5, go.obj6, go.obj7, go.obj8, go.obj9, go.obj10, go.obj11, go.obj12, go.obj13, go.obj14, go.obj15, go.obj16, go.obj17, go.obj18, go.obj19, go.obj20, go.obj21, go.obj22, go.obj23, go.obj24, go.obj25, go.obj26, go.obj27, go.obj28, go.obj29, go.obj30, go.obj31, go.obj32, go.obj33, go.obj34, go.obj35, go.obj36, go.obj37, go.obj38, go.obj39, go.obj40, go.obj41, go.obj42, go.obj43, go.obj44)
library(RJSONIO)
exportJSON <- toJSON(AllCombine)
write(exportJSON,"SNPDensityAll10percAfterRemoveZeros.out")
ls -1 *list |paste - <(less SNPDensityAll10percAfterRemoveZeros.out|sed 's/\[/\t/g' |awk 'NF<15' |grep "\"" | sed '0~5 s/$/\nwait/g' |tr "\n" "\t" |sed 's/wait/\n/g' |sed 's/,/\t/g' |sed 's/"//g' |cut -f 1-7,9- |sed 's/}//g') |less -S
```
### Results of overlap of low snp density
```
121EffectorOnlyNewGeneNames.list           notA: 27138            inA: 117         notA: 2512             inA: 4           ]    pval: 0.99319   Jaccard: 0.0015192      is.tested: true
AllGenes.list              notA: 27254            inA: 1           notA: 2516             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
AllHGTRevised.list                 notA: 25866            inA: 1389        notA: 2347             inA: 169      odds.ratio: 1.3409      pval: 0.00042846        Jaccard: 0.043278       is.tested: true
AllPredictedEffectorLTR.list               notA: 27238            inA: 17          notA: 2511             inA: 5        odds.ratio: 3.1902      pval: 0.033545  Jaccard: 0.0019739      is.tested: true
AllPredictedEffectorTIR.list               notA: 27205            inA: 50          notA: 2508             inA: 8        odds.ratio: 1.7355      pval: 0.11367   Jaccard: 0.0031177      is.tested: true
AllPredictedEffectorTandem.list            notA: 27151            inA: 104         notA: 2504             inA: 12       odds.ratio: 1.2511      pval: 0.27473   Jaccard: 0.0045802      is.tested: true
AllPredictedEffectors.list                 notA: 26859            inA: 396         notA: 2481             inA: 35       odds.ratio: 0.95683     pval: 0.62347   Jaccard: 0.012019       is.tested: true
DorsalEffectorLTR.list             notA: 27251            inA: 4           notA: 2514             inA: 2        odds.ratio: 5.4188      pval: 0.085166  Jaccard: 0.00079365     is.tested: true
DorsalEffectorTIR.list             notA: 27240            inA: 15          notA: 2516             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
DorsalEffectorTandem.list                  notA: 27221            inA: 34          notA: 2512             inA: 4        odds.ratio: 1.2749      pval: 0.40203   Jaccard: 0.0015686      is.tested: true
Dorsal-likeGenesFixed.list                 notA: 27107            inA: 148         notA: 2504             inA: 12       odds.ratio: 0.87774     pval: 0.70861   Jaccard: 0.0045045      is.tested: true
Family976.list             notA: 27160            inA: 95          notA: 2510             inA: 6        odds.ratio: 0.68345     pval: 0.86488   Jaccard: 0.002298       is.tested: true
Family1265.list            notA: 27152            inA: 103         notA: 2507             inA: 9        odds.ratio: 0.94635     pval: 0.6123    Jaccard: 0.0034364      is.tested: true
HGT.list                   notA: 26319            inA: 936         notA: 2439             inA: 77       odds.ratio: 0.8877      pval: 0.85289   Jaccard: 0.022306       is.tested: true
HGTNewNames.list                   notA: 27120            inA: 135         notA: 2504             inA: 12       odds.ratio: 0.96273     pval: 0.59319   Jaccard: 0.0045266      is.tested: true
HGTRevised.list            notA: 27120            inA: 135         notA: 2500             inA: 16       odds.ratio: 1.2857      pval: 0.20657   Jaccard: 0.0060355      is.tested: true
KnownEffectorLTR.list              notA: 27250            inA: 5           notA: 2516             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
KnownEffectorTIR.list              notA: 27240            inA: 15          notA: 2514             inA: 2        odds.ratio: 1.4447      pval: 0.42735   Jaccard: 0.0007902      is.tested: true
KnownEffectorTandem.list                   notA: 27220            inA: 35          notA: 2513             inA: 3        odds.ratio: 0.92843     pval: 0.63377   Jaccard: 0.001176       is.tested: true
NonRedundantLtrRetroelementMerge.list              notA: 25977            inA: 1278        notA: 2393             inA: 123      odds.ratio: 1.0448      pval: 0.33958   Jaccard: 0.03242        is.tested: true
NotAllPredictedEffectorLTR.list            notA: 25994            inA: 1261        notA: 2398             inA: 118      odds.ratio: 1.0144      pval: 0.45726   Jaccard: 0.031242       is.tested: true
NotAllPredictedEffectorTIR.list            notA: 25568            inA: 1687        notA: 2346             inA: 170      odds.ratio: 1.0982      pval: 0.13991   Jaccard: 0.040447       is.tested: true
NotAllPredictedEffectorTandem.list                 notA: 21618            inA: 5637        notA: 1932             inA: 584      odds.ratio: 1.1592      pval: 0.001698  Jaccard: 0.07163        is.tested: true
NotDorsalEffectorLTR.list                  notA: 25981            inA: 1274        notA: 2395             inA: 121      odds.ratio: 1.0303      pval: 0.39417   Jaccard: 0.031926       is.tested: true
NotDorsalEffectorTIR.list                  notA: 25533            inA: 1722        notA: 2338             inA: 178      odds.ratio: 1.1288      pval: 0.075997  Jaccard: 0.042001       is.tested: true
NotDorsalEffectorTandem.list               notA: 21548            inA: 5707        notA: 1924             inA: 592      odds.ratio: 1.1617      pval: 0.0014094 Jaccard: 0.071993       is.tested: true
NotKnownEffectorLTR.list                   notA: 25982            inA: 1273        notA: 2393             inA: 123      odds.ratio: 1.0491      pval: 0.32437   Jaccard: 0.032462       is.tested: true
NotKnownEffectorTIR.list                   notA: 25533            inA: 1722        notA: 2340             inA: 176      odds.ratio: 1.1152      pval: 0.10006   Jaccard: 0.041529       is.tested: true
NotKnownEffectorTandem.list                notA: 21549            inA: 5706        notA: 1923             inA: 593      odds.ratio: 1.1645      pval: 0.001196  Jaccard: 0.072124       is.tested: true
NotPredictedEffectorLTR.list               notA: 25993            inA: 1262        notA: 2396             inA: 120      odds.ratio: 1.0316      pval: 0.38998   Jaccard: 0.031763       is.tested: true
NotPredictedEffectorTIR.list               notA: 25562            inA: 1693        notA: 2346             inA: 170      odds.ratio: 1.0941      pval:   0.15    Jaccard: 0.04039        is.tested: true
NotPredictedEffectorTandem.list            notA: 21597            inA: 5658        notA: 1927             inA: 589      odds.ratio: 1.1667      pval: 0.0010875 Jaccard: 0.072058       is.tested: true
PredictedEffector.list             notA: 26985            inA: 270         notA: 2494             inA: 22       odds.ratio: 0.88162     pval: 0.74404   Jaccard: 0.0078966      is.tested: true
PredictedEffectorLTR.list                  notA: 27239            inA: 16          notA: 2513             inA: 3        odds.ratio: 2.0323      pval: 0.21327   Jaccard: 0.0011848      is.tested: true
PredictedEffectorTIR.list                  notA: 27211            inA: 44          notA: 2508             inA: 8        odds.ratio: 1.9726      pval: 0.069072  Jaccard: 0.003125       is.tested: true
PredictedEffectorTandem.list               notA: 27172            inA: 83          notA: 2509             inA: 7        odds.ratio: 0.91336     pval: 0.64605   Jaccard: 0.0026933      is.tested: true
RepeatAffectedGenesNoSimple.list                   notA: 15733            inA: 11522       notA: 1681             inA: 835      odds.ratio: 0.67828     pval:      1    Jaccard: 0.059481       is.tested: true
SNPDensityAll0SNPS.list            notA: 22643            inA: 4612        notA: 2516             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
SNPDensityAll10Percentile.list             notA: 24285            inA: 2970        notA: 2516             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
SNPDensityAll10percAfterRemoveZeros.list                   notA: 27255            inA: 0           notA: 0                inA: 2516     odds.ratio:   Infinity  pval:      0    Jaccard:      1 is.tested: true
SNPDensityAll90Percentile.list             notA: 24286            inA: 2969        notA: 2516             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
Secretome.gene.list                notA: 25189            inA: 2066        notA: 2360             inA: 156      odds.ratio: 0.80593     pval: 0.99556   Jaccard: 0.034046       is.tested: true
SupportedIRFMergeClassified.list                   notA: 25518            inA: 1737        notA: 2338             inA: 178      odds.ratio: 1.1184      pval: 0.092974  Jaccard: 0.041853       is.tested: true
tandem.gene.list                   notA: 21514            inA: 5741        notA: 1920             inA: 596      odds.ratio: 1.1632      pval: 0.0012647 Jaccard: 0.072181       is.tested: true
```
# Gene overlap for high snp density (90 percentile)
```
ls -1 *list |awk '{print "go.obj <- newGeneOverlap("$1"$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)"}' |awk '{print $0"\ngo.obj"NR" <- testGeneOverlap(go.obj)"}' |less

go.obj <- newGeneOverlap(x121EffectorOnlyNewGeneNames.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj1 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllGenes.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj2 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllHGTRevised.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj3 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj4 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj5 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj6 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectors.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj7 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj8 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj9 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj10 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsallikeGenesFixed.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj11 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family976.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj12 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family1265.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj13 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGT.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj14 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTNewNames.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj15 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTRevised.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj16 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj17 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj18 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj19 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NonRedundantLtrRetroelementMerge.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj20 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj21 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj22 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj23 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj24 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj25 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj26 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj27 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj28 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj29 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj30 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj31 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj32 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffector.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj33 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorLTR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj34 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTIR.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj35 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTandem.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj36 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(RepeatAffectedGenesNoSimple.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj37 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll0SNPS.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj38 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll10Percentile.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj39 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll10percAfterRemoveZeros.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj40 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll90Percentile.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj41 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Secretome.gene.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj42 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SupportedIRFMergeClassified.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj43 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(tandem.gene.list$V1, SNPDensityAll90Percentile.list$V1, genome.size=All.genes)
go.obj44 <- testGeneOverlap(go.obj)
###############################################

#commands that perform the test
AllCombine <- c(go.obj1, go.obj2, go.obj3, go.obj4, go.obj5, go.obj6, go.obj7, go.obj8, go.obj9, go.obj10, go.obj11, go.obj12, go.obj13, go.obj14, go.obj15, go.obj16, go.obj17, go.obj18, go.obj19, go.obj20, go.obj21, go.obj22, go.obj23, go.obj24, go.obj25, go.obj26, go.obj27, go.obj28, go.obj29, go.obj30, go.obj31, go.obj32, go.obj33, go.obj34, go.obj35, go.obj36, go.obj37, go.obj38, go.obj39, go.obj40, go.obj41, go.obj42, go.obj43, go.obj44)
library(RJSONIO)
exportJSON <- toJSON(AllCombine)
write(exportJSON,"SNPDensityAll90Percentile.out")

#take the results and reformat into a table.
ls -1 *list |paste - <(less SNPDensityAll90Percentile.out|sed 's/\[/\t/g' |awk 'NF<15' |grep "\"" | sed '0~5 s/$/\nwait/g' |tr "\n" "\t" |sed 's/wait/\n/g' |sed 's/,/\t/g' |sed 's/"//g' |cut -f 1-7,9- |sed 's/}//g') |less -S
```
### Results from High snp density overlap analysis
```
121EffectorOnlyNewGeneNames.list           notA: 26698            inA: 104         notA: 2952             inA: 17          ]    pval: 0.093069  Jaccard: 0.0055321      is.tested: true
AllGenes.list              notA: 26801            inA: 1           notA: 2969             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
AllHGTRevised.list                 notA: 25389            inA: 1413        notA: 2824             inA: 145      odds.ratio: 0.92258     pval: 0.82732   Jaccard: 0.03309        is.tested: true
AllPredictedEffectorLTR.list               notA: 26781            inA: 21          notA: 2968             inA: 1        odds.ratio: 0.42969     pval: 0.90095   Jaccard: 0.00033445     is.tested: true
AllPredictedEffectorTIR.list               notA: 26751            inA: 51          notA: 2962             inA: 7        odds.ratio: 1.2396      pval: 0.35702   Jaccard: 0.0023179      is.tested: true
AllPredictedEffectorTandem.list            notA: 26700            inA: 102         notA: 2955             inA: 14       odds.ratio: 1.2402      pval: 0.26566   Jaccard: 0.0045588      is.tested: true
AllPredictedEffectors.list                 notA: 26411            inA: 391         notA: 2929             inA: 40       odds.ratio: 0.92246     pval: 0.70865   Jaccard: 0.011905       is.tested: true
DorsalEffectorLTR.list             notA: 26796            inA: 6           notA: 2969             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
DorsalEffectorTIR.list             notA: 26790            inA: 12          notA: 2966             inA: 3        odds.ratio:  2.258      pval: 0.18297   Jaccard: 0.0010064      is.tested: true
DorsalEffectorTandem.list                  notA: 26771            inA: 31          notA: 2962             inA: 7        odds.ratio: 2.0408      pval: 0.078911  Jaccard: 0.0023333      is.tested: true
Dorsal-likeGenesFixed.list                 notA: 26660            inA: 142         notA: 2951             inA: 18       odds.ratio: 1.1452      pval: 0.33097   Jaccard: 0.0057859      is.tested: true
Family976.list             notA: 26711            inA: 91          notA: 2959             inA: 10       odds.ratio: 0.99198     pval: 0.55838   Jaccard: 0.003268       is.tested: true
Family1265.list            notA: 26698            inA: 104         notA: 2961             inA: 8        odds.ratio: 0.6936      pval: 0.88089   Jaccard: 0.0026033      is.tested: true
HGT.list                   notA: 25898            inA: 904         notA: 2860             inA: 109      odds.ratio: 1.0919      pval: 0.21114   Jaccard: 0.028144       is.tested: true
HGTNewNames.list                   notA: 26673            inA: 129         notA: 2951             inA: 18       odds.ratio: 1.2612      pval: 0.21235   Jaccard: 0.0058102      is.tested: true
HGTRevised.list            notA: 26655            inA: 147         notA: 2965             inA: 4        odds.ratio: 0.24463     pval: 0.99988   Jaccard: 0.0012837      is.tested: true
KnownEffectorLTR.list              notA: 26797            inA: 5           notA: 2969             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
KnownEffectorTIR.list              notA: 26788            inA: 14          notA: 2966             inA: 3        odds.ratio: 1.9354      pval: 0.2369    Jaccard: 0.0010057      is.tested: true
KnownEffectorTandem.list                   notA: 26770            inA: 32          notA: 2963             inA: 6        odds.ratio:  1.694      pval: 0.17308   Jaccard: 0.0019993      is.tested: true
NonRedundantLtrRetroelementMerge.list              notA: 25661            inA: 1141        notA: 2709             inA: 260      odds.ratio: 2.1584      pval: 1.0265e-23        Jaccard: 0.06326        is.tested: true
NotAllPredictedEffectorLTR.list            notA: 25682            inA: 1120        notA: 2710             inA: 259      odds.ratio: 2.1916      pval: 2.1326e-24        Jaccard: 0.063341       is.tested: true
NotAllPredictedEffectorTIR.list            notA: 25209            inA: 1593        notA: 2705             inA: 264      odds.ratio: 1.5444      pval: 1.2563e-09        Jaccard: 0.057869       is.tested: true
NotAllPredictedEffectorTandem.list                 notA: 21396            inA: 5406        notA: 2154             inA: 815      odds.ratio: 1.4975      pval: 1.5573e-19        Jaccard: 0.097313       is.tested: true
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
PredictedEffectorTandem.list               notA: 26725            inA: 77          notA: 2956             inA: 13       odds.ratio: 1.5264      pval: 0.11055   Jaccard: 0.0042679      is.tested: true
RepeatAffectedGenesNoSimple.list                   notA: 16334            inA: 10468       notA: 1080             inA: 1889     odds.ratio: 2.7291      pval: 2.0423e-144       Jaccard: 0.14058        is.tested: true
SNPDensityAll0SNPS.list            notA: 22190            inA: 4612        notA: 2969             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
SNPDensityAll10Percentile.list             notA: 23832            inA: 2970        notA: 2969             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
SNPDensityAll10percAfterRemoveZeros.list                   notA: 24286            inA: 2516        notA: 2969             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
SNPDensityAll90Percentile.list             notA: 26802            inA: 0           notA: 0                inA: 2969     odds.ratio:   Infinity  pval:      0    Jaccard:      1 is.tested: true
Secretome.gene.list                notA: 24756            inA: 2046        notA: 2793             inA: 176      odds.ratio: 0.76247     pval: 0.99976   Jaccard: 0.035095       is.tested: true
SupportedIRFMergeClassified.list                   notA: 25158            inA: 1644        notA: 2698             inA: 271      odds.ratio: 1.5371      pval: 1.1741e-09        Jaccard: 0.058747       is.tested: true
tandem.gene.list                   notA: 21294            inA: 5508        notA: 2140             inA: 829      odds.ratio: 1.4976      pval: 9.439e-20 Jaccard: 0.097794       is.tested: true
```
# Gene overlap analyses for zero SNPs
```
ls -1 *list |awk '{print "go.obj <- newGeneOverlap("$1"$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)"}' |awk '{print $0"\ngo.obj"NR" <- testGeneOverlap(go.obj)"}' |less
####################################################

go.obj <- newGeneOverlap(x121EffectorOnlyNewGeneNames.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj1 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllGenes.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj2 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllHGTRevised.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj3 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj4 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj5 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj6 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(AllPredictedEffectors.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj7 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj8 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj9 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsalEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj10 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(DorsallikeGenesFixed.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj11 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family976.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj12 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Family1265.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj13 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGT.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj14 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTNewNames.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj15 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(HGTRevised.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj16 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj17 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj18 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(KnownEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj19 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NonRedundantLtrRetroelementMerge.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj20 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj21 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj22 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotAllPredictedEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj23 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj24 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj25 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotDorsalEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj26 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj27 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj28 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotKnownEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj29 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj30 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj31 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(NotPredictedEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj32 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffector.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj33 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorLTR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj34 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTIR.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj35 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(PredictedEffectorTandem.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj36 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(RepeatAffectedGenesNoSimple.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj37 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll0SNPS.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj38 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll10Percentile.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj39 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll10percAfterRemoveZeros.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj40 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SNPDensityAll90Percentile.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj41 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(Secretome.gene.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj42 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(SupportedIRFMergeClassified.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj43 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(tandem.gene.list$V1, SNPDensityAll0SNPS.list$V1, genome.size=All.genes)
go.obj44 <- testGeneOverlap(go.obj)
####################################################

#generate the commands to perform the analysis
AllCombine <- c(go.obj1, go.obj2, go.obj3, go.obj4, go.obj5, go.obj6, go.obj7, go.obj8, go.obj9, go.obj10, go.obj11, go.obj12, go.obj13, go.obj14, go.obj15, go.obj16, go.obj17, go.obj18, go.obj19, go.obj20, go.obj21, go.obj22, go.obj23, go.obj24, go.obj25, go.obj26, go.obj27, go.obj28, go.obj29, go.obj30, go.obj31, go.obj32, go.obj33, go.obj34, go.obj35, go.obj36, go.obj37, go.obj38, go.obj39, go.obj40, go.obj41, go.obj42, go.obj43, go.obj44)
library(RJSONIO)
exportJSON <- toJSON(AllCombine)
write(exportJSON,"SNPDensityAll0SNPS.list.out")
ls -1 *list |paste - <(less SNPDensityAll0SNPS.list.out|sed 's/\[/\t/g' |awk 'NF<15' |grep "\"" | sed '0~5 s/$/\nwait/g' |tr "\n" "\t" |sed 's/wait/\n/g' |sed 's/,/\t/g' |sed 's/"//g' |cut -f 1-7,9- |sed 's/}//g') |less -S
```
### results from zero snp analysis
```
121EffectorOnlyNewGeneNames.list           notA: 25052            inA: 107         notA: 4598             inA: 14          ]    pval: 0.91097   Jaccard: 0.0029667      is.tested: true
AllGenes.list              notA: 25158            inA: 1           notA: 4612             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
AllHGTRevised.list                 notA: 23861            inA: 1298        notA: 4352             inA: 260      odds.ratio: 1.0983      pval: 0.096817  Jaccard: 0.043993       is.tested: true
AllPredictedEffectorLTR.list               notA: 25145            inA: 14          notA: 4604             inA: 8        odds.ratio: 3.1207      pval: 0.01373   Jaccard: 0.0017294      is.tested: true
AllPredictedEffectorTIR.list               notA: 25112            inA: 47          notA: 4601             inA: 11       odds.ratio: 1.2774      pval: 0.28128   Jaccard: 0.002361       is.tested: true
AllPredictedEffectorTandem.list            notA: 25070            inA: 89          notA: 4585             inA: 27       odds.ratio: 1.6587      pval: 0.017746  Jaccard: 0.0057435      is.tested: true
AllPredictedEffectors.list                 notA: 24806            inA: 353         notA: 4534             inA: 78       odds.ratio: 1.2089      pval: 0.077142  Jaccard: 0.01571        is.tested: true
DorsalEffectorLTR.list             notA: 25153            inA: 6           notA: 4612             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
DorsalEffectorTIR.list             notA: 25144            inA: 15          notA: 4612             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
DorsalEffectorTandem.list                  notA: 25123            inA: 36          notA: 4610             inA: 2        odds.ratio: 0.30277     pval: 0.98675   Jaccard: 0.00043029     is.tested: true
Dorsal-likeGenesFixed.list                 notA: 25007            inA: 152         notA: 4604             inA: 8        odds.ratio: 0.28588     pval: 0.99999   Jaccard: 0.0016793      is.tested: true
Family976.list             notA: 25082            inA: 77          notA: 4588             inA: 24       odds.ratio: 1.7039      pval: 0.019185  Jaccard: 0.0051184      is.tested: true
Family1265.list            notA: 25081            inA: 78          notA: 4578             inA: 34       odds.ratio: 2.3881      pval: 5.4716e-05        Jaccard: 0.0072495      is.tested: true
HGT.list                   notA: 24326            inA: 833         notA: 4432             inA: 180      odds.ratio:  1.186      pval: 0.02454   Jaccard: 0.033058       is.tested: true
HGTNewNames.list                   notA: 25035            inA: 124         notA: 4589             inA: 23       odds.ratio: 1.0119      pval: 0.51441   Jaccard: 0.0048564      is.tested: true
HGTRevised.list            notA: 25034            inA: 125         notA: 4586             inA: 26       odds.ratio: 1.1354      pval: 0.31032   Jaccard: 0.0054887      is.tested: true
KnownEffectorLTR.list              notA: 25156            inA: 3           notA: 4610             inA: 2        odds.ratio: 3.6377      pval: 0.1739    Jaccard: 0.00043337     is.tested: true
KnownEffectorTIR.list              notA: 25144            inA: 15          notA: 4610             inA: 2        odds.ratio: 0.72724     pval: 0.76468   Jaccard: 0.00043225     is.tested: true
KnownEffectorTandem.list                   notA: 25129            inA: 30          notA: 4604             inA: 8        odds.ratio: 1.4555      pval: 0.22699   Jaccard: 0.0017234      is.tested: true
NonRedundantLtrRetroelementMerge.list              notA: 24055            inA: 1104        notA: 4315             inA: 297      odds.ratio: 1.4997      pval: 3.8962e-09        Jaccard: 0.051959       is.tested: true
NotAllPredictedEffectorLTR.list            notA: 24069            inA: 1090        notA: 4323             inA: 289      odds.ratio: 1.4762      pval: 2.0179e-08        Jaccard: 0.050684       is.tested: true
NotAllPredictedEffectorTIR.list            notA: 23752            inA: 1407        notA: 4162             inA: 450      odds.ratio: 1.8251      pval: 2.4807e-24        Jaccard: 0.074763       is.tested: true
NotAllPredictedEffectorTandem.list                 notA: 20280            inA: 4879        notA: 3270             inA: 1342     odds.ratio: 1.7058      pval: 2.5775e-47        Jaccard: 0.1414 is.tested: true
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
PredictedEffectorTandem.list               notA: 25095            inA: 64          notA: 4586             inA: 26       odds.ratio:  2.223      pval: 0.00090271        Jaccard: 0.0055603      is.tested: true
RepeatAffectedGenesNoSimple.list                   notA: 14920            inA: 10239       notA: 2494             inA: 2118     odds.ratio: 1.2375      pval: 2.3104e-11        Jaccard: 0.14262        is.tested: true
SNPDensityAll0SNPS.list            notA: 25159            inA: 0           notA: 0                inA: 4612     odds.ratio:   Infinity  pval:      0    Jaccard:      1 is.tested: true
SNPDensityAll10Percentile.list             notA: 25159            inA: 0           notA: 1642             inA: 2970     odds.ratio:   Infinity  pval:      0    Jaccard: 0.64397        is.tested: true
SNPDensityAll10percAfterRemoveZeros.list                   notA: 22643            inA: 2516        notA: 4612             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
SNPDensityAll90Percentile.list             notA: 22190            inA: 2969        notA: 4612             inA: 0        odds.ratio:      0      pval:      1    Jaccard:      0 is.tested: true
Secretome.gene.list                notA: 23166            inA: 1993        notA: 4383             inA: 229      odds.ratio: 0.60732     pval:      1    Jaccard: 0.034671       is.tested: true
SupportedIRFMergeClassified.list                   notA: 23705            inA: 1454        notA: 4151             inA: 461      odds.ratio: 1.8105      pval: 2.7293e-24        Jaccard: 0.075997       is.tested: true
tandem.gene.list                   notA: 20191            inA: 4968        notA: 3243             inA: 1369     odds.ratio: 1.7156      pval: 6.9682e-49        Jaccard: 0.1429 is.tested: true
```
# Actually plot the data
```
ggplot(data=test, aes(y=V2,x =V1))+geom_violin(size=2,fill = "#FF6666") +labs(title="Gene Expression of Genomic Strata",y="log(Fold Change)",x="") +coord_flip(ylim=c(-10,11)) + xlim( "Genes with Repeat-overlapping Exons", "All Predicted Effectors in DNA Transposons", "DOG-box Genes in DNA Transposons", "Predicted Effectors in DNA Transposons", "Known Effectors in DNA Transposons", "Genes in DNA Transposons", "All Predicted Effectors in LTR Retroelements", "DOG-box Genes in LTR Retroelements", "Predicted Effectors in LTR Retroelements", "Known Effectors in LTR Retroelements", "Genes in LTR Retroelements", "All Predicted Effectors in TDR", "DOG-box Genes in TDR", "Predicted Effectors in TDR", "Known Effectors in TDR", "Genes in TDR", "High Confidence HGT", "All Predicted Effectors", "DOG-box Genes",  "Predicted Effectors", "Known Effectors", "All Genes") + scale_colour_manual(values = cols) +  theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size = 46),title = element_text(size = 60),axis.ticks = element_line(colour = "black", size = 3),text = element_text( size = 30))


ggplot(data=SNPDensity4BoxPlot, aes(y=V2,x =V1))+geom_violin(size=2,fill = "#FF6666") +labs(title="SNP Density of Genomic Strata",y="log(SNP/CDS Length)",x="") +coord_flip(ylim=c(-10,0)) + xlim("Genes with Repeat-overlapping Exons", "All Predicted Effectors in DNA Transposons", "DOG-box Genes in DNA Transposons", "Predicted Effectors in DNA Transposons", "Known Effectors in DNA Transposons", "Genes in DNA Transposons","All Predicted Effectors in LTR Retroelements", "DOG-box Genes in LTR Retroelements", "Predicted Effectors in LTR Retroelements", "Known Effectors in LTR Retroelements", "Genes in LTR Retroelements", "All Predicted Effectors in TDR", "DOG-box Genes in TDR", "Predicted Effectors in TDR", "Known Effectors in TDR", "Genes in TDR", "High Confidence HGT", "All Predicted Effectors", "DOG-box Genes",  "Predicted Effectors", "Known Effectors", "All Genes") + scale_colour_manual(values = cols) +  theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size = 46),title = element_text(size = 60),axis.ticks = element_line(colour = "black", size = 3),text = element_text( size = 30))
```
