# Need to get a pca or DAPC analysis done for scn pop data

```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/41_PCA
#make the scaffold names numeric only
sed 's/scaffold_//g' ../14_GATK/Sortedrenamed.combined_variants.vcf >numbers4scaffs.vcf &

#get rid of samples with lots of missing data
module load vcftools
vcftools --vcf numbers4scaffs.vcf  --missing-indv --out missingness
#all samples were fine


#impute and phase snps so they are of high quality
module show GIF2/beagle/4.1
module load java
 java -jar -Xmx100g /shared/software/GIF/programs/beagle/4.1/beagle.jar gt=numbers4scaffs.vcf out=numbers4scaffsPhasedImputed nthreads=16 impute=true ibd=true

 gunzip numbers4scaffsPhasedImputed.vcf.gz


#make population names list
 less numbers4scaffsPhasedImputed.vcf |awk 'NR==11' |cut -f 10- |tr "\t" "\n" >pop_code_lines.txt
 G3
LY1
OP20
OP25
OP50
PA3
TN1
TN13
TN15
TN16
TN19
TN21
TN22
TN7
TN8



#make population virulence list
vi pop_code_vir.txt
0
1234567
123567
157
123567
0
1257
1367
123567
1257
1234567
1234567
1257
12567
1367
```

### Run SNPRelate in R
```
library(gdsfmt)
library(SNPRelate)
vcf.fn <- "numbers4scaffsPhasedImputed.vcf"
snpgdsVCF2GDS(vcf.fn, "numbers4scaffsPhasedImputed.gds", method="biallelic.only")
genofile <- snpgdsOpen("numbers4scaffsPhasedImputed.gds")
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, autosome.only=FALSE, num.thread=4)
pc.percent <- pca$varprop*100
tab <- data.frame(sample.id = pca$sample.id, EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
pca$sample.id
colfunc <- colorRampPalette(c("#FF0000", "#00FF00"))
pop_code <- scan("pop_code_vir.txt", what=character())
head(cbind(sample.id, pop_code))
tab <- data.frame(sample.id = pca$sample.id, pop = factor(pop_code)[match(pca$sample.id, sample.id)], EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
```

### Results

Data to be plotted
```
sample.id     pop          EV1         EV2
1         G3       0  0.215575508  0.26478131
2        LY1       0 -0.256114130 -0.26360100
3       OP20  123567 -0.503563794  0.13260288
4       OP25     157  0.005588751  0.72957492
5       OP50  123567 -0.478710305  0.15576553
6        PA3       0  0.238225464  0.26600616
7        TN1    1257  0.211046653 -0.09095782
8       TN13    1367  0.210646333 -0.12261938
9       TN15  123567  0.186324650 -0.18147198
10      TN16    1257  0.099395250 -0.11023895
11      TN19 1234567 -0.229543074 -0.19917007
12      TN21 1234567 -0.247652410 -0.23102874
13      TN22    1257  0.166140397 -0.08526872
14       TN7   12567  0.197232729 -0.05516117
15       TN8    1367  0.185407978 -0.20921298
```
Plot PCAVirulence <brk/>
![PCA Plot](assets/PCAVirulence.jpg)

<brk/>
Interpretation: <brk/>
```
It seems interesting that two of the super-virulent genotypes cluster together with an avirulent LY1.  Interestingly, the other two avirulent genotypes cluster together, but far away from LY1.  The upper left cluster has every the hg-type represented and cluster tightly.  This may indicate that all of these pathotypes have a common origin, plus they all have TN in their name (I have no idea what that means!).  OP50 and OP20 are clustered together, but the same pathotype is found with TN15 in the upper left cluster.  This may indicate two independent origins of this pathotype.  Then there is TN1 off by its lonesome.  So a possible 5 independent origins within these 15 populations.  
```
