# Have the updated effector fasta list from Tom now, so need to gmap and figure up all overlap.
### Gmap of effector genes
```
#Had to modify the effector.fa file so that GLAND## would appear rather than the indistinct Esophageal name.
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/29_effectorMapping

#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 04:00:00
#SBATCH -J gmap.genome.738_0
#SBATCH -o gmap.genome.738_0.o%j
#SBATCH -e gmap.genome.738_0.e%j
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module use /work/GIF/software/modules
module load gsnap
gmap_build -d genome738sl.polished.mitoFixed  -D /work/GIF/remkv6/Baum/CamTechGenomeComparison/29_effectorMapping/ genome738sl.polished.mitoFixed.fa
gmap -D /work/GIF/remkv6/Baum/CamTechGenomeComparison/29_effectorMapping/ -d genome738sl.polished.mitoFixed -B 5 -t 16 --input-buffer-size=1000000 --output-buffer-size=1000000 -f 2  effector.fa >effectors.gmapped.gff3
qstat -f ${PBS_JOBID} |head
```
### Bedtools comparisons
```
Some of the mapping of the effectors seems erroneous with some genes spanning 30-50kb. I removed those in the following two calculations
How many
less effectors.gmapped.gff3 |awk '$3=="gene"' |awk '$5-$4>200 && $5-$4<10000' |wc
    130    1170   16113
#there were 26 super large alignments that I decided to eliminate because they are less likely, and blur any associations.  
less effectors.gmapped.gff3 |awk '$5-$4>10000'|wc
     26     234    4650





#number of effectors that map to genes in merged gene models
bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../52_functional/augustusFunctionalAnnotation.gff3 |awk '$12=="mRNA" && $3=="gene" {print $9}' |sed 's/\.path/\t/g' |cut -f 1 |sort|uniq|wc
     75      75    3195





#number of genes that the above 75 effectors map to.  
bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../52_functional/augustusFunctionalAnnotation.gff3 |awk '$12=="mRNA" && $3=="gene"' |cut -f 18 |sed 's/\./\t/g'|cut -f 1 |sort|uniq |wc
    111     111    1057




#How many genes overlapped with each effector?  
bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../50_tandemMummer/augustus_renamed_with_putative_function.gff3 |awk '$12=="mRNA" && $3=="gene"' |cut -f 9,18- |cut -f 1 |uniq -c |sort -k1,1nr >uniqcNumGenesPerEffector.txt
#
      9 ID=GLAND5|Pioneer,11A06/2D01/24A12/16B09/30E03/22C12/4D06/29D09family
      6 ID=lcl|esthggfha15A10|Pioneer(missing5')
      5 ID=lcl|flhggfha21E12|Pioneer,30D08/16A01family
      4 ID=GLAND13|Pioneer,invertase
      4 ID=GLAND9|Pioneer
      3 ID=GLAND7|Pioneer,15A10family
      3 ID=lcl|esthggfha8C06|Ran-bindingprotein
      3 ID=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family
      3 ID=lcl|flhggfha30D08|Pioneer,16A01/21E12family
      3 ID=lcl|flhggfha5D08|Pioneer
      2 ID=GLAND10|Cellulosebindingprotein
      2 ID=GLAND1|GNATprotein
      2 ID=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family
      2 ID=lcl|flhggfha17G06|Pioneer
      2 ID=lcl|flhggfha27D09|Pioneer,10A07/20G04/10A07/6E07/13A06/32E03family
      2 ID=lcl|flhggfha2A05|VAP
      2 ID=lcl|flhggfha30C02|Pioneer
      2 ID=lcl|flhggfha33E05|Pioneer,34B08/23G12/6E07/10A07family
      2 ID=lcl|flhggfha3B05|Cellulosebindingprotein
      2 ID=lcl|flhggfha3H07|Ubiquitinextension
      2 ID=lcl|flhggfha4E02|Pioneer,10C02family
      2 ID=lcl|flhggfha4F01|Annexin
      2 ID=lcl|flhggfha4G05|Pioneer,30G12/25A01family
      2 ID=lcl|flhggfha4G12|CLE-likepeptide

#

#How many repeats overlap with the new effector genes
bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '$3=="gene"' |cut -f 10,13,14,18 |sed 's/Target "Motif://g' |sed 's/".*//g' |sort|uniq -c |sort -k1,1nr |cut -f 4 |sort|uniq -c |sort -k1,1nr|wc
    164     328    3687


#what repeat types overlap with the new effector genes
bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '$3=="gene"' |cut -f 10,13,14,18 |sed 's/Target "Motif://g' |sed 's/".*//g' |sort|uniq -c |sort -k1,1nr |cut -f 4 |sort|uniq -c |sort -k1,1nr|less
#
     22 A-rich
     22 rnd-4_family-976
     12 rnd-5_family-1918 -- all overlap with 6 effector alignments here these are intronic insertions
     10 rnd-3_family-228 -- this overlaps 1 effector with 2 positions in the genome
      9 rnd-2_family-24 -- these are intronic.
      9 rnd-4_family-1265 -- LINEs, but not all are overlapping with effector exons
      8 rnd-4_family-1152 -- 8C06 has three alignments, all of which have exons that are this repeat (DNA/TcMar-Tc2, a grab bag of retro and dna transposons).  the other two effector overlaps are intronic
      7 (TAAT)n
      6 (AATT)n
      6 rnd-4_family-1155
      5 rnd-4_family-138
      5 rnd-4_family-71
      5 rnd-5_family-535
#

#976 is interesting, which effectors is it found in?
bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '$3=="gene"' |cut -f 10,13,14,18 |sed 's/Target "Motif://g' |sed 's/".*//g' |sort|uniq -c |sort -k1,1nr |grep "rnd-4_family-976" |cut -f 1,2,3 |sort|wc
22

#how do these repeats overlap with effector gene functions
 bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '$3=="gene"' |cut -f 9,18- |grep "rnd-4_family-976" |awk '{print $1}' |sed 's/,.*//g' |uniq| grep -f - <(bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../52_functional/augustusFunctionalAnnotation.gff3 |awk '$12=="mRNA" && $3=="gene"' |cut -f 9,18-) |less
#
ID=lcl|flhggfha2D01|Pioneer,11A06/24A12/16B09/30E03/22C12/4D06/29D09family.path2;Name=lcl|flhggfha2D01|Pioneer,11A06/24A12/16B09/30E03/22C12/4D06/29D09family   ID=g8082.t1;Parent=g8082        739
ID=lcl|flhggfha2D01|Pioneer,11A06/24A12/16B09/30E03/22C12/4D06/29D09family.path2;Name=lcl|flhggfha2D01|Pioneer,11A06/24A12/16B09/30E03/22C12/4D06/29D09family   ID=g8082.t1;Parent=g8082; Note=NA; Note=Similar to Q86DE9_HETGL: Putative gland protein G24A12 n%3D1 (Heterodera glycines); ID=g8082.t1;Parent=g8082;   739
ID=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family.path1;Name=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family   ID=g20414.t1;Parent=g20414      737
ID=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family.path1;Name=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family   ID=g20414.t1;Parent=g20414; Note=NA; Note=NA; ID=g20414.t1;Parent=g20414;       737
ID=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family.path2;Name=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family   ID=g20413.t1;Parent=g20413      739
ID=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family.path2;Name=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family   ID=g20413.t1;Parent=g20413; Note=NA; Note=NA; ID=g20413.t1;Parent=g20413;       739
ID=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family.path3;Name=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family   ID=g14326.t1;Parent=g14326      744
ID=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family.path3;Name=lcl|flhggfha16B09|Pioneer,30E03/22C12/11A06/2D01/24A12/4D06/29D09family   ID=g14326.t1;Parent=g14326; Note=NA; Note=NA; ID=g14326.t1;Parent=g14326;       744
ID=lcl|flhggfha24A12|Pioneer,2D01/11A06/16B09/22C12/30E03/4D06/29D09family.path1;Name=lcl|flhggfha24A12|Pioneer,2D01/11A06/16B09/22C12/30E03/4D06/29D09family   ID=g8110.t1;Parent=g8110        752
ID=lcl|flhggfha24A12|Pioneer,2D01/11A06/16B09/22C12/30E03/4D06/29D09family.path1;Name=lcl|flhggfha24A12|Pioneer,2D01/11A06/16B09/22C12/30E03/4D06/29D09family   ID=g8110.t1;Parent=g8110; Note=NA; Note=Similar to Q86DG8_HETGL: Putative gland protein G11A06 n%3D2 (Heterodera glycines); ID=g8110.t1;Parent=g8110;   752
ID=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family.path1;Name=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family   ID=g14325.t1;Parent=g14325      744
ID=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family.path1;Name=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family   ID=g14325.t1;Parent=g14325; Note=NA; Note=NA; ID=g14325.t1;Parent=g14325;       744
ID=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family.path2;Name=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family   ID=g14326.t1;Parent=g14326      744
ID=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family.path2;Name=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family   ID=g14326.t1;Parent=g14326; Note=NA; Note=NA; ID=g14326.t1;Parent=g14326;       744
ID=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family.path3;Name=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family   ID=g20414.t1;Parent=g20414      737
ID=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family.path3;Name=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family   ID=g20414.t1;Parent=g20414; Note=NA; Note=NA; ID=g20414.t1;Parent=g20414;       737
ID=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family.path4;Name=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family   ID=g20413.t1;Parent=g20413      739
ID=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family.path4;Name=lcl|flhggfha22C12|Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family   ID=g20413.t1;Parent=g20413; Note=NA; Note=NA; ID=g20413.t1;Parent=g20413;       739
ID=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family.path1;Name=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family   ID=g26790.t1;Parent=g26790      746
ID=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family.path1;Name=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family   ID=g26790.t1;Parent=g26790; Note=NA; Note=Similar to Q86DG8_HETGL: Putative gland protein G11A06 n%3D2 (Heterodera glycines); ID=g26790.t1;Parent=g26790;       746
ID=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family.path2;Name=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family   ID=g26791.t1;Parent=g26791      746
ID=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family.path2;Name=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family   ID=g26791.t1;Parent=g26791; Note=NA; Note=Similar to Q86DG8_HETGL: Putative gland protein G11A06 n%3D2 (Heterodera glycines); ID=g26791.t1;Parent=g26791;       746
ID=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family.path4;Name=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family   ID=g8110.t1;Parent=g8110        752
ID=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family.path4;Name=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family   ID=g8110.t1;Parent=g8110; Note=NA; Note=Similar to Q86DG8_HETGL: Putative gland protein G11A06 n%3D2 (Heterodera glycines); ID=g8110.t1;Parent=g8110;   752
ID=GLAND6|Pioneer,4D06/29D09/11A06/2D01/16B09/30E03/22C12/24A12family.path2;Name=GLAND6|Pioneer,4D06/29D09/11A06/2D01/16B09/30E03/22C12/24A12family     ID=g26009.t1;Parent=g26009      686
ID=GLAND6|Pioneer,4D06/29D09/11A06/2D01/16B09/30E03/22C12/24A12family.path2;Name=GLAND6|Pioneer,4D06/29D09/11A06/2D01/16B09/30E03/22C12/24A12family     ID=g26009.t1;Parent=g26009; Note=NA; Note=NA; ID=g26009.t1;Parent=g26009;       686
incomplete
#

#new effectors against the new gene models
bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |awk '$12=="gene" && $3=="gene"' |cut -f 18 |sort |uniq |wc
    111     111    1168


#get the borders of the effectors to find repeat enrichment.
bedtools flank -b 500 -i <(awk '$3=="gene"' effector.gmapped.filtered.gff3) -g ../50_tandemMummer/chr.length.genome >effector.500bp.borders.gff3

#Which repeatmodeler repeats are found in the 500bp borders?
bedtools intersect -wo -a effector.500bp.borders.gff3 -b ../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '$3=="gene"' |cut -f 10,13,14,18 |sed 's/Target "Motif://g' |sed 's/".*//g' |sort|uniq -c |sort -k1,1nr |cut -f 4 |sort|uniq -c |sort -k1,1nr|less
#
     22 rnd-4_family-976
     13 A-rich
     10 rnd-4_family-352
      8 rnd-4_family-1265
      6 (T)n
      5 rnd-3_family-315
      5 rnd-4_family-115
      5 rnd-4_family-149
      5 rnd-5_family-623
      4 rnd-3_family-146
      4 rnd-4_family-268
      4 rnd-4_family-424
      4 rnd-4_family-650
      4 rnd-4_family-791
      4 rnd-5_family-4142
#


#What is going on with family 976?
bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '$3=="gene"' |grep "Motif:rnd-4_family-976" |cut -f 1,4,5,9 |uniq|less

#
000008  542724  543617  ID=lcl|fl_hggfha2D01|Pioneer,.path2;Name=lcl|fl_hggfha2D01|Pioneer, -- completely contained within transposon 976
```
![000008_541000_2d01](assets/000008_541000_2d01.png)

```
000543  19773   20628   ID=lcl|fl_hggfha16B09|Pioneer,.path1;Name=lcl|fl_hggfha16B09|Pioneer,
000543  18375   19218   ID=lcl|fl_hggfha16B09|Pioneer,.path2;Name=lcl|fl_hggfha16B09|Pioneer, -- both of these are both almost wholly contained within transposons.  I have a picture of this already.  shows two transposons overlapping two effector genes.  both effector genes are expressed.
53syntmer       223208  224057  ID=lcl|fl_hggfha16B09|Pioneer,.path3;Name=lcl|fl_hggfha16B09|Pioneer, -- transposons comprising large percentage of effector.   
```
![53syntmer-220000_22c12](assets/53syntmer-220000_22c12.png)
```

53syntmer       221877  222621  ID=lcl|fl_hggfha22C12|Pioneer,.path1;Name=lcl|fl_hggfha22C12|Pioneer,
53syntmer       223238  223982  ID=lcl|fl_hggfha22C12|Pioneer,.path2;Name=lcl|fl_hggfha22C12|Pioneer,
000543  19803   20540   ID=lcl|fl_hggfha22C12|Pioneer,.path3;Name=lcl|fl_hggfha22C12|Pioneer,
000543  18450   19189   ID=lcl|fl_hggfha22C12|Pioneer,.path4;Name=lcl|fl_hggfha22C12|Pioneer,
000899  881602  882465  ID=lcl|fl_hggfha11A06|Pioneer,.path1;Name=lcl|fl_hggfha11A06|Pioneer,  -- this and the one below are completely enclosed within mule family 976 transposons
```
![000899-881000_11a06.png](assets/000899-881000_11a06.png)
```
000899  882686  883550  ID=lcl|fl_hggfha11A06|Pioneer,.path2;Name=lcl|fl_hggfha11A06|Pioneer,
000008  657951  658792  ID=lcl|fl_hggfha11A06|Pioneer,.path4;Name=lcl|fl_hggfha11A06|Pioneer,
000300  93021   93852   ID=Esophageal.path2;Name=Esophageal -- smaller fragment of transposon is in the exon of an effector gene that is expressed.
000297K 42068   42809   ID=Esophageal.path1;Name=Esophageal --If this follows the same pattern as the other tandemly duplicated mules/effectors, then gene 13620 may be an effector also.The two large repeatmasker tracks above are 95% identical with 82% coverage, yet one of them does not have an effector mapping to it.Gene13620 has 99% coverage with 77% identity to esophageal gland 5.  SO this is a double mule/effector combo also.
```
![000297k_38918..45011](assets/000297k_38918..45011.png)
```

42syntmer       70434   71193   ID=Esophageal.path2;Name=Esophageal
```
![42syntmer_69440..73440_esophageal](assets/42syntmer_69440..73440_esophageal.png)
```
000543  19773   20625   ID=lcl|fl_hggfha30E03|Pioneer,.path1;Name=lcl|fl_hggfha30E03|Pioneer, -- this one is not completely overlapped with a mule, but a considerable portion of the effector is.  
```
![000543_19001..22001.png](assets/000543_19001..22001.png)
```
000450K 85733   86683   ID=lcl|fl_hggfha4D06|Pioneer,.path1;Name=lcl|fl_hggfha4D06|Pioneer, -- small mule fragment in exon
```
![000450k_84000..88000.png](assets/000450k_84000..88000.png)

```
#what is going on with family 1265 -- LINE transposon
000257K 60051   67845   ID=lcl|fl_hggfha4G05|Pioneer,.path1;Name=lcl|fl_hggfha4G05|Pioneer,  -- repeat in intron
000037  319860  320376  ID=lcl|fl_hggfha21E12|Pioneer,.path1;Name=lcl|fl_hggfha21E12|Pioneer, -- ~70bp effector contained within transposon,expressed as part of gene2878
000252  162371  162889  ID=lcl|fl_hggfha21E12|Pioneer,.path2;Name=lcl|fl_hggfha21E12|Pioneer, -- 700bp effector contained within transposon, expressed as part of gene12368
{{:people:remkv6:000252-162371_21e12.png?direct&600|}}
000059  144640  145158  ID=lcl|fl_hggfha21E12|Pioneer,.path3;Name=lcl|fl_hggfha21E12|Pioneer, -- 700bp effector contained within transposon, expressed as part of gene4199

000059  151352  151870  ID=lcl|fl_hggfha21E12|Pioneer,.path4;Name=lcl|fl_hggfha21E12|Pioneer,
000491K 29572   30093   ID=lcl|fl_hggfha21E12|Pioneer,.path5;Name=lcl|fl_hggfha21E12|Pioneer,
000257K 60100   67727   ID=lcl|fl_hggfha25A01|Pioneer,.path1;Name=lcl|fl_hggfha25A01|Pioneer, -- repeat in intron same as above
000252  162115  162896  ID=lcl|fl_hggfha30D08|Pioneer,.path1;Name=lcl|fl_hggfha30D08|Pioneer,
000491K 38281   39065   ID=lcl|fl_hggfha30D08|Pioneer,.path2;Name=lcl|fl_hggfha30D08|Pioneer,These three are an artifact of a rediculously long effector gene mapping.When the giant intron is removed, a ~100bp exon overlaps with a repeat.   Doubt these are real.
000491K 29565   30349   ID=lcl|fl_hggfha30D08|Pioneer,.path3;Name=lcl|fl_hggfha30D08|Pioneer,
000491K 47761   48545   ID=lcl|fl_hggfha30D08|Pioneer,.path4;Name=lcl|fl_hggfha30D08|Pioneer,
35syntmer       166630  167402  ID=lcl|fl_hggfha30D08|Pioneer,.path5;Name=lcl|fl_hggfha30D08|Pioneer, -- 150bp exons within repeat.  not sure if it is 1265 wihtout condo up.
000257K 60088   67839   ID=lcl|fl_hggfha30G12|Pioneer,.path1;Name=lcl|fl_hggfha30G12|Pioneer,
000059  144594  145459  ID=lcl|fl_hggfha16A01|Pioneer,.path1;Name=lcl|fl_hggfha16A01|Pioneer,
000059  151051  151916  ID=lcl|fl_hggfha16A01|Pioneer,.path2;Name=lcl|fl_hggfha16A01|Pioneer,
35syntmer       161865  167444  ID=lcl|fl_hggfha16A01|Pioneer,.path3;Name=lcl|fl_hggfha16A01|Pioneer,
000252  162070  162941  ID=lcl|fl_hggfha16A01|Pioneer,.path4;Name=lcl|fl_hggfha16A01|Pioneer,
000491K 29523   30394   ID=lcl|fl_hggfha16A01|Pioneer,.path5;Name=lcl|fl_hggfha16A01|Pioneer,
```

![family976overlapwitheffectorsexample.png](assets/family976overlapwitheffectorsexample.png)

```
the repeatmasker family here is round4 family 976, and they are the two overlapping transposons with the effector genes/track. Interproscan classified these genes as MULE transposases. How cool is that?
?
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52

#where, how big, and how many are there of family 976 in the genome.
grep "rnd-4_family-976" ../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '{if($4>$5) print $1,$2,$3,$5,$4,$6,$7,$8,$9; else print $0}'|awk '{print $5-$4, $0}' |sort -k1,1nr |awk '{print $1}' |summary.sh
Total:  100,950
Count:  324
Mean:   311
Median: 103
Min:    5
Max:    2,335

#this is the location of known effectors that have repeatfamily 976.  They appear to be ~700bp.  
bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '$3=="gene"' |grep "Motif:rnd-4_family-976" |cut -f 1,4,5 |uniq|awk '{print $3-$2,$0}' |sort|less
#
620 000045      208772  209392
737 000543      19803   20540
739 000543      18450   19189
741 000297K     12398   13139
741 000297K     42068   42809
744 53syntmer   221877  222621
744 53syntmer   223238  223982
752 000008      658032  658784
759 42syntmer   70434   71193
831 000300      93021   93852
841 000008      657951  658792
843 000543      18375   19218
8447 000045     315937  324384
849 53syntmer   223208  224057
849 53syntmer   223208  224057
852 000543      19773   20625
855 000543      19773   20628
863 000899      881602  882465
864 000899      882686  883550
893 000008      542724  543617
950 000450K     85733   86683
#
that is all of them.
#the above represents the alignments for 9 effectors
bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '$3=="gene"' |grep "Motif:rnd-4_family-976" |cut -f 9 |sed 's/,/\t/g' |cut -f 1 |sort|uniq|wc
      9       9     237

#supp file of 976 loci in genome.
bedtools intersect -wo -a effector.gmapped.filtered.gff3 -b ../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '$3=="gene"' |grep "Motif:rnd-4_family-976" |cut -f 1,4,5 |uniq|awk '{print $3-$2,$0}' |sort >length976InKnownEffectors.txt


#length stats of overlap to known effectors that have repeatfamily 976.  This is not repeat size, just the size of the overlap
#note duplicate from above was removed.
less length976InKnownEffectors.txt |awk '{print $1}'|summary.sh
Total:  24,514
Count:  21
Mean:   1,167
Median: 841
Min:    620
Max:    8,447
```

### running meme on effector genes
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/29_effectorMapping
?
1
2
3
4
5

meme effector.fa -dna -nmotifs 10 -minsites 5  -maxsites 1000000 -minw 7 -maxw 300 -revcomp -maxsize 100000000 -oc effector.meme

#Tom asked to have this done with the proteins
TransDecoder.LongOrfs -t ../effector.fa
meme longest_orfs.pep  -protein -nmotifs 20 -minsites 5  -maxsites 1000000 -minw 7 -maxw 300 -maxsize 100000000 -oc effector.meme
```

### Restart Effector Meme analysis to network
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/29_effectorMapping/effectorMemeProt/TomTranslations/effector.known.meme

#reorders the meme.txt format to show two columns, motif and effectorName.
grep -A 11 "BLOCKS" meme.txt |cut -d " " -f 1,2 | grep -v BL |grep -v "-" |grep -v // |sed 's/\t//g' |sed 's/,//g'|sed 's/(.*//g'|sed '/^$/d' |awk '{if ($1 =="Motif") {a=$0;next} else {print a"\t"$0}}'  >knownEffectorNetwork.txt

mkdir FIMO
#copied and pasted excel formatted fimo file.
awk '{print "Motif "$1"\t"$3}' FimoMotifs.list |awk  'NR>1' |sed 's/\./\t/g' | cut -f 1,2 |sort|uniq >FimoMotifs.network
cat FimoMotifs.network ../knownEffectorNetwork.txt  |tr "\t" "," >Intermixed.network

awk '{print $9}' Secretome.gff3|sed 's/;/\t/g'|sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort|uniq >Secretome.gene.list
```
### reanalyzing gene model overlap for effector alignments
```
Exon to exon overlap is the only suitable measure to identify genes that are not just in effector introns.
?
1
2
3
4
5
6

Hmm, could only get 121 genes that overlapped effector alignments at exons, so likely some of the gene calls were missed or were absent/low coverage in rna-seq.
bedtools intersect -wo -a effectors.gmapped.gff3 -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$3=="exon" && $15=="exon" {print $9,$21}' |awk '{print $2}' |sed 's/;/\t/g' |cut -f 1 |sed 's/T/G/g' |sed 's/\./\t/2' |cut -f 1 |sort|uniq|wc
   121     121    2541
bedtools intersect -wo -a effectors.gmapped.gff3 -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$3=="exon" && $15=="exon" {print $21,$9}' |sed 's/ID=//g'|sed 's/lcl|//g'|sed 's/\./\t/2'| sed 's/T/G/1' |awk '{print $1,$3}' |sed 's/|/\t/1' |sed 's/\./\t/2' |awk '{print $1,$2,$3}' |sort |uniq |wc -l
166
bedtools intersect -wo -a effectors.gmapped.gff3 -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$3=="exon" && $15=="exon" {print $21,$9}' |sed 's/ID=//g'|sed 's/lcl|//g'|sed 's/\./\t/2'| sed 's/T/G/1' |awk '{print $1,$3}' |sed 's/|/\t/1' |sed 's/\./\t/2' |awk '{print $1,$2,$3}' |sort |uniq >EffectorsVsGenes.list
```

### meme fimo effector retest

```
 Need to test if removing some redundant transcript isoforms will have a significant impact on meme motif call, and thus the fimo/effectorome analysis.
 #/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/21_MemeFimoEffectorome/01_cdHit/

ln -s ../../../29_effectorMapping/effectorMemeProt/TomTranslations/effector.pep
cd-hit -i effector.pep -d 100 -o effector.pep.cdhit
===================================
total seq: 80
longest and shortest : 638 and 49
Total letters: 18955
Sequences have been sorted

80 finished     69 clusters
writing new database
writing clustering information
program completed !

Total CPU time 0
===================================
 less effector.pep.cdhit.clstr |awk '$1==0' |cut -d " " -f 2- |sed 's/,/\t/g' |sed 's/>//g'|awk '{print $1}' |sed 's/\.\.\.//g'|sort|uniq|grep -f - effector.pep |sed 's/>//g'| cdbyank effector.pep.cidx >Reduced
EffectorPep.fa

meme ReducedEffectorPep.fa -protein -nmotifs 20 -minsites 5 -maxsites 1000000 -minw 7 -maxw 300 -maxsize 100000000 -oc effector.meme


It makes the motifs slightly different, but the same sequences are still covered in motifs. I do not think it will make much of a difference to redo the analysis. If I had to guess the result of this, it would provide smaller motifs for fimo, thus identifying more genes whilst raising background noise. 3 of the 4 signal motifs are still there. I ran FIMO with this html file on the predicted genes, and got a fairly poor overlap of 88 genes when comparing to the old analysis 88/314 old 88/274 new.

Deciding to do meme with a background file of the predicted genes to see if the meme motifs change, and if fimo results get better.
```

### Meme Fimo with background
```
#filtered on Q value and number of motifs/gene must be at least 2.
  less fimo.txt |awk '$9<2e-3' |awk '{print $2,$3}' |sort|uniq|awk '{print $2}' |grep -v "\.t2" |grep -v "\.t3" |sort|uniq -c |awk '$1>1 {print $2}' |grep -w -f - <(less fimo.txt |awk '$9<2e-3' |awk '{print $2,$3}') |sort|uniq >FilteredNetwork
#combines with motifs identified in known effectors
tail -n 113 IntermixedNetwork2 |cat FilteredNetwork - >FilteredNetworkIntermixed

#The above network was lacking single connected nodes, and so many of the meme motifs were lacking connected genes.  

#making a network that includes single motif genes, needed higher q value cutoff
 less fimo.txt |awk '$9<1e-4' |awk '{print $2,$3}' |sort|uniq|awk '{print $2}' |grep -v "\.t2" |grep -v "\.t3" |sort|uniq -c |awk ' {print $2}' |grep -w -f - <(less fimo.txt |awk '$9<2e-3' |awk '{print $2,$3}' ) |sort|uniq >SinglesAllowedFilteredNetwork
 tail -n 113 IntermixedNetwork2 |cat SinglesAllowedFilteredNetwork - >SinglesAllowedFilteredNetworkIntermixed
awk '{print $2}' SinglesAllowedFilteredNetwork >SinglesAllowedFilteredNetwork.list
```
