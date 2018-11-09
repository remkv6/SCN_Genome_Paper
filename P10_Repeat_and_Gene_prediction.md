# Repeat and Gene prediction
```
This started off as a masking of the genome with a C. elegans library, then again with a Nematoda library. Analysis became more comprehensive when RepeatModeler was included for de-novo identification. I did not get a consensii.classified file for my first two runs, so Repeatmasking was unclassified. Then obtained a classified file by running RepeatModeler to a finish in round 5. This was then used to mask the genome first with all repeats and then with all repeats/but simple. This gave me two different files for masking and running gene model prediction with Braker. I obtained three files of with differing gene models: Braker with an unmasked genome, Braker with a completely masked genome, and Braker with all repeats masked excluding simple repeats.
RepeatModeler

#This is the most up to date version for the genome738sl.polished.mitoFixed.fa.
#/data021/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler/ #file setup


ln -s ../18_mitochondria/swapMitoScaffolds/genome738sl.polished.mitoFixed.fa
module load repeatmodeler/1.0.8
/shared/software/GIF/programs/repeatmodeler/1.0.8/BuildDatabase -name test738mitoswapolished -engine ncbi -dir /data021/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler

#pbs script to run repeatmodeler


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N RepeatModeler
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
module load repeatmodeler
RepeatModeler -database  test738mitoswapolished -engine ncbi -pa 16
ssh condo "qstat -f ${PBS_JOBID} |head"

RepeatMasker

Repeatmasker.pbs script to mask all repeats and classify


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N RepeatMasker
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module load parallel
module load repeatmasker/4.0.6

RepeatMasker -pa 16  -gff -lib consensi.fa.classified ../genome738sl.polished.mitoFixed.fa

#in case you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"
Repeats table all repeats masked and classified


==================================================
file name: genome738sl.polished.mitoFixed.fa
sequences:           738
total length:  123846405 bp  (123846405 bp excl N/X-runs)
GC level:         37.74 %
bases masked:   42120134 bp ( 34.01 %)
==================================================
              number of      length   percentage
              elements*    occupied  of sequence
--------------------------------------------------
SINEs:              453        53028 bp    0.04 %
     ALUs            0            0 bp    0.00 %
     MIRs            0            0 bp    0.00 %

LINEs:             8757      2240944 bp    1.81 %
     LINE1        1216       160488 bp    0.13 %
     LINE2         184        27457 bp    0.02 %
     L3/CR1       3547      1468392 bp    1.19 %

LTR elements:      6413      3544979 bp    2.86 %
     ERVL            0            0 bp    0.00 %
     ERVL-MaLRs      0            0 bp    0.00 %
     ERV_classI      0            0 bp    0.00 %
     ERV_classII     0            0 bp    0.00 %

DNA elements:     50211      9216582 bp    7.44 %
    hAT-Charlie    282        30570 bp    0.02 %
    TcMar-Tigger     0            0 bp    0.00 %

Unclassified:    112290     23424794 bp   18.91 %

Total interspersed repeats: 38480327 bp   31.07 %


Small RNA:            0            0 bp    0.00 %

Satellites:        2362       433513 bp    0.35 %
Simple repeats:   47493      2621602 bp    2.12 %
Low complexity:   14422      1095821 bp    0.88 %
==================================================

* most repeats fragmented by insertions or deletions
 have been counted as one element


The query species was assumed to be homo
RepeatMasker version open-4.0.6 , default mode

run with rmblastn version 2.2.27+
The query was compared to classified sequences in "consensi.fa.classified"
RepBase Update 20160829, RM database version 20160829
RepeatMasker to mask all repeats except simple repeats.
#/data021/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler/RM_3002.WedNov301615082016/RepeatMasknosimple


ln -s ../consensi.fa.classified
ln -s ../../genome738sl.polished.mitoFixed.fa
PBS script


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N RepeatMasker
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module load parallel
module load repeatmasker/4.0.6

RepeatMasker -pa 16  -nolow -gff -lib consensi.fa.classified genome738sl.polished.mitoFixed.fa

#sh RepeatMasker.sh
# in case you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"
#Repeats Table masking no simple


==================================================
file name: genome738sl.polished.mitoFixed.fa
sequences:           738
total length:  123846405 bp  (123846405 bp excl N/X-runs)
GC level:         37.74 %
bases masked:   38724062 bp ( 31.27 %)
==================================================
              number of      length   percentage
              elements*    occupied  of sequence
--------------------------------------------------
SINEs:              453        52973 bp    0.04 %
     ALUs            0            0 bp    0.00 %
     MIRs            0            0 bp    0.00 %

LINEs:             8778      2246114 bp    1.81 %
     LINE1        1212       160488 bp    0.13 %
     LINE2         184        27457 bp    0.02 %
     L3/CR1       3550      1469470 bp    1.19 %

LTR elements:      6405      3544545 bp    2.86 %
     ERVL            0            0 bp    0.00 %
     ERVL-MaLRs      0            0 bp    0.00 %
     ERV_classI      0            0 bp    0.00 %
     ERV_classII     0            0 bp    0.00 %

DNA elements:     50255      9246459 bp    7.47 %
    hAT-Charlie    284        30694 bp    0.02 %
    TcMar-Tigger     0            0 bp    0.00 %

Unclassified:    112230     23438083 bp   18.93 %

Total interspersed repeats: 38528174 bp   31.11 %


Small RNA:            0            0 bp    0.00 %

Satellites:        2363       436084 bp    0.35 %
Simple repeats:    1290       264369 bp    0.21 %
Low complexity:       0            0 bp    0.00 %
==================================================

* most repeats fragmented by insertions or deletions
 have been counted as one element


The query species was assumed to be homo
RepeatMasker version open-4.0.6 , default mode

run with rmblastn version 2.2.27+
The query was compared to classified sequences in "consensi.fa.classified"
RepBase Update 20160829, RM database version 20160829

Braker on each of these repeatmasked genomes, and an unmasked genome

There are some major changes in gene models when repeats are masked, and so we tried predicting gene models with three different genomes. completely masked, all masked but simple, and completely unmasked.

Unmasked
#https://intranet.gif.biotech.iastate.edu/doku.php/people:remkv6:genome368-genome2692_synteny:braker_redone_differential_transcript_mapping


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N Braker738
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
sh ~/common_scripts/runBraker.sh /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R1.fq.gz /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R2.fq.gz /data021/GIF/remkv6/Baum/CamTechGenomeComparison/15_transcriptsTo738/Braker/genome.738.fa

# in case you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"

Result: 31245 genes
Completely masked\
#/data021/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler/braker/braker/genome738sl.polished.mitoFixed.noquiver.masked Braker Script


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N Braker738
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
sh ./runBraker.sh /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R1.fq.gz /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R2.fq.gz /data021/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler/braker/genome738sl.polished.mitoFixed.noquiver.masked.fa

ase you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"
Result: 20064 genes

Repeatmasked genome without simple repeats masked
#/data021/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler/RM_3002.WedNov301615082016/RepeatMasknosimple/braker


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N Braker738
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
sh ./runBraker.sh /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R1.fq.gz /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R2.fq.gz /data021/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler/RM_3002.WedNov301615082016/RepeatMasknosimple/braker/738.polish.mito.masked.4braker.fa

#ase you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"
Result: 20797 genes

Reset, get table classified for consensii repeats

I now have merged gene models for the ESTs, melissa's rnaseq, gland RNA-seq, and isoseq. I still need the actual repeat content of the genome though, since repeatmasker did not seem to associate the repeats types to the proper categories. Current output places most repeats in the unclassified interspersed repeats category.

For 738 genome, with slow prediction.

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler/masked.consensii.mod
ln -s genome738sl.polished.mitoFixed.fa
ln -s ../RM_3002.WedNov301615082016/consensi.fa.classified


==================================================
file name: genome738sl.polished.mitoFixed.fa
sequences:           738
total length:  123846405 bp  (123846405 bp excl N/X-runs)
GC level:         37.74 %
bases masked:   42614616 bp ( 34.41 %)
==================================================
              number of      length   percentage
              elements*    occupied  of sequence
--------------------------------------------------
SINEs:              476        54434 bp    0.04 %
     ALUs            0            0 bp    0.00 %
     MIRs            0            0 bp    0.00 %

LINEs:             9035      2266374 bp    1.83 %
     LINE1        1270       163370 bp    0.13 %
     LINE2         193        27800 bp    0.02 %
     L3/CR1       3599      1480092 bp    1.20 %

LTR elements:      6584      3621372 bp    2.92 %
     ERVL            0            0 bp    0.00 %
     ERVL-MaLRs      0            0 bp    0.00 %
     ERV_classI      0            0 bp    0.00 %
     ERV_classII     0            0 bp    0.00 %

DNA elements:     51378      9323978 bp    7.53 %
    hAT-Charlie    274        30115 bp    0.02 %
    TcMar-Tigger     0            0 bp    0.00 %

Unclassified:    114734     23700013 bp   19.14 %

Total interspersed repeats: 38966171 bp   31.46 %


Small RNA:            0            0 bp    0.00 %

Satellites:        2416       438254 bp    0.35 %
Simple repeats:   47437      2620401 bp    2.12 %
Low complexity:   14401      1099446 bp    0.89 %
==================================================

* most repeats fragmented by insertions or deletions
 have been counted as one element


The query species was assumed to be homo
RepeatMasker version open-4.0.6 , sensitive mode

run with rmblastn version 2.2.27+
The query was compared to classified sequences in "consensi.fa.classified"
RepBase Update 20160829, RM database version 20160829


For the JGI genome to compare


RepeatMasker -pa 16 -a RM_alignments -s -gc -noisy -gff -lib consensii.name.short.classified Hetgly_JGI.fasta

==================================================
file name: Hetgly_JGI.fasta
sequences:         10356
total length:  174460779 bp  (160029005 bp excl N/X-runs)
GC level:         35.80 %
bases masked:   56042736 bp ( 32.12 %)
==================================================
              number of      length   percentage
              elements*    occupied  of sequence
--------------------------------------------------
SINEs:              585        56120 bp    0.03 %
     ALUs            0            0 bp    0.00 %
     MIRs            0            0 bp    0.00 %

LINEs:            12146      2587370 bp    1.48 %
     LINE1        1673       187558 bp    0.11 %
     LINE2         186        23013 bp    0.01 %
     L3/CR1       5261      1765786 bp    1.01 %

LTR elements:      9867      4595676 bp    2.63 %
     ERVL            0            0 bp    0.00 %
     ERVL-MaLRs      0            0 bp    0.00 %
     ERV_classI      0            0 bp    0.00 %
     ERV_classII     0            0 bp    0.00 %

DNA elements:     68460     12173747 bp    6.98 %
    hAT-Charlie    345        26698 bp    0.02 %
    TcMar-Tigger     0            0 bp    0.00 %

Unclassified:    161318     32129467 bp   18.42 %

Total interspersed repeats: 51542380 bp   29.54 %


Small RNA:            0            0 bp    0.00 %

Satellites:        3742       560994 bp    0.32 %
Simple repeats:   57617      3146740 bp    1.80 %
Low complexity:   19702      1405900 bp    0.81 %
==================================================

* most repeats fragmented by insertions or deletions
 have been counted as one element


The query species was assumed to be homo
RepeatMasker version open-4.0.6 , sensitive mode

run with rmblastn version 2.2.27+
The query was compared to classified sequences in "consensii.name.short.classified"
RepBase Update 20160829, RM database version 20160829
```
### gene model comparisons
```


masked genome gene prediction
#/data013/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler/braker


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N Braker738
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
sh ./runBraker.sh /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R1.fq.gz /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R2.fq.gz /data021/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler/braker/genome738sl.polished.mitoFixed.noquiver.masked.fa

ssh condo "qstat -f ${PBS_JOBID} |head"
Only simple sequence repeats unmasked
#/data013/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler/RM_3002.WedNov301615082016/RepeatMasknosimple/braker


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N Braker738
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
sh ./runBraker.sh /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R1.fq.gz /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R2.fq.gz /data021/GIF/remkv6/Baum/CamTechGenomeComparison/22_RepeatModeler/RM_3002.WedNov301615082016/RepeatMasknosimple/braker/738.polish.mito.masked.4braker.fa

#ase you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"

Completely unmasked genome gene predictions
#/data013/GIF/remkv6/Baum/CamTechGenomeComparison/32_genePredictionComp/UnmaskedBraker


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N Braker738
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
sh ./runBraker.sh /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R1.fq.gz /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R2.fq.gz genome738sl.polished.mitoFixed.fa

# in case you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"
```
### Comparisons
```

#unmasked genome
[remkv6@condo genome738sl.polished.mitoFixed]$ awk '$3=="gene"' augustus.gff3 |wc
  29959  269631 1608855
#simple sequence repeats unmasked
[remkv6@condo 32_genePredictionComp]$ awk '$3=="gene"' maskedNosimp.gff3 |wc
  20797  187173 1112927
#all repeats masked
[remkv6@condo 32_genePredictionComp]$ awk '$3=="gene"' maskedall.gff3 |wc
  20064  180576 1075307
#unmasked genome with isoseq+rnaseq
[remkv6@condo isoseq_geneModels]$ awk '$3=="gene"' augustus.gff3 |wc
  31135  280215 1682150

We thought we would stick with the unmasked genome gene predictions. We chose this route because there seemed to be a lot of exon/transposon overlap, even though 20,000 is close to the C.elegans gene count.
```
  ![Gene Prediction Comparisons](assets/GenePredictionComparisons.png)


### Manual modification of gene gff
  ```
  #/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/1_genomeNgff
Manually changed gene 29880 and 29881 in the gff to represent one gene as 29880.  These genes made an alignment to another nematode species single protein, pat-12.

Subsequent stats


#how many genes were called in the unmasked genome compared to masked?
bedtools intersect -v -a  MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3  -b maskedall.gff3 |awk '$3=="gene"' |wc
   8215   73935  442420
  ```
