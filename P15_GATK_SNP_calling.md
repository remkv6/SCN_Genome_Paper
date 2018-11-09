# Have 15 populations of H. glycines, and can call snps to identify associations with virulence and population.
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/63_GATK
module purge

###GATK_00 modified###

#!/bin/bash
# Prepares the Reference Genome for mapping as well as for using it with GATK pipeline
# You need to supply the referece genome as REF below or as:
# ./GATK_00_PrepareRef.sh your_genome.fasta
module load GIF2/picard
module load samtools
module load bwa
module load bedtools2
module load parallel
module load python
REF="$1"
#index genome for (a) picard, (b) samtools and (c) bwa
#parallel <<FIL
java -Xmx100G -jar $PICARD_HOME/picard.jar CreateSequenceDictionary \
  REFERENCE=${REF} \
  OUTPUT=${REF%.*}.dict
samtools faidx ${REF}
bwa index -a bwtsw ${REF}
#FIL
#Create interval list (here 100 kb intervals)
fasta_length.py ${REF} > ${REF%.*}_length.txt
bedtools makewindows -w 100000 -g ${REF%.*}_length.txt > ${REF%.*}_100kb_coords.bed
java -Xmx100G -jar $PICARD_HOME/picard.jar BedToIntervalList \
  INPUT=${REF%.*}_100kb_coords.bed \
  SEQUENCE_DICTIONARY=${REF%.*}.dict \
  OUTPUT=${REF%.*}_100kb_gatk_intervals.list

#running the script in interactive node
sh GATK_00_PrepareRef.sh 1tandemcopyonly.masked.genome738.fasta
```
### Mapping and cleanup
```
#softlinked files below
G3_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/G3/FCC7GPJANXX-CHKPEI15080015_L1_1.fq.gz
G3_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/G3/FCC7GPJANXX-CHKPEI15080015_L1_2.fq.gz
LY1_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/LY1/LY1_S6_L006_R1_001.fastq.gz
LY1_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/LY1/LY1_S6_L006_R2_001.fastq.gz
OP20_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP20/OP20_S5_L006_R1_001.fastq.gz
OP20_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP20/OP20_S5_L006_R2_001.fastq.gz
OP25_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP25/FCC7GPJANXX-CHKPEI15080011_L1_1.fq.gz
OP25_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP25/FCC7GPJANXX-CHKPEI15080011_L1_2.fq.gz
OP50_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP50/FCC7GPJANXX-CHKPEI15080012_L1_1.fq.gz
OP50_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/OP50/FCC7GPJANXX-CHKPEI15080012_L1_2.fq.gz
PA3_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN1/FCC7GPJANXX-wHAIPI022662-14_L1_1.fq.gz
PA3_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN1/FCC7GPJANXX-wHAIPI022662-14_L1_2.fq.gz
TN1_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN1/FCC7GPJANXX-wHAIPI022662-14_L1_1.fq.gz
TN1_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN1/FCC7GPJANXX-wHAIPI022662-14_L1_2.fq.gz
TN7_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN7/FCC7GPJANXX-wHAIPI022663-35_L1_1.fq.gz
TN7_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN7/FCC7GPJANXX-wHAIPI022663-35_L1_2.fq.gz
TN8_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN8/TN8_S1_L006_R1_001.fastq.gz
TN8_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN8/TN8_S1_L006_R2_001.fastq.gz
TN13_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN13/FCC7GPJANXX-wHAIPI022664-37_L1_1.fq.gz
TN13_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN13/FCC7GPJANXX-wHAIPI022664-37_L1_2.fq.gz
TN15_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN15/TN15_S2_L006_R1_001.fastq.gz
TN15_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN15/TN15_S2_L006_R2_001.fastq.gz
TN16_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN16/FCC7GPJANXX-wHAIPI022665-46_L1_1.fq.gz
TN16_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN16/FCC7GPJANXX-wHAIPI022665-46_L1_2.fq.gz
TN19_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN19/FCC7GPJANXX-CHKPEI15080013_L1_1.fq.gz
TN19_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN19/FCC7GPJANXX-CHKPEI15080013_L1_2.fq.gz
TN21_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN21/TN21_S3_L006_R1_001.fastq.gz
TN21_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN21/TN21_S3_L006_R2_001.fastq.gz
TN22_1.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN22/TN22_S4_L006_R1_001.fastq.gz
TN22_2.fq.gz -> /work/GIF/archive1/Baum/101415_SNPSCN/RawReads/clean_reads/TN22/TN22_S4_L006_R2_001.fastq.gz

#changing names to work with Arun's script
mv FCC7GPJANXX-CHKPEI15080011_L1_1.fq.gz   OP25_1.fq.gz  
mv FCC7GPJANXX-CHKPEI15080011_L1_2.fq.gz     OP25_2.fq.gz
mv FCC7GPJANXX-CHKPEI15080012_L1_1.fq.gz    OP50_1.fq.gz
mv FCC7GPJANXX-CHKPEI15080012_L1_2.fq.gz    OP50_2.fq.gz
mv FCC7GPJANXX-CHKPEI15080013_L1_1.fq.gz    TN19_1.fq.gz
mv FCC7GPJANXX-CHKPEI15080013_L1_2.fq.gz    TN19_2.fq.gz
mv FCC7GPJANXX-CHKPEI15080014_L1_1.fq.gz    PA3_1.fq.gz
mv FCC7GPJANXX-CHKPEI15080014_L1_2.fq.gz    PA3_2.fq.gz
mv FCC7GPJANXX-CHKPEI15080015_L1_1.fq.gz G3_1.fq.gz
mv FCC7GPJANXX-CHKPEI15080015_L1_2.fq.gz    G3_2.fq.gz
mv FCC7GPJANXX-wHAIPI022662-14_L1_1.fq.gz   TN1_1.fq.gz
mv FCC7GPJANXX-wHAIPI022662-14_L1_2.fq.gz   TN1_2.fq.gz
mv FCC7GPJANXX-wHAIPI022663-35_L1_1.fq.gz   TN7_1.fq.gz
mv FCC7GPJANXX-wHAIPI022663-35_L1_2.fq.gz   TN7_2.fq.gz
mv FCC7GPJANXX-wHAIPI022664-37_L1_1.fq.gz   TN13_1.fq.gz
mv FCC7GPJANXX-wHAIPI022664-37_L1_2.fq.gz   TN13_2.fq.gz
mv FCC7GPJANXX-wHAIPI022665-46_L1_1.fq.gz   TN16_1.fq.gz
mv FCC7GPJANXX-wHAIPI022665-46_L1_2.fq.gz   TN16_2.fq.gz
mv LY1_S6_L006_R1_001.fastq.gz  LY1_1.fq.gz
mv LY1_S6_L006_R2_001.fastq.gz  LY1_2.fq.gz
mv OP20_S5_L006_R1_001.fastq.gz OP20_1.fq.gz
mv OP20_S5_L006_R2_001.fastq.gz OP20_2.fq.gz
mv TN8_S1_L006_R1_001.fastq.gz  TN8_1.fq.gz
mv TN8_S1_L006_R2_001.fastq.gz TN8_2.fq.gz
mv TN15_S2_L006_R1_001.fastq.gz TN15_1.fq.gz
mv TN15_S2_L006_R2_001.fastq.gz TN15_2.fq.gz
mv TN21_S3_L006_R1_001.fastq.gz TN21_1.fq.gz
mv TN21_S3_L006_R2_001.fastq.gz TN21_2.fq.gz
mv TN22_S4_L006_R1_001.fastq.gz TN22_1.fq.gz
mv TN22_S4_L006_R2_001.fastq.gz TN22_2.fq.gz
```
#### Copied and Aruns run file runBWA_and_Cleanup.sh
```
#!/bin/bash
module use /work/GIF/software/modules
module load bwa
module load samtools
module load GIF2/gatk
module load GIF2/picard

genome="/work/GIF/remkv6/Baum/CamTechGenomeComparison/63_GATK/1tandemcopyonly.masked.genome738.fasta"
R1="$1"
R2="$2"
OUT=$(echo $R1 |cut -f 1 -d "_");
SAM="${OUT}_Hg738_slp.sam"
bwa mem -M -t 16 $genome $R1 $R2 > ${SAM};
samtools view --threads 16 -b -o ${SAM%.*}.bam ${SAM}
#samtools sort -m 6G -o ${SAM%.*}_sorted.bam -T ${SAM%.*}_temp --threads 16 ${SAM%.*}.bam

FILE="${SAM%.*}.bam"
REF="${genome}"

SAMPLE=$(echo ${FILE} | cut -d "_" -f 1)
UNIT=$(echo ${FILE} | cut -d "_" -f 2)
RGLB="onetandem"
GATK=$GATK_HOME/GenomeAnalysisTK.jar

#mkdir -p /local/scratch/${USER}/${PBS_JOBID}
echo $TMPDIR

echo "Sorting BAM of ${FILE}"

if [ ! -f $PBS_O_WORKDIR/${FILE%.*}_picsort.bam ]; then
java -Djava.io.tmpdir=$TMPDIR -Xmx100G -jar $PICARD_HOME/picard.jar SortSam \
  TMP_DIR=${TMPDIR}\
  INPUT=${FILE} \
  OUTPUT=${FILE%.*}_picsort.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM='null' || {
  echo >&2 sorting failed for $FILE
  exit 1
}
fi

echo "Cleaning Alignment file of ${FILE}"
if [ ! -f $PBS_O_WORKDIR/${FILE%.*}_picsort_cleaned.bam ]; then
java -Djava.io.tmpdir=$TMPDIR -Xmx100G -jar $PICARD_HOME/picard.jar CleanSam \
  TMP_DIR=${TMPDIR} \
  INPUT=${FILE%.*}_picsort.bam \
  OUTPUT=${FILE%.*}_picsort_cleaned.bam \
  MAX_RECORDS_IN_RAM='null' || {
  echo >&2 cleaning failed for $FILE
  exit 1
}
fi

echo "Marking Duplicates of ${FILE}"
if [ ! -f $PBS_O_WORKDIR/${FILE%.*}_dedup.bam ]; then
java -Djava.io.tmpdir=$TMPDIR -Xmx100G -jar $PICARD_HOME/picard.jar MarkDuplicates \
  TMP_DIR=${TMPDIR} \
  INPUT=${FILE%.*}_picsort_cleaned.bam \
  OUTPUT=${FILE%.*}_dedup.bam \
  METRICS_FILE=${FILE%.*}_metrics.txt \
  ASSUME_SORTED=true \
  REMOVE_DUPLICATES=true \
  MAX_RECORDS_IN_RAM=5000000 || {
  echo >&2 deduplicating failed for $FILE
  exit 1
}
fi

echo "Adding RG info of ${FILE}"
if [ ! -f $PBS_O_WORKDIR/${FILE%.*}_dedup_RG.bam ]; then
java -Djava.io.tmpdir=$TMPDIR -Xmx100G -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups \
  TMP_DIR=${TMPDIR} \
  INPUT=${FILE%.*}_dedup.bam \
  OUTPUT=${FILE%.*}_dedup_RG.bam \
  RGID=${SAMPLE} RGLB=${RGLB} \
  RGPL=illumina \
  RGPU=${UNIT} \
  RGSM=${SAMPLE} \
  MAX_RECORDS_IN_RAM='null' \
  CREATE_INDEX=true || {
  echo >&2 RG adding failed for $FILE
  exit 1
}
fi

echo "Indel Realigner: create intervals of ${FILE}"
if [ ! -f $PBS_O_WORKDIR/${FILE%.*}_target_intervals.list ]; then
samtools index ${FILE%.*}_dedup_RG.bam

java -Djava.io.tmpdir=$TMPDIR -Xmx100G -jar $GATK \
  -T RealignerTargetCreator \
  -R ${REF} \
  -I ${FILE%.*}_dedup_RG.bam \
  -o ${FILE%.*}_target_intervals.list || {
echo >&2 Target intervels list generation failed for $FILE
exit 1
}
fi

echo "Indel Realigner: write realignments of ${FILE}"
if [ ! -f $PBS_O_WORKDIR/${FILE%.*}_realigned.bam ]; then
java -Djava.io.tmpdir=$TMPDIR -Xmx100G -jar $GATK \
  -T IndelRealigner \
  -R ${REF} \
  -I ${FILE%.*}_dedup_RG.bam \
  -targetIntervals ${FILE%.*}_target_intervals.list \
  -o ${FILE%.*}_realigned.bam || {
echo >&2 Indel realignment failed for $FILE
exit 1
}
fi

samtools index ${FILE%.*}_realigned.bam

echo "cleaning up of ${FILE%.*}"
#if your job stops midway move all the intermediate files into the main directory and comment out this section
#rewrite to check this folder instead of the main folder for restart

mkdir -p OriginalBAMfile IntermediateBAMfiles
mv ${SAM%.*}.bam OriginalBAMfile
mv ${FILE%.*}_picsort.bam IntermediateBAMfiles
mv ${FILE%.*}_picsort_cleaned.bam IntermediateBAMfiles
mv ${FILE%.*}_dedup.bam IntermediateBAMfiles
mv ${FILE%.*}_dedup_RG.bam IntermediateBAMfiles
mv ${FILE%.*}_target_intervals.list IntermediateBAMfiles
mv ${FILE%.*}_metrics.txt IntermediateBAMfiles

echo "All done!"
```
make list of commands
```
paste  <( ls -1 *_1.* )  <( ls -1 *_2* )  |sed 's/^/sh runBWA_and_Cleanup.sh /g' >runAndClean.list
```
### SNP density
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/63_GATK/GATK
Calculate the upper and lower
#Number of snps per 10kb region
less all_sorted.vcf |awk 'NF>10 {print $1,$2}' |grep -v "#"|awk '{print $1, $2/10000}' |sed 's/\./\t/g'|awk '{print $1,($2+1)*10000}' |uniq -c |sort -k1,1nr|awk '{print $1}' |summary.sh
Total:  2,258,506
Count:  12,165
Mean:   185
Median: 167
Min:    1
Max:    1,125

#This is for a standard deviation? calculation.
less all_sorted.vcf |awk 'NF>10 {print $1,$2}' |grep -v "#"|awk '{print $1, $2/10000}' |sed 's/\./\t/g'|awk '{print $1,($2+1)*10000}' |uniq -c |sort -k1,1nr|awk '{print ($1-185)*($1-185)}' |summary.sh
Total:  181,917,871
Count:  12,165
Mean:   14,954
Median: 5,476
Min:    0
Max:    883,600

#this was the number I used to distinguish sparse and dense regions
square root of 14,954 =~ 122.28

#snp dense regions
122+185==307
less all_sorted.vcf |awk 'NF>10 {print $1,$2}' |grep -v "#"|awk '{print $1, $2/10000}' |sed 's/\./\t/g'|awk '{print $1,($2+1)*10000}' |uniq -c |sort -k1,1nr|awk '$1>307' |wc
   1673    5019   45266

#snp sparse regions
185-122=63   
less all_sorted.vcf |awk 'NF>10 {print $1,$2}' |grep -v "#"|awk '{print $1, $2/10000}' |sed 's/\./\t/g'|awk '{print $1,($2+1)*10000}' |uniq -c |sort -k1,1nr|awk '$1<63' |wc
   1649    4947   44682


#genes in snp dense regions
bedtools intersect -wo -a <(sed 's/SCAFFOLD/scaffold/g' SnpDense.gff) -b ../../58_Renamatorium/1_genomeNgff/fixed.augustus.gff3 |awk '$11=="gene"{print $17}' |sort|uniq|wc
   4484    4484  182082
#genes in snp sparse regions   
bedtools intersect -wo -a <(sed 's/SCAFFOLD/scaffold/g' SnpSparse.gff) -b ../../58_Renamatorium/1_genomeNgff/fixed.augustus.gff3 |awk '$11=="gene"{print $17}' |sort|uniq|wc
   4094    4094  166430

less all_sorted.vcf |awk 'NF>10 {print $1,$2}' |grep -v "#"|awk '{print $1, $2/10000}' |sed 's/\./\t/g'|awk '{print $1,($2+1)*10000}' |uniq -c |sort -k1,1nr|awk '$1<63'  |awk '{print
$2,"SNPsparse","10kbregion",$3-9999,$3,"+",".","."}' |tr " " "\t" >SnpSparse.gff

less all_sorted.vcf |awk 'NF>10 {print $1,$2}' |grep -v "#"|awk '{print $1, $2/10000}' |sed 's/\./\t/g'|awk '{print $1,($2+1)*10000}' |uniq -c |sort -k1,1nr|awk '$1>307'  |awk '{print
 $2,"SNPsparse","10kbregion",$3-9999,$3,"+",".","."}' |tr " " "\t" >SnpDense.gff


#Extracting genes for ontologizer analysis
bedtools intersect -wo -a <(sed 's/SCAFFOLD/scaffold/g' SnpSparse.gff) -b ../../58_Renamatorium/1_genomeNgff/fixed.augustus.gff3 |awk '$11=="gene"{print $17}' |sed 's/ID=//'g|sed 's/;
/\t/g'|cut -f 1 >SnpSparseGene.list
[remkv6@condo020 GATK]$ bedtools intersect -wo -a <(sed 's/SCAFFOLD/scaffold/g' SnpDense.gff) -b ../../58_Renamatorium/1_genomeNgff/fixed.augustus.gff3 |awk '$11=="gene"{print $17}' |sed 's/ID=//'g|sed 's/;/
\t/g'|cut -f 1 >SnpDenseGene.list
```
### Ontologizer analysis
```
#

ln -s ../../../57_secretome/ontologenizer/Ontologizer.jar
ln -s ../../../57_secretome/ontologenizer/population
ln -s ../../../57_secretome/ontologenizer/simpleformat.ids
ln -s ../../../57_secretome/ontologenizer/Secretome2/go.obo

cp ../../../58_Renamatorium/1_genomeNgff/geneRenamer.sh .
less geneRenamer.sh |awk -F "/" '{print $1,$3,$2,$4}' |sed 's| |/|5' |sed 's| |/|4' |sed 's| |/|3' |sed 's/augustus.aa/SnpDenseGene.list/g' >geneRenamerSnpDense.sh
less geneRenamer.sh |awk -F "/" '{print $1,$3,$2,$4}' |sed 's| |/|5' |sed 's| |/|4' |sed 's| |/|3' |sed 's/augustus.aa/SnpSparseGene.list/g' >geneRenamerSnpSparse.sh

 cp ../SnpDenseGene.list .
 cp ../SnpSparseGene.list .
 sh geneRenamerSnpDense.sh &
 sh geneRenamerSnpSparse.sh &

java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s SnpDenseGene.list
java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s SnpSparseGene.list
SNP DENSE ENRICHMENT

ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0006259      13669   406     4484    182     1924    311     2       false   1.3808539309159711E-58  1.7743973012270228E-55  0.0     "DNA metabolic process"
GO:0071897      13669   268     4484    138     1169    240     3       false   9.6749246796453E-41     1.243227821334421E-37   1.6933808542033813E-272 "DNA biosynthetic process"
GO:0006260      13669   259     4484    135     1052    223     3       false   6.564384927731032E-40   8.435234632134376E-37   3.566462500540798E-254  "DNA replication"
GO:0016779      13669   368     4484    153     854     195     1       false   2.041267961464588E-30   2.6230293304819957E-27  1.0763066966234276E-252 "nucleotidyltransferase activity"
GO:0033554      13669   242     4484    104     839     154     2       false   5.192975266468513E-29   6.672973217412039E-26   4.1781108156920284E-218 "cellular response to stress"
GO:0006950      13669   275     4484    108     852     154     1       false   6.98578283788228E-27    8.97673094667873E-24    6.526442148199629E-232  "response to stress"
GO:0003676      13669   1969    4484    470     3206    600     2       false   1.4851873514230628E-22  1.9084657465786358E-19  0.0     "nucleic acid binding"
GO:0090305      13669   208     4484    97      1090    234     1       false   3.70638937296031E-20    4.7627103442539984E-17  5.9966330597453525E-230 "nucleic acid phosphodiester bond hydrolysis"
GO:0004518      13669   206     4484    97      352     108     2       false   5.908898406883229E-17   7.592934452844949E-14   4.3091954321275545E-103 "nuclease activity"
GO:0034061      13669   243     4484    134     388     155     2       false   2.4904790976879567E-16  3.2002656405290244E-13  1.0290624705873474E-110 "DNA polymerase activity"
GO:0046483      13669   1258    4484    255     2376    367     1       false   2.0007351080598627E-12  2.5709446138569235E-9   0.0     "heterocycle metabolic process"
GO:0006725      13669   1252    4484    254     2376    367     1       false   2.2090176591499994E-12  2.8385876920077494E-9   0.0     "cellular aromatic compound metabolic process"
GO:0090304      13669   1090    4484    234     2546    402     2       false   9.825044929236727E-12   1.2625182734069195E-8   0.0     "nucleic acid metabolic process"
GO:0008152      13669   3414    4484    513     4820    628     1       false   1.8181255344778808E-11  2.336291311804077E-8    0.0     "metabolic process"
GO:0016788      13669   350     4484    108     1372    254     1       false   2.9294924997966974E-11  3.764397862238756E-8    0.0     "hydrolase activity, acting on ester bonds"
GO:1901363      13669   3206    4484    600     4351    723     1       false   7.65267319423347E-11    9.833685054590008E-8    0.0     "heterocyclic compound binding"
GO:0097159      13669   3206    4484    600     4351    723     1       false   7.65267319423347E-11    9.833685054590008E-8    0.0     "organic cyclic compound binding"
GO:0005488      13669   4351    4484    723     6535    964     1       false   5.77125202285055E-10    7.416058849362957E-7    0.0     "binding"
GO:1901360      13669   1273    4484    256     2948    460     1       false   3.3367297752615196E-9   4.287697761211053E-6    0.0     "organic cyclic compound metabolic process"
GO:0006139      13669   1233    4484    251     2903    459     5       false   6.552333328093669E-9    8.419748326600364E-6    0.0     "nucleobase-containing compound metabolic process"
GO:0034641      13669   1435    4484    276     2874    444     2       false   1.2254558405937353E-8   1.57471075516295E-5     0.0     "cellular nitrogen compound metabolic process"
GO:0016740      13669   1117    4484    227     3085    485     1       false   1.1584782072900878E-7   1.488644496367763E-4    0.0     "transferase activity"
GO:0006508      13669   599     4484    101     1391    163     1       false   1.8947346053685188E-7   2.4347339678985468E-4   0.0     "proteolysis"
GO:0065008      13669   143     4484    30      1177    101     1       false   5.299859697080495E-7    6.810319710748436E-4    2.356573527676257E-188  "regulation of biological quality"
GO:0050896      13669   852     4484    154     4820    628     1       false   2.0838200314718424E-6   0.0026777087404413176   0.0     "response to stimulus"
GO:0003674      13669   6535    4484    964     7396    1043    1       false   2.3810347324718354E-6   0.0030596296312263085   0.0     "molecular_function"
GO:0051716      13669   807     4484    150     3152    433     2       false   4.0068090244665254E-6   0.005148749596439485    0.0     "cellular response to stimulus"
GO:0004386      13669   113     4484    28      410     50      1       false   5.339511856391947E-6    0.006861272735463652    3.319843551520327E-104  "helicase activity"
GO:0006298      13669   167     4484    89      233     104     1       false   1.5059530744956897E-5   0.019351497007269614    8.448366374772887E-60   "mismatch repair"
GO:0044249      13669   1199    4484    221     2416    371     2       false   1.932901404749459E-5    0.024837783051030547    0.0     "cellular biosynthetic process"
GO:0032200      13669   64      4484    25      138     32      1       false   3.6980535502864804E-5   0.04751998812118127     6.067366163410429E-41   "telomere organization"
SNP SPARSE ENRICHMENT


ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0044848      13490   5       4094    5       4820    602     1       false   2.995134538554313E-5    0.0392362624550615      4.622185398361942E-17   "biological phase"
GO:0098763      13490   5       4094    5       22      5       2       false   3.79737221842484E-5     0.049745576061365405    3.79737221842484E-5     "mitotic cell cycle phase"
```
### Effector overlap
```
bedtools intersect -wo -a SnpSparse.gff -b ../../58_Renamatorium/4_effectors/effector.gmapped.filtered.gff3 |less

 231416  233797  .       +       .       ID=lcl|esthggfha15A10|Pioneer(missing5').path1;Name=lcl|esthggfha15A10|Pioneer(missing5')       2381
 233069  234060  .       +       .       ID=GLAND7|Pioneer,15A10family.path4;Name=GLAND7|Pioneer,15A10family     991
 208772  209392  .       -       .       ID=GLAND5|Pioneer,11A06/2D01/24A12/16B09/30E03/22C12/4D06/29D09family.path4;Name=GLAND5|Pioneer,11A06/2D01/24A12/16B09/30E03/22C12/4D06/29D09family     620
 315937  324384  .       +       .       ID=GLAND5|Pioneer,11A06/2D01/24A12/16B09/30E03/22C12/4D06/29D09family.path5;Name=GLAND5|Pioneer,11A06/2D01/24A12/16B09/30E03/22C12/4D06/29D09family     4063
 93200   96555   .       +       .       ID=GLAND7|Pioneer,15A10family.path2;Name=GLAND7|Pioneer,15A10family     3355
 260250  260621  .       -       .       ID=GLAND9|Pioneer.path1;Name=GLAND9|Pioneer     371
 47761   48545   .       -       .       ID=lcl|flhggfha30D08|Pioneer,16A01/21E12family.path4;Name=lcl|flhggfha30D08|Pioneer,16A01/21E12family   784
 155991  163878  .       -       .       ID=lcl|esthggfha15A10|Pioneer(missing5').path3;Name=lcl|esthggfha15A10|Pioneer(missing5')       3877
 72474   73361   .       -       .       ID=lcl|esthggfha8C06|Ran-bindingprotein.path1;Name=lcl|esthggfha8C06|Ran-bindingprotein 887
 111122  112009  .       -       .       ID=lcl|esthggfha8C06|Ran-bindingprotein.path2;Name=lcl|esthggfha8C06|Ran-bindingprotein 887
 29572   30093   .       +       .       ID=lcl|flhggfha21E12|Pioneer,30D08/16A01family.path5;Name=lcl|flhggfha21E12|Pioneer,30D08/16A01family   428
 29565   30349   .       +       .       ID=lcl|flhggfha30D08|Pioneer,16A01/21E12family.path3;Name=lcl|flhggfha30D08|Pioneer,16A01/21E12family   435
 29523   30394   .       +       .       ID=lcl|flhggfha16A01|Pioneer,30D08/21E12family.path5;Name=lcl|flhggfha16A01|Pioneer,30D08/21E12family   477
 9938    10737   .       -       .       ID=lcl|flhggfha4G12|CLE-likepeptide.path2;Name=lcl|flhggfha4G12|CLE-likepeptide 62
 9892    10702   .       -       .       ID=lcl|flhggfha2B10|SYV46,CLE-like.path2;Name=lcl|flhggfha2B10|SYV46,CLE-like   108


bedtools intersect -wo -a SnpDense.gff -b ../../58_Renamatorium/4_effectors/effector.gmapped.filtered.gff3 |less
79130   80815   .       +       .       ID=GLAND8|Pioneer.path1;Name=GLAND8|Pioneer     870
150933  152580  .       +       .       ID=lcl|flhggfha4D09|Pioneer.path1;Name=lcl|flhggfha4D09|Pioneer 1647
124     5064    .       -       .       ID=lcl|flhggfha5D06|Pioneer.path1;Name=lcl|flhggfha5D06|Pioneer 4940
2267    4002    .       +       .       ID=lcl|flhggfha12H04|Pioneer.path1;Name=lcl|flhggfha12H04|Pioneer       1735
141546  143565  .       -       .       ID=lcl|flhggfha17G01|Pioneer.path1;Name=lcl|flhggfha17G01|Pioneer       2019
658032  658784  .       +       .       ID=lcl|flhggfha24A12|Pioneer,2D01/11A06/16B09/22C12/30E03/4D06/29D09family.path1;Name=lcl|flhggfha24A12|Pioneer,2D01/11A06/16B09/22C12/30E03/4D06/29D09family   752
657951  658792  .       +       .       ID=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family.path4;Name=lcl|flhggfha11A06|Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family   841
147828  151070  .       -       .       ID=lcl|flhggfha13C08|Cellulase(ENG-5).path1;Name=lcl|flhggfha13C08|Cellulase(ENG-5)     2172
45870   50816   .       +       .       ID=GLAND13|Pioneer,invertase.path2;Name=GLAND13|Pioneer,invertase       815
142780  147884  .       -       .       ID=lcl|flGSB3|28B03variant.path1;Name=lcl|flGSB3|28B03variant   5104
142780  147840  .       -       .       ID=lcl|flhggfha28B03|Skp-1like,8H07family.path1;Name=lcl|flhggfha28B03|Skp-1like,8H07family     5060
119787  120740  .       -       .       ID=GLAND10|Cellulosebindingprotein.path1;Name=GLAND10|Cellulosebindingprotein   739
33441   34284   .       -       .       ID=GLAND4|Pioneer,32E03/10A07/6E07/13A06/27D09/20G04family,also11A06inG.rostoc.path1;Name=GLAND4|Pioneer,32E03/10A07/6E07/13A06/27D09/20G04family,also11A06inG.rostoc   843
60051   67845   .       +       .       ID=lcl|flhggfha4G05|Pioneer,30G12/25A01family.path1;Name=lcl|flhggfha4G05|Pioneer,30G12/25A01family     7794
60100   67727   .       +       .       ID=lcl|flhggfha25A01|Pioneer,30G12/4G05family.path1;Name=lcl|flhggfha25A01|Pioneer,30G12/4G05family     7627
60088   67839   .       +       .       ID=lcl|flhggfha30G12|Pioneer,4G06/25A01family.path1;Name=lcl|flhggfha30G12|Pioneer,4G06/25A01family     7751
247132  247921  .       +       .       ID=lcl|flhggfha7E05|Pioneer.path1;Name=lcl|flhggfha7E05|Pioneer 789
52983   53858   .       -       .       ID=lcl|flhggfha5D08|Pioneer.path1;Name=lcl|flhggfha5D08|Pioneer 875
240597  241910  .       +       .       ID=lcl|flhggfha31A08|Pioneer(insitusuggestsnervering).path1;Name=lcl|flhggfha31A08|Pioneer(insitusuggestsnervering)     1313
863649  865974  .       -       .       ID=lcl|flhggfha8A07|Pioneer.path1;Name=lcl|flhggfha8A07|Pioneer 2325
113776  114728  .       +       .       ID=lcl|flhggfha33A09|Pioneer.path1;Name=lcl|flhggfha33A09|Pioneer       952
490980  492249  .       -       .       ID=lcl|flhggfha27D09|Pioneer,10A07/20G04/10A07/6E07/13A06/32E03family.path1;Name=lcl|flhggfha27D09|Pioneer,10A07/20G04/10A07/6E07/13A06/32E03family     1269
490993  492293  .       -       .       ID=lcl|flhggfha10A07|Pioneer,20G04/13A06/6E07/32E03family.path2;Name=lcl|flhggfha10A07|Pioneer,20G04/13A06/6E07/32E03family     1300
490993  492259  .       -       .       ID=lcl|flhggfha20G04|Pioneer,10A07/27D09/13A06/6E07/32E03family.path1;Name=lcl|flhggfha20G04|Pioneer,10A07/27D09/13A06/6E07/32E03family 1266
```
