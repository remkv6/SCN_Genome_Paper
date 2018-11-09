#  Is there any contamination in the 738 genome?

```
#/data021/GIF/remkv6/Baum/CamTechGenomeComparison/17_738contamRemoval
We are attempting to publish this genome and so we want to make sure that there are not any contaminating scaffolds present in the assembly.
The idea is to use blobtools to identify this contamination, which requires a bam file, a blast output file, and a genome fasta file.
This was done on the unpolished 738.genome, because I felt polishing would not make much difference and the bam files were already created.

```
Megablast
```
ln ../14_mergegenomes/genome.738.fa
fasta-splitter.pl --n-parts 24 --measure count genome.738.fa
cp ~/common_scripts/runMegablast.sh .
#delete the evalue parameter  from the above script
vi runMegablast.sh

#Blast pbs script
#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N 738blastNR
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module load parallel
module load ncbi-blast
parallel  <<EOF
./runMegablast.sh genome.738.part-01.fa
./runMegablast.sh genome.738.part-02.fa
./runMegablast.sh genome.738.part-03.fa
./runMegablast.sh genome.738.part-04.fa
./runMegablast.sh genome.738.part-05.fa
./runMegablast.sh genome.738.part-06.fa
./runMegablast.sh genome.738.part-07.fa
./runMegablast.sh genome.738.part-08.fa
./runMegablast.sh genome.738.part-09.fa
./runMegablast.sh genome.738.part-10.fa
./runMegablast.sh genome.738.part-11.fa
./runMegablast.sh genome.738.part-12.fa
./runMegablast.sh genome.738.part-13.fa
./runMegablast.sh genome.738.part-14.fa
./runMegablast.sh genome.738.part-15.fa
./runMegablast.sh genome.738.part-16.fa
./runMegablast.sh genome.738.part-17.fa
./runMegablast.sh genome.738.part-18.fa
./runMegablast.sh genome.738.part-19.fa
./runMegablast.sh genome.738.part-20.fa
./runMegablast.sh genome.738.part-21.fa
./runMegablast.sh genome.738.part-22.fa
./runMegablast.sh genome.738.part-23.fa
./runMegablast.sh genome.738.part-24.fa
EOF


# in case you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"


cat *.megablast.out > all.genome.nt.cul5.megablast.out
```

Blobtools
```
mkdir blob_tools
cd blob_tools
ln -s ../../15_transcriptsTo738/Braker/genome.738.fa
ln -s ../../15_transcriptsTo738/Braker/genome.738_sorted_rnaseq.bam
awk '$3>100' ../all.genome.nt.cul5.megablast.out >all.genome.nt.cul5.100score.megablast.out
#the pbs script to run blobtools, provided by arun

#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=4:00:00
#PBS -N BlobTools
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
module load python
module load parallel
module load bwa
module load samtools
module load blobtools
BAM=genome.738_sorted_rnaseq.bam
GENOME=genome.738.fa
BLAST=all.genome.nt.cul5.100score.megablast.out
NODES=/shared/software/GIF/programs/blobtools/nodes.dmp
NAMES=/shared/software/GIF/programs/blobtools/names.dmp

blobtools create \
  -i $GENOME \
  -b $BAM \
  -t $BLAST \
  --nodes $NODES \
  --names $NAMES \
  -o blobplot_out

mkdir -p blobplot_files

blobtools view \
  -i blobplot_out.blobDB.json \
  -o blobplot_files/

blobtools blobplot -i blobplot_out.blobDB.json -o blobplot_files/

grep -v '^#' blobplot_files/blobplot_out.blobDB.table.txt | cut -f 1,3 > blobDB.id.gc.txt
awk '$2 < 0.25' blobDB.id.gc.txt |   cut -f1 |   perl -lne 'print $_.",<20%"'   > blobDB.id.gc.catcolour.txt
awk '$2 >= 0.20 && $2 < 0.30' blobDB.id.gc.txt |   cut -f1 |   perl -lne 'print $_.",20-29%"'   >> blobDB.id.gc.catcolour.txt
awk '$2 >= 0.30 && $2 < 0.40' blobDB.id.gc.txt |   cut -f1 |   perl -lne 'print $_.",30-39%"'   >> blobDB.id.gc.catcolour.txt
awk '$2 >= 0.40 && $2 < 0.50' blobDB.id.gc.txt |   cut -f1 |   perl -lne 'print $_.",40-49%"'   >> blobDB.id.gc.catcolour.txt
awk '$2 >= 0.50 && $2 < 0.60' blobDB.id.gc.txt |   cut -f1 |   perl -lne 'print $_.",50-59%"'   >> blobDB.id.gc.catcolour.txt
awk '$2 >= 0.60 && $2 < 0.70' blobDB.id.gc.txt |   cut -f1 |   perl -lne 'print $_.",60-69%"'   >> blobDB.id.gc.catcolour.txt
awk '$2 >= 0.70 && $2 < 0.80' blobDB.id.gc.txt |   cut -f1 |   perl -lne 'print $_.",70-79%"'   >> blobDB.id.gc.catcolour.txt
awk '$2 >= 0.80 && $2 < 0.90' blobDB.id.gc.txt |   cut -f1 |   perl -lne 'print $_.",80-89%"'   >> blobDB.id.gc.catcolour.txt
awk '$2 >= 0.90 && $2 < 1.00' blobDB.id.gc.txt |   cut -f1 |   perl -lne 'print $_.",90-99%"'   >> blobDB.id.gc.catcolour.txt

blobtools covplot \
  -i blobplot_out.blobDB.json \
  -c all_reads.bam.cov \
  --catcolour blobDB.id.gc.catcolour.txt \
  --notitle \
  --ylabel WGA-resequencing-library \
  --xlabel WGS-resequencing-library \

issh condo "qstat -f ${PBS_JOBID} | head"

I tried blobplot with e-25 cutoff and with a 200score cutoff. the 100score cutoff gave the highest ratio of nematode to everything else. All else will have to be checked manually to see if the calls are real contamination or false.
```

  ![Blobtools output](assets/Blobtools.png)
