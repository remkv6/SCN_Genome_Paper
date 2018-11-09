# Mapping Melissa's transcripts to the two SCN genomes modified with synteny

```
#/data021/GIF/remkv6/CamTechGenomeComparison/12_differentialTranscriptMapping/2692
ln -s ../10_scaffold_extend/finalsynt.2692.namemod.nodup.fa-test
#/data021/GIF/remkv6/CamTechGenomeComparison/12_differentialTranscriptMapping/368
ln -s ../10_scaffold_extend/finalsynt.368.genome.fa
#/data021/GIF/remkv6/CamTechGenomeComparison/12_differentialTranscriptMapping/
ln -s ../../SplicedLeaders/TrinityGguidedAll/nematode_transcripts_full_Trinity_genome_guided.fasta

creating gmap databases and running
#############################368genome##########################
#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -N GMAP-SCN-368
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
module load gmap-gsnap/20160404
gmap_build -d 368.genome  -D /data021/GIF/remkv6/CamTechGenomeComparison/12_differentialTranscriptMapping/368/ finalsynt.368.genome.fa
gmap -D /data021/GIF/remkv6/CamTechGenomeComparison/12_differentialTranscriptMapping/368/368.genome -d 368.genome -B 5 -t 16 --input-buffer-size=1000000 --output-buffer-size=1000000 -f psl --split-output=368 ../nematode_transcripts_full_Trinity_genome_guided.fasta

ssh condo "qstat -f ${PBS_JOBID} |head"



#####################2692###################
#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -N GMAP-SCN-269
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
#module use /shared/software/GIF/modules
module purge
module load gmap-gsnap/20160404
gmap_build -d 2692.genome  -D /data021/GIF/remkv6/CamTechGenomeComparison/12_differentialTranscriptMapping/2692/ finalsynt.2692.namemod.nodup.fa-test
gmap -D /data021/GIF/remkv6/CamTechGenomeComparison/12_differentialTranscriptMapping/2692/2692.genome -d 2692.genome -B 5 -t 16 --input-buffer-size=1000000 --output-buffer-size=1000000 -f psl --split-output=2692 ../nematode_transcripts_full_Trinity_genome_guided.fasta

ssh condo "qstat -f ${PBS_JOBID} |head"
```

### Post processing of mapping data using isoforms

```
grep ">" nematode_transcripts_full_Trinity_genome_guided.fasta |awk '{print $1}'>trans.names
cat 368.mult 368.transloc 368.uniq |awk '{print $10}'|sort|uniq >transcripts-mapped.368
cat trans.names 368/transcripts-mapped.368 |sed 's/>//g'|sort|uniq -c |sort -k 1,1|awk '$1=="1" {print $2}'>transcripts.notfound.in368
cat 2692.mult 2692.transloc 2692.uniq |awk '{print $10}'|sort|uniq >transcripts-mapped.2692
while read line; do echo "awk '"\$10"=="\"$line\""' transcripts-mapped.2692 >>transcriptsmissingin368";done < ../transcripts.notfound.in368 >test.sh
test.sh
awk '{print $10,$14}' transcripts-mapped.2692 |sort|uniq|tr " " "\t" |awk '{print $2}'|sort|uniq -c |awk '{print $2,$1}'|sort -k 1,1 -V -r>gene_counts_scaffold.2692only
awk '{print $10,$14}' transcriptsmissingin368 |sort|uniq|tr " " "\t" |awk '{print $2}'|sort|uniq -c |sort -k 1,1 -V -r|awk '{print $2, $1}'|sort -k 1,1 -V -r>gene_counts_scaffold.368only
join -a 1 -1 1 -2 1 <(sort gene_counts_scaffold.2692only) <(sort gene_counts_scaffold.368only) |awk 'NF==3'|sort -V -r -k 3,2>genesin269scaffadd
ln -s ../../05_iadhore/Nematode/output/forgffsubtraction/genome.2692.masked.chr.len
join -a 1 -2 1 -2 1 <(sort genesin269scaffadd) <(sort genome.2692.masked.chr.len) >scaf.2692.368.length
grep -v -f ../../11finalsynth2692selfiadhore/Nematode/output/STILLsyntenic  scaf.2692.368.length >unsyntenicscaffs
```

### genes, no isoforms considered

```
mv transcripts-mapped.368 isoforms-mapped.368
sed 's/_i.*//g' isoforms-mapped.368 >transcripts-mapped.368
#just changes the names for mapping
grep ">" nematode_transcripts_full_Trinity_genome_guided.fasta |awk '{print $1}' |sed -e 's/_i.*//g' -e 's/>//g' |sort|uniq>trans.names
cat trans.names transcripts-mapped.368 |sort|uniq -c|awk '$1=="1" {print $2}' >transcriptsNotin368
cat 2692.mult 2692.transloc 2692.uniq |awk '{print $10}'|sed 's/_i.*//g'>transcripts-mapped.2692
cat transcriptsNotin368 transcripts-mapped.2692 |sort|uniq -c|awk '$1=="2" {print $2}' >transcriptsNot368In2692
wc -l transcriptsNot368In2692
6599

#making table --> scaffold, total genes on 2692 scaffold, total genes on 2692scaffold missing from 368, 2692 scaffold length
awk '{print $14}' 2692.uniq.genesonly |sort |uniq -c |awk '{print $2,$1}'>2692.scaffold.names.totalct
 awk '{print $14}' transcriptsNot368In2692.alldata |sort|uniq -c |awk '{print $2,$1}'>2692.scaffold.names.missing368ct
join -1 1 -2 1 <(sort 2692.scaffold.names.totalct) <(sort 2692.scaffold.names.missing368ct)>scaf.2692.368
join -1 1 -2 1 <(sort scaf.2692.368) <(sort genome.2692.masked.chr.len) >scaf.2692.368.chrlen

#total  MB added if scaffolds
awk 'BEGIN {FS=" "}; {i+=$4} END {print i}' scaf.2692.368.chrlen
24635762
#total number of scaffolds that would be added
wc -l scaf.2692.368.chrlen
371
```
 I am going to run bwa mem first to remove duplicate scaffolds from the synteny extended 2692 genome, prior to adding these 371 scaffolds to the 368-synt-extended genome.   This should be mostly effective for repetitive scaffolds.

 ### Removing identical scaffolds

 ```
 ln -s ../10_scaffold_extend/finalsynt.2692.namemod.nodup.fa

 ##PBS script
 #!/bin/bash
 #PBS -l nodes=1:ppn=16
 #PBS -l walltime=1:00:00
 #PBS -N bwa-mem-2692
 #PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
 cd $PBS_O_WORKDIR
 ulimit -s unlimited
 module use /shared/software/GIF/modules
 module load bwa

 bwa index /data021/GIF/remkv6/Baum/CamTechGenomeComparison/13_bwaDupRemoval/finalsynt.2692.namemod.nodup.fa
 bwa mem -t 16 -a /data021/GIF/remkv6/Baum/CamTechGenomeComparison/13_bwaDupRemoval/finalsynt.2692.namemod.nodup.fa /data021/GIF/remkv6/Baum/CamTechGenomeComparison/13_bwaDupRemoval/finalsynt.2692.namemod.nodup.fa >out.bwa

 # in case you need stats after job completion retain this as last line
 ssh condo "qstat -f ${PBS_JOBID} |head"




 #just get the scaffold to scaffold information, but remove self primary alignments
 grep K out.bwa|awk '$10="";{print $0}'|awk '$2==0 && $1!=$3' |sort -V -k 1,1>nonself-primary-alignments
 #attach chr length to first column scaffold and swap cols 1 with 3
 join -a 1 -1 1 -2 1 nonself-primary-alignments finalsynt.2692.len |sort -k 2,2 -V|awk '{print $3,$2,$1,substr($0,index($0,$4))}' |sort -k 1,1 -V>join.1
 #append chr length to first column scaffold again, and swap cols 1 with 3 -- remove duplicates caused by substr error.
 join -a 1 -1 1 -2 1 join.1 finalsynt.2692.len |awk '{print $3,$2,$1,substr($0,index($0,$4))}' |awk '{print $1, $2, $3, $8,$9
 ,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}' >join.2

 #print the scaffolds to col 1 that are identical and the smaller of the two
 cat <(awk '$15>$16 {print $3,$0}' join.2) <(awk '$16>$15 {print $1,$0}' join.2) >scaff.to.remove.all
 less scaff.to.remove.all
 002110K 002110K 0 000055K 1 0 20435M * 0 0 * NM:i:0 MD:Z:20435 AS:i:20435 XS:i:20435 20435 359167
 002073K 002073K 0 000122K 1 0 21321M * 0 0 * NM:i:0 MD:Z:21321 AS:i:21321 XS:i:21321 21321 235111
 002074K 002074K 0 000122K 1 0 21321M * 0 0 * NM:i:0 MD:Z:21321 AS:i:21321 XS:i:21321 21321 235111
 001822K 001822K 0 000123K 1 0 27983M * 0 0 * NM:i:0 MD:Z:27983 AS:i:27983 XS:i:27983 27983 294888
 001642K 001642K 0 000151K 1 0 31121M * 0 0 * NM:i:0 MD:Z:31121 AS:i:31121 XS:i:31121 31121 204560
 001526K 001526K 0 000227K 1 0 33695M * 0 0 * NM:i:0 MD:Z:33695 AS:i:33695 XS:i:33695 33695 182592
 001862K 001862K 0 000231K 1 0 27022M * 0 0 * NM:i:0 MD:Z:27022 AS:i:27022 XS:i:27022 27022 169888
 001415K 001415K 0 000324K 1 0 36261M * 0 0 * NM:i:0 MD:Z:36261 AS:i:36261 XS:i:36261 36261 140201
 001782K 001782K 0 000414K 1 0 28598M * 0 0 * NM:i:0 MD:Z:28598 AS:i:28598 XS:i:28598 28598 130357
 001465K 001465K 0 000423K 1 0 34777M * 0 0 * NM:i:0 MD:Z:34777 AS:i:34777 XS:i:34777 34777 110320
 001447K 001447K 0 000461K 1 0 35598M * 0 0 * NM:i:0 MD:Z:35598 AS:i:35598 XS:i:35598 35598 102563
 001170K 001170K 0 000540K 1 0 43391M * 0 0 * NM:i:0 MD:Z:43391 AS:i:43391 XS:i:43391 43391 72508
 001430K 001430K 0 000543K 1 0 35834M * 0 0 * NM:i:0 MD:Z:35834 AS:i:35834 XS:i:35834 35834 88778
 001808K 001808K 0 000548K 1 0 28158M * 0 0 * NM:i:0 MD:Z:28158 AS:i:28158 XS:i:28158 28158 93428
 001809K 001809K 0 000548K 1 0 28158M * 0 0 * NM:i:0 MD:Z:28158 AS:i:28158 XS:i:28158 28158 93428
 002129K 002129K 0 000681K 1 0 20027M * 0 0 * NM:i:0 MD:Z:20027 AS:i:20027 XS:i:20027 20027 72921
 001935K 001935K 0 000688K 1 0 25410M * 0 0 * NM:i:0 MD:Z:25410 AS:i:25410 XS:i:25410 25410 73438
 002441K 002441K 0 000692K 1 0 9974M * 0 0 * NM:i:0 MD:Z:9974 AS:i:9974 XS:i:9974 9974 72619
 001599K 001599K 0 000693K 1 0 32262M * 0 0 * NM:i:0 MD:Z:32262 AS:i:32262 XS:i:32262 32262 72943
 002017K 002017K 0 000706K 1 0 23141M * 0 0 * NM:i:0 MD:Z:23141 AS:i:23141 XS:i:23141 23141 71584
 001445K 001445K 0 000735K 1 0 35602M * 0 0 * NM:i:0 MD:Z:35602 AS:i:35602 XS:i:35602 35602 68542
 001819K 001819K 0 000746K 1 0 27994M * 0 0 * NM:i:0 MD:Z:27994 AS:i:27994 XS:i:27994 27994 67883
 001338K 001338K 0 000780K 1 0 37797M * 0 0 * NM:i:0 MD:Z:37797 AS:i:37797 XS:i:37797 37797 64540
 001339K 001339K 0 000780K 1 0 37797M * 0 0 * NM:i:0 MD:Z:37797 AS:i:37797 XS:i:37797 37797 64540
 001237K 001237K 0 000798K 1 0 41141M * 0 0 * NM:i:0 MD:Z:41141 AS:i:41141 XS:i:41141 41141 63469
 001078K 001078K 0 000807K 1 0 47300M * 0 0 * NM:i:0 MD:Z:47300 AS:i:47300 XS:i:47300 47300 62705
 001502K 001502K 0 000827K 1 0 34406M * 0 0 * NM:i:0 MD:Z:34406 AS:i:34406 XS:i:34406 34406 40853
 002194K 002194K 0 000835K 1 0 18015M * 0 0 * NM:i:0 MD:Z:18015 AS:i:18015 XS:i:18015 18015 60747
 001450K 001450K 0 000915K 1 0 35490M * 0 0 * NM:i:0 MD:Z:35490 AS:i:35490 XS:i:35490 35490 55647
 001764K 001764K 0 000934K 1 0 28741M * 0 0 * NM:i:0 MD:Z:28741 AS:i:28741 XS:i:28741 28741 54396
 001008K 001008K 0 000991K 1 0 51149M * 0 0 * NM:i:0 MD:Z:51149 AS:i:51149 XS:i:51149 51149 51883
 002163K 002163K 0 001025K 1 0 19284M * 0 0 * NM:i:0 MD:Z:19284 AS:i:19284 XS:i:19284 19284 50275
 001900K 001900K 0 001029K 1 0 26594M * 0 0 * NM:i:0 MD:Z:26594 AS:i:26594 XS:i:26594 26594 50128
 001334K 001334K 0 001032K 1 0 37916M * 0 0 * NM:i:0 MD:Z:37916 AS:i:37916 XS:i:37916 37916 49831
 002323K 002323K 0 001254K 1 0 13597M * 0 0 * NM:i:0 MD:Z:13597 AS:i:13597 XS:i:13597 13597 40635
 001651K 001651K 0 001357K 1 0 30981M * 0 0 * NM:i:0 MD:Z:30981 AS:i:30981 XS:i:30981 30981 37256
 002047K 002047K 0 001409K 1 0 22456M * 0 0 * NM:i:0 MD:Z:22456 AS:i:22456 XS:i:22456 22456 36436
 002369K 002369K 0 001523K 1 0 11993M * 0 0 * NM:i:0 MD:Z:11993 AS:i:11993 XS:i:11993 11993 33713
 002683K 002683K 0 001532K 1 0 764M * 0 0 * NM:i:0 MD:Z:764 AS:i:764 XS:i:764 764 33364
 002684K 002684K 0 001532K 1 0 764M * 0 0 * NM:i:0 MD:Z:764 AS:i:764 XS:i:764 764 33364
 001658K 001658K 0 001609K 1 0 30956M * 0 0 * NM:i:0 MD:Z:30956 AS:i:30956 XS:i:30956 30956 31821
 001742K 001742K 0 001699K 1 0 29394M * 0 0 * NM:i:0 MD:Z:29394 AS:i:29394 XS:i:29394 29394 30510
 002584K 002584K 0 001928K 1 0 5401M * 0 0 * NM:i:0 MD:Z:5401 AS:i:5401 XS:i:5401 5401 25619
 002602K 002602K 0 002156K 1 0 4871M * 0 0 * NM:i:0 MD:Z:4871 AS:i:4871 XS:i:4871 4871 19333
 002667K 002667K 0 002156K 1 0 1589M * 0 0 * NM:i:0 MD:Z:1589 AS:i:1589 XS:i:1589 1589 19333

 #the first column will now be removed from the final synth fasta
 grep -v  -f scaff.to.remove.onlyscaf finalsynt.2692.len |awk '{print $1}'|cdbyank finalsynt.2692.namemod.nodup.fa.cidx -o dups.removed.finalsynt.2692.nodup.fa
 #scaffold reduction from 2632 to 2587.
 ```

 ### Combining 368 with scaffolds that have missing genes from the above dataset
 ```
scaffolds that do not add genes will be placed in a separate file
##/data021/GIF/remkv6/Baum/CamTechGenomeComparison/14_mergegenomes

#extract and concatenate the genomes
awk '{print $1}' scaf.2692.368.chrlen|cdbyank dups.removed.finalsynt.2692.nodup.fa.cidx -o 2692.genic.subset.fa
cat 2692.genic.subset.fa finalsynt.368.genome.fa >genome.738.fa

#the following scripts use the duplicated removed scaffold list, and removed scaffolds that were added to the genome.738.fa.
grep ">" dups.removed.finalsynt.2692.nodup.fa|sed 's/>//g'| grep -v -f <(awk '{print $1}' scaf.2692.368.chrlen) -|cdbyank dups.removed.finalsynt.2692.nodup.fa.cidx -o repeatscafs.738.fa

#This further removes duplicate scaffolds that were identified in the last iadhore test.
grep ">" repeatscafs.738.fa|sed 's/>//g'| grep -f ../11finalsynth2692selfiadhore/2692_subtraction/scaffold_list_2692_done.list -|cdbyank repeatscafs.738.fa.cidx -o final.repeat.scaffs.738.fa
 ```

 #### Repeat scaffolds assembly statistics
 ```
 Number of scaffolds       1664
Total size of scaffolds   53438982
    Longest scaffold     209206
   Shortest scaffold        213
Number of scaffolds > 1K nt       1655  99.5%
Number of scaffolds > 10K nt       1431  86.0%
Number of scaffolds > 100K nt         26   1.6%
Number of scaffolds > 1M nt          0   0.0%
Number of scaffolds > 10M nt          0   0.0%
  Mean scaffold size      32115
Median scaffold size      29854
 N50 scaffold length      39698
  L50 scaffold count        446
         scaffold %A      31.07
         scaffold %C      18.94
         scaffold %G      18.95
         scaffold %T      31.04
         scaffold %N       0.00
 scaffold %non-ACGTN       0.00
Number of scaffold non-ACGTN nt          0

Percentage of assembly in scaffolded contigs       0.0%
Percentage of assembly in unscaffolded contigs     100.0%
Average number of contigs per scaffold        1.0

Average length of break (>25 Ns) between contigs in scaffold 0

   Number of contigs       1664
Number of contigs in scaffolds          0
Number of contigs not in scaffolds       1664
Total size of contigs   53438982
      Longest contig     209206
     Shortest contig        213
Number of contigs > 1K nt       1655  99.5%
Number of contigs > 10K nt       1431  86.0%
Number of contigs > 100K nt         26   1.6%
Number of contigs > 1M nt          0   0.0%
Number of contigs > 10M nt          0   0.0%
    Mean contig size      32115
  Median contig size      29854
   N50 contig length      39698
    L50 contig count        446
           contig %A      31.07
           contig %C      18.94
           contig %G      18.95
           contig %T      31.04
           contig %N       0.00
   contig %non-ACGTN       0.00
Number of contig non-ACGTN nt          0
 ```
 #### genome 738 statistics
 ```
 Number of scaffolds        738
Total size of scaffolds  124177886
    Longest scaffold    2007323
   Shortest scaffold       3804
Number of scaffolds > 1K nt        738 100.0%
Number of scaffolds > 10K nt        728  98.6%
Number of scaffolds > 100K nt        346  46.9%
Number of scaffolds > 1M nt          7   0.9%
Number of scaffolds > 10M nt          0   0.0%
  Mean scaffold size     168263
Median scaffold size      91698
 N50 scaffold length     303693
  L50 scaffold count        109
         scaffold %A      31.15
         scaffold %C      18.84
         scaffold %G      18.84
         scaffold %T      31.17
         scaffold %N       0.00
 scaffold %non-ACGTN       0.00
Number of scaffold non-ACGTN nt          0

Percentage of assembly in scaffolded contigs       0.0%
Percentage of assembly in unscaffolded contigs     100.0%
Average number of contigs per scaffold        1.0

Average length of break (>25 Ns) between contigs in scaffold 0

   Number of contigs        738
Number of contigs in scaffolds          0
Number of contigs not in scaffolds        738
Total size of contigs  124177886
      Longest contig    2007323
     Shortest contig       3804
Number of contigs > 1K nt        738 100.0%
Number of contigs > 10K nt        728  98.6%
Number of contigs > 100K nt        346  46.9%
Number of contigs > 1M nt          7   0.9%
Number of contigs > 10M nt          0   0.0%
    Mean contig size     168263
  Median contig size      91698
   N50 contig length     303693
    L50 contig count        109
           contig %A      31.15
           contig %C      18.84
           contig %G      18.84
           contig %T      31.17
           contig %N       0.00
   contig %non-ACGTN       0.00
Number of contig non-ACGTN nt          0
 ```
