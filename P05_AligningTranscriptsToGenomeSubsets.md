# Aligning the genome guided transcripts to the new 738 genome, as well as the repeat.738.genome

```
ln -s ../14_mergegenomes/final.repeat.scaffs.738.fa
 ln -s ../14_mergegenomes/genome.738.fa
 ln -s ../../../SplicedLeaders/TrinityGguidedAll/nematode_transcripts_full_Trinity_genome_guided.fasta

creating gmap databases and running
#############################738genome##########################
#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=12:00:00
#PBS -N GMAP-SCN-738
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
module load gmap-gsnap/20160404
gmap_build -d 738.genome  -D /data021/GIF/remkv6/Baum/CamTechGenomeComparison/15_transcriptsTo738/genome.738.fa
gmap -D /data021/GIF/remkv6/Baum/CamTechGenomeComparison/15_transcriptsTo738/ -d 738.genome -B 5 -t 16 --input-buffer-size=1000000 --output-buffer-size=1000000 -f psl --split-output=738genome nematode_transcripts_full_Trinity_genome_guided.fasta

ssh condo "qstat -f ${PBS_JOBID} |head"



#####################repeat.738.genome###################
#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=12:00:00
#PBS -N GMAP-SCN-repeat738
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
#module use /shared/software/GIF/modules
module purge
module load gmap-gsnap/20160404
gmap_build -d repeat.738.genome  -D /data021/GIF/remkv6/Baum/CamTechGenomeComparison/15_transcriptsTo738/final.repeat.scaffs.738.fa
gmap -D /data021/GIF/remkv6/Baum/CamTechGenomeComparison/15_transcriptsTo738/ -d repeat.738.genome -B 5 -t 16 --input-buffer-size=1000000 --output-buffer-size=1000000 -f psl --split-output=repeat.738.genome nematode_transcripts_full_Trinity_genome_guided.fasta

ssh condo "qstat -f ${PBS_JOBID} |head"
```

###  Blasting transcripts to nr database
```
#/data021/GIF/remkv6/Baum/CamTechGenomeComparison/15_transcriptsTo738/nomapIDs
Andrew pulled out this list of gene IDS with the following script
?
1

more GMAP-SCN-738.e110112.condo | grep "No paths" | awk '{print $NF}' > 738genome.nomapIds
softlinked the above file in the new folder
pulling out geneID fasta sequences

cat 738genome.nomapIds |cdbyank nematode_transcripts_full_Trinity_genome_guided.fasta.cidx  -o 738genome.nomapIDS

###PBS Script###
#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=4:00:00
#PBS -N MasterBlaster
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module load parallel
module load ncbi-blast
blastx -db /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/nr/nr -num_threads 16 -max_target_seqs 4 -query 738genome.nomapIDS -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles' -out missing.transcripts.blast.out

# in case you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"

Pulling out those hits that were greater than 80% identical


awk '$3>79.9 {print $1,substr($0,index($0,$14))}' missing.transcripts.blast.out >80PercIDHits.list
Results of this –> only 16/5392 queries had a hit above 80%

TRINITY_GG_438_c78_g1_i1 family protein [Brugia malayi]
TRINITY_GG_438_c78_g1_i1 family protein [Wuchereria bancrofti]
TRINITY_GG_438_c78_g1_i1 protein LOAG_02674 [Loa loa]<>hypothetical protein LOAG_02674 [Loa loa]
TRINITY_GG_438_c78_g1_i1 [Brugia malayi]
TRINITY_GG_438_c78_g3_i1 protein LOAG_02674 [Loa loa]<>hypothetical protein LOAG_02674 [Loa loa]
TRINITY_GG_438_c78_g3_i1 family protein [Brugia malayi]
TRINITY_GG_438_c78_g3_i1 family protein [Wuchereria bancrofti]
TRINITY_GG_438_c78_g3_i1 [Brugia malayi]
TRINITY_GG_438_c78_g2_i1 family protein [Brugia malayi]
TRINITY_GG_438_c78_g2_i1 protein LOAG_02674 [Loa loa]<>hypothetical protein LOAG_02674 [Loa loa]
TRINITY_GG_438_c78_g2_i1 family protein [Wuchereria bancrofti]
TRINITY_GG_438_c78_g2_i1 [Brugia malayi]
TRINITY_GG_438_c139_g1_i2 7 26  2   1371    487 61  336 7e-137    426   Band 7 protein domain containing protein [Haemonchus contortus]
TRINITY_GG_438_c139_g1_i2 7908|emb|CEF62830.1|  87.28   283 15  1   1353    505 173 434 1e-136    429   Band 7 protein family and Stomatin family-containing protein, partial [Strongyloides ratti]
TRINITY_GG_438_c139_g1_i2 protein CELE_F14D12.4 [Caenorhabditis elegans]<>Uncharacterized protein CELE_F14D12.4 [Caenorhabditis elegans]
TRINITY_GG_438_c139_g1_i2 CBR-MEC-2 [Caenorhabditis briggsae]
TRINITY_GG_438_c139_g1_i4 7 61  336 3e-133    417   Band 7 protein domain containing protein [Haemonchus contortus]
TRINITY_GG_438_c139_g1_i4 protein CELE_F14D12.4 [Caenorhabditis elegans]<>Uncharacterized protein CELE_F14D12.4 [Caenorhabditis elegans]
TRINITY_GG_438_c139_g1_i4 7908|emb|CEF62830.1|  84.59   292 15  2   1380    505 173 434 6e-133    420   Band 7 protein family and Stomatin family-containing protein, partial [Strongyloides ratti]
TRINITY_GG_438_c139_g1_i4 CBR-MEC-2 [Caenorhabditis briggsae]
TRINITY_GG_438_c139_g1_i5 7908|emb|CEF62830.1|  81.29   310 36  2   1431    505 146 434 4e-139    431   Band 7 protein family and Stomatin family-containing protein, partial [Strongyloides ratti]
TRINITY_GG_438_c139_g1_i5 briggsae CBR-MEC-2 protein [Caenorhabditis briggsae]
TRINITY_GG_438_c139_g1_i5 7 487 69  336 2e-138    425   Band 7 protein domain containing protein [Haemonchus contortus]
TRINITY_GG_438_c139_g1_i6 7908|emb|CEF62830.1|  82.62   305 36  2   1416    505 146 434 7e-140    433   Band 7 protein family and Stomatin family-containing protein, partial [Strongyloides ratti]
TRINITY_GG_438_c139_g1_i6 briggsae CBR-MEC-2 protein [Caenorhabditis briggsae]
TRINITY_GG_438_c139_g1_i6 7 69  336 5e-139    427   Band 7 protein domain containing protein [Haemonchus contortus]
TRINITY_GG_625_c199_g2_i3 [Loa loa]
TRINITY_GG_625_c199_g2_i3 unc-1 [Toxocara canis]
TRINITY_GG_625_c199_g2_i3 unc-1 [Ascaris suum]
TRINITY_GG_625_c199_g2_i3 protein [Caenorhabditis brenneri]
TRINITY_GG_625_c199_g2_i4 [Loa loa]
TRINITY_GG_625_c199_g2_i4 unc-1 [Ascaris suum]
TRINITY_GG_625_c199_g2_i4 unc-1 [Toxocara canis]
TRINITY_GG_625_c199_g2_i4 protein [Caenorhabditis brenneri]
TRINITY_GG_625_c199_g2_i1 [Loa loa]
TRINITY_GG_625_c199_g2_i1 unc-1 [Ascaris suum]
TRINITY_GG_625_c199_g2_i1 unc-1 [Toxocara canis]
TRINITY_GG_625_c199_g2_i1 protein [Caenorhabditis brenneri]
TRINITY_GG_625_c199_g2_i8 [Loa loa]
TRINITY_GG_625_c199_g2_i8 unc-1 [Ascaris suum]
TRINITY_GG_625_c199_g2_i8 unc-1 [Toxocara canis]
TRINITY_GG_625_c199_g2_i8 protein [Caenorhabditis brenneri]
TRINITY_GG_625_c199_g2_i7 [Loa loa]
TRINITY_GG_625_c199_g2_i7 [Loa loa]
TRINITY_GG_625_c199_g2_i7 unc-1 [Toxocara canis]
TRINITY_GG_625_c199_g2_i7 unc-1 [Toxocara canis]
TRINITY_GG_625_c199_g2_i7 unc-1 [Ascaris suum]
TRINITY_GG_625_c199_g2_i7 unc-1 [Ascaris suum]
TRINITY_GG_625_c199_g2_i7 protein [Caenorhabditis brenneri]
TRINITY_GG_625_c199_g2_i7 protein [Caenorhabditis brenneri]
TRINITY_GG_961_c3_g2_i1 protein PHYSODRAFT_535525, partial [Phytophthora sojae]<>hypothetical protein PHYSODRAFT_535525, partial [Phytophthora sojae]
TRINITY_GG_961_c3_g2_i1 protein SERLA73DRAFT_67532, partial [Serpula lacrymans var. lacrymans S7.3]
TRINITY_GG_961_c3_g2_i1 protein METBIDRAFT_48166 [Metschnikowia bicuspidata var. bicuspidata NRRL YB-4993]<>hypothetical protein HANVADRAFT_28238 [Hanseniaspora valbyensis NRRL Y-1626]
TRINITY_GG_961_c3_g2_i1 protein Bm1_11025 [Brugia malayi]
TRINITY_GG_1117_c3_g1_i1 dehydrogenase, partial [Heterodera glycines]
TRINITY_GG_1309_c1_g1_i1 protein [Caenorhabditis remanei]<>CRE-TAG-308 protein [Caenorhabditis remanei]
TRINITY_GG_1309_c1_g1_i1 [Oesophagostomum dentatum]
TRINITY_GG_1309_c1_g1_i1 protein CAEBREN_19778, partial [Caenorhabditis brenneri]
TRINITY_GG_1309_c1_g1_i1 protein CBG12974 [Caenorhabditis briggsae]<>Protein CBG12974 [Caenorhabditis briggsae]
TRINITY_GG_1353_c35_g1_i1 uncharacterized protein LOC107439839 [Parasteatoda tepidariorum]

Because of this result, we decided to start on the synteny analysis between the JGI genome and our current genome to find missing genes.
```

### Braker run for each genome subset
```
ln -s ../final.repeat.scaffs.738.fa
ln -s ../genome.738.fa



#!/bin/bash #PBS -l nodes=1:ppn=16 #PBS -l walltime=24:00:00 #PBS -N Braker738 #PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID} cd $PBS_O_WORKDIR ulimit -s unlimited module use /shared/software/GIF/modules module purge #not sure why the purge is necessary, but ends up with a perl handshake error if it isnt there

sh ~/common_scripts/runBraker.sh /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R1.fq.gz /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R2.fq.gz /data021/GIF/remkv6/Baum/CamTechGenomeComparison/15_transcriptsTo738/Braker/genome.738.fa

# in case you need stats after job completion retain this as last line ssh condo “qstat -f ${PBS_JOBID} |head”




#!/bin/bash #PBS -l nodes=1:ppn=16 #PBS -l walltime=12:00:00 #PBS -N Braker738 #PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID} cd $PBS_O_WORKDIR ulimit -s unlimited module use /shared/software/GIF/modules module purge #not sure why the purge is necessary, but ends up with a perl handshake error if it isnt there

sh ~/common_scripts/runBraker.sh /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R1.fq.gz /data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/01_Data/20151217_RNAseq_Maker/all_R2.fq.gz /data021/GIF/remkv6/Baum/CamTechGenomeComparison/15_transcriptsTo738/Braker/final.repeat.scaffs.738.fa

# in case you need stats after job completion retain this as last line ssh condo “qstat -f ${PBS_JOBID} |head”

```
### Pertinent stats

Number of genes
```
genome.738.fa
grep -c gene augustus.gff3
31245

final.repeat.scaffs.738.fa
grep -c gene augustus.gff3
14537
```

RNA-seq alignment stats
```
genome.738.fa

301790458 reads; of these:
  301790458 (100.00%) were paired; of these:
    50295631 (16.67%) aligned concordantly 0 times
    163513494 (54.18%) aligned concordantly exactly 1 time
    87981333 (29.15%) aligned concordantly >1 times
    ----
    50295631 pairs aligned concordantly 0 times; of these:
      1682792 (3.35%) aligned discordantly 1 time
    ----
    48612839 pairs aligned 0 times concordantly or discordantly; of these:
      97225678 mates make up the pairs; of these:
        71248041 (73.28%) aligned 0 times
        14997149 (15.43%) aligned exactly 1 time
        10980488 (11.29%) aligned >1 times
88.20% overall alignment rate


final.repeat.scaffs.738.fa

	301790458 reads; of these:
  301790458 (100.00%) were paired; of these:
    239574861 (79.38%) aligned concordantly 0 times
    33926164 (11.24%) aligned concordantly exactly 1 time
    28289433 (9.37%) aligned concordantly >1 times
    ----
    239574861 pairs aligned concordantly 0 times; of these:
      246028 (0.10%) aligned discordantly 1 time
    ----
    239328833 pairs aligned 0 times concordantly or discordantly; of these:
      478657666 mates make up the pairs; of these:
        468685321 (97.92%) aligned 0 times
        4776689 (1.00%) aligned exactly 1 time
        5195656 (1.09%) aligned >1 times
22.35% overall alignment rate

Andrew and I both agree that there are too many genes. Either there is contaminating scaffolds, a pan genome (likely), or there are approximately double the number of genes in HG than in C. elegans (22k). There were 46K in the 2692(368)genome and 24855 in the 368 genome.
```
