# Isoseq analysis of 2 hgtype at three different stages

Email from nico regarding isoseq data
```
The 24857 represented the isoforms called without doing demultiplexing.
The files associated with that analysis would have also been called consensus_isoform.fasta, though they would not have any reference to a specific sample or index.
I have reuploaded the files to the sftp server incase they were not part of the original download,
you will find three files under the Q1055 directory, one simply called consensus_isoforms.fasta which contains all 24857 clusters,
one called polished_high_qv_consensus_isoforms.fasta, which contains quiver polished versions with consensus quality scores greater than .99
and one called polished_low_qv_consensus_isoforms.fasta which contains the reads with a quality lower than .99

Below is a list relating sample name to index number.
idx1 = W82_PA3_Eggs
idx2 = W82_PA3_J2/J3
idx3 = W82_PA3_J4/Adult
idx4 = W82_TN19_Eggs
idx5 = W82_TN19_J4/Adult
idx6 = W82_TN19_J2/J3

As for the methods, first reads of inserts were generated using the pacbio smrtanalysis RS_ReadsOfInsert protocol.
This was done on the pool as a whole, then reads were segregated based on exact matched to the index sequences, this yeilded the idx*_reads_of_insert_sorted.fasta files.
The indexed results then each went into the isoseq protocol,
which will indentify full length transcripts by indentifying the 5' and 3' adapters as well as the polyA tail,
this yeilds the reads_of_insert_forward_and_reverse_idx*_isoseq_flnc.fasta files.
 Next the Full length nonchimeric reads are cluster for correction using the ICE algorithm that is part of the isoseq protocol, and this step creates the reads_of_insert_forward_and_reverse_idx*_isoseq_consensus.fasta files.
We expect a reduction in total reads after each step.  

The oligos for each index are listed below


dT_BC1
tcagacgatgcgtcat
dT_BC2
ctatacatgactctgc
dT_BC3
tactagagtagcactc
dT_BC4
tgtgtatcagtacatg
dT_BC5
gatctctactatatgc
dT_BC6
acagtctatactgctg
```

### Alignment of idx*reads_of_insert_sorted.fasta to consensus_isoforms.fasta
```
Working directory: /data006b/GIF_2b/usha/projects/Baum/01_camtech/09_isoseq_blast/

blast commands


#!/bin/tcsh

#PBS  -o BATCH_OUTPUT.max2hitsblast
#PBS  -e BATCH_ERRORS.max2hitsblast

#PBS -lvmem=256GB,pmem=8Gb,mem=256GB,nodes=1:ppn=64:ib,walltime=48:00:00

# Change to directory from which qsub command was issued
  cd $PBS_O_WORKDIR

module use /data004/software/GIF/modules
module load ncbi-blast/2.2.30+

blastn -db consensus_isoforms.fasta -query idx1_reads_of_insert_sorted.fasta -num_threads 64 -perc_identity 90 -evalue 0.0001 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue frames' > idx1_consensus_blast.txt &

wait

blastn -db consensus_isoforms.fasta -query idx2_reads_of_insert_sorted.fasta -num_threads 64 -perc_identity 90 -evalue 0.0001 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue frames' > idx2_consensus_blast.txt &

wait
blastn -db consensus_isoforms.fasta -query idx3_reads_of_insert_sorted.fasta -num_threads 64 -perc_identity 90 -evalue 0.0001 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue frames' > idx3_consensus_blast.txt &

wait

blastn -db consensus_isoforms.fasta -query idx4_reads_of_insert_sorted.fasta -num_threads 64 -perc_identity 90 -evalue 0.0001 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue frames' > idx4_consensus_blast.txt &

wait

blastn -db consensus_isoforms.fasta -query idx5_reads_of_insert_sorted.fasta -num_threads 64 -perc_identity 90 -evalue 0.0001 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue frames' > idx5_consensus_blast.txt &

wait

blastn -db consensus_isoforms.fasta -query idx6_reads_of_insert_sorted.fasta -num_threads 64 -perc_identity 90 -evalue 0.0001 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue frames' > idx6_consensus_blast.txt &

wait

For presence/absence, get a list of isoforms that are present in each life stage. We used the above blast results and used the below query to extract the ids


awk '{if((($7-$8)/$9 < -.90) && ($9>49)) print $2}' idx6_consensus_blast.txt | sort | uniq -c | sort -r -n -k1 > idx6_isoforms_present.txt

awk '{print $2,$1}' idx1_isoforms_present.txt | sort -k  1,1b> idx1.txt
awk '{print $2,$1}' idx2_isoforms_present.txt | sort -k  1,1b> idx2.txt

join -a1 -a2 -o 0,1.2,2.2 -e "0" idx1.txt idx2.txt | join -a1 -a2 -o 0,1.2,1.3,2.2 -e "0" - idx3.txt  | join -a1 -a2 -o 0,1.2,1.3,1.4,2.2 -e "0" - idx4.txt | join -a1 -a2 -o 0,1.2,1.3,1.4,1.5,2.2 -e "0" - idx5.txt | join -a1 -a2 -o 0,1.2,1.3,1.4,1.5,1.6,2.2 -e "0" - idx6.txt > join.txt
awk '$5==0 && $6==0 && $7==0 && $2!=0 && $3!=0 && $4!=0' join.txt  > only_avir.txt
## 10

awk '$2==0 && $3==0 && $4==0 && $5!=0 && $6!=0 && $7!=0' join.txt > only_vir.txt
## 34
awk '$4==0 && $6!=0' join.txt > adult_noavir_yesvir.txt
## 1644
awk '$3==0 && $7!=0' join.txt > j2j3_noavir_yesvir.txt
## 646
```
### General stats
```
#how many fasta transcripts
grep -c ">" consensus_isoforms.fasta
24857
#number of loci that sequences mapped to
awk '$3=="gene"' genome738sl.polished.consensus_isoforms.gff |cut -f 1,2,3,4,5,6,7,8,9 |sort|uniq|wc
  32569  293121 3526910

#How many of the above fasta sequences mapped to the genome?
awk '$3=="gene"' genome738sl.polished.consensus_isoforms.gff |cut -f 1,2,3,4,5,6,7,8,9 |sort|uniq|cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |awk '{print $1}' |sort |uniq -c |sort -k 1,1 |wc
  22842   45684  491008
Wow, 2005 isoforms did not map to the genome at all...
The same above shows that at max some isoforms mapped 5 times.  

#What is the distribution of multiple mapping?
awk '$3=="gene"' genome738sl.polished.consensus_isoforms.gff |cut -f 1,2,3,4,5,6,7,8,9 |sort|uniq|cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |awk '{print $1}' |sort |uniq -c |awk '{print $1}' |sort |uniq -c |less
 15854 1
   5442 2
    839 3
    221 4
    486 5

#How many of the multiple mapping isoseq isoforms are repeats?
 awk '$3=="gene"' genome738sl.polished.consensus_isoforms.gff |cut -f 1,2,3,4,5,6,7,8,9 |sort|uniq|cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |awk '{print $1}' |sort |uniq -c |awk '$1>1 {print $2}' |grep -f - genome738sl.polished.consensus_isoforms.gff |bedtools intersect -wo -a - -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq |wc
   4862    4862   65764
So, 4862/6988=69.6% have repeats within the gene, but not necessarily exon overlap.

#How many isoforms total, whether multiple mapping or not have repeats within the gene?
bedtools intersect -wo -a <(awk '$3=="gene"' genome738sl.polished.consensus_isoforms.gff) -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq |wc
  15578   15578  211589
So, 15578/32569=  47.8% of the isoforms overlap with a gene that has repeats within, but not necessarily exon overlap.

#How many of these have repeats that overlap exons?
bedtools intersect -wo -a <(awk '$3=="exon"' genome738sl.polished.consensus_isoforms.gff) -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq |wc
  10346   10346  140601
So, 10346/32569= 31.8% of the isoforms overlap with a gene that has a repeat overlapping an exon.

#So I made a list of consensus isoforms that are not repeats, excluding family-976 which seems to play a functional effector role.  
bedtools intersect -wo -a <(awk '$3=="exon"' genome738sl.polished.consensus_isoforms.gff) -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "family-976" |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq |grep -v -f - genome738sl.polished.consensus_isoforms.gff >ConsensusIsoformsNotRepeats

#So how many are there for me to play with?
bedtools intersect -wo -a <(awk '$3=="exon"' genome738sl.polished.consensus_isoforms.gff) -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "family-976" |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq |grep -v -f - genome738sl.polished.consensus_isoforms.gff |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq|wc
  12513   12522  167929
```
### Looking at individual samples
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/6_isoseq/idx1

ln -s /work/GIF4/usha/L3/Baum/01_camtech/09_isoseq_blast/idx1_isoforms_present.txt
ln -s ../ConsensusIsoformsNotRepeats
wget http://purl.obolibrary.org/obo/go.obo
#Performing GO enrichment
 less ConsensusIsoformsNotRepeats |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq |grep -f - idx1_isoforms_present.txt |awk '{print $2}' |grep -f - ../genome738sl.polished.consensus_isoforms.gff |bedtools intersect -wo -a - -b ../../12_functional/augustus_swissprot_iprscan.gff3 |cut -f 18- |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sed 's/;/\t/g' |grep -v "Parent" |cut -f 1 |sort|uniq >idx1.gene.list
ln -s ../../../57_secretome/ontologenizer/simpleformat.ids
ln -s ../../../57_secretome/ontologenizer/population
cp ../../../57_secretome/ontologenizer/Ontologizer.jar .

java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s idx1.gene.list

There were no significant go terms enriched when using the p adjusted values.



#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/6_isoseq/idx2

ln -s /work/GIF4/usha/L3/Baum/01_camtech/09_isoseq_blast/idx2_isoforms_present.txt
ln -s ../ConsensusIsoformsNotRepeats
wget http://purl.obolibrary.org/obo/go.obo
#Performing GO enrichment
 less ConsensusIsoformsNotRepeats |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq |grep -f - idx2_isoforms_present.txt |awk '{print $2}' |grep -f - ../genome738sl.polished.consensus_isoforms.gff |bedtools intersect -wo -a - -b ../../12_functional/augustus_swissprot_iprscan.gff3 |cut -f 18- |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sed 's/;/\t/g' |grep -v "Parent" |cut -f 1 |sort|uniq >idx2.gene.list
ln -s ../../../57_secretome/ontologenizer/simpleformat.ids
ln -s ../../../57_secretome/ontologenizer/population
cp ../../../57_secretome/ontologenizer/Ontologizer.jar .

java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s idx2.gene.list


Some significant go terms were found
ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0043604      10960   213     489     26      1069    39      2       false   1.1535438131965266E-10  9.078389809856664E-8    4.7911230994740565E-231 "amide biosynthetic process"
GO:0043603      10960   229     489     27      1435    49      1       false   1.4480709243805636E-10  1.1396318174875036E-7   9.249512420609186E-273  "cellular amide metabolic process"
GO:0006575      10960   88      489     15      2873    84      2       false   1.5755861612994052E-8   1.2399863089426318E-5   3.3066636862877092E-170 "cellular modified amino acid metabolic process"
GO:0006518      10960   218     489     26      1658    64      2       false   1.9126268094072635E-8   1.5052372990035163E-5   1.9599105144441444E-279 "peptide metabolic process"
GO:1901566      10960   401     489     32      2429    74      2       false   3.1107290799885825E-8   2.4481437859510143E-5   0.0     "organonitrogen compound biosynthetic process"
GO:0006414      10960   127     489     12      917     18      2       false   2.748845365017362E-7    2.1633413022686637E-4   1.7167452646656717E-159 "translational elongation"
GO:0016874      10960   175     489     20      3076    108     1       false   1.5785208553727008E-6   0.0012422959131783155   6.994592154040206E-291  "ligase activity"
GO:0006790      10960   126     489     15      2376    76      1       false   5.788458198252365E-6    0.0045555166020246115   3.0475658379282293E-213 "sulfur compound metabolic process"
GO:0044272      10960   108     489     14      1208    43      2       false   7.891355408834837E-6    0.006210496706753017    2.515737838301016E-157  "sulfur compound biosynthetic process"
GO:0044763      10960   786     489     41      3559    99      2       false   9.392076688689725E-6    0.0073915643539988135   0.0     "single-organism cellular process"
GO:0006082      10960   158     489     16      3214    99      4       false   1.768522996295724E-5    0.013918275980847349    7.223121488300162E-273  "organic acid metabolic process"
GO:0016774      10960   12      489     4       854     15      1       false   2.8245350450576936E-5   0.02222909080460405     3.439995984260115E-27   "phosphotransferase activity, carboxyl group as acceptor"
GO:0044267      10960   895     489     32      2393    48      2       false   3.136710511432141E-5    0.02468591172497095     0.0     "cellular protein metabolic process"
GO:1901265      10960   1566    489     60      3206    86      2       false   5.569463578219476E-5    0.043831678360587274    0.0     "nucleoside phosphate binding"



#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/6_isoseq/idx3

ln -s /work/GIF4/usha/L3/Baum/01_camtech/09_isoseq_blast/idx3_isoforms_present.txt
ln -s ../ConsensusIsoformsNotRepeats
wget http://purl.obolibrary.org/obo/go.obo
#Performing GO enrichment
 less ConsensusIsoformsNotRepeats |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq |grep -f - idx3_isoforms_present.txt |awk '{print $2}' |grep -f - ../genome738sl.polished.consensus_isoforms.gff |bedtools intersect -wo -a - -b ../../12_functional/augustus_swissprot_iprscan.gff3 |cut -f 18- |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sed 's/;/\t/g' |grep -v "Parent" |cut -f 1 |sort|uniq >idx3.gene.list
ln -s ../../../57_secretome/ontologenizer/simpleformat.ids
ln -s ../../../57_secretome/ontologenizer/population
cp ../../../57_secretome/ontologenizer/Ontologizer.jar .

java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s idx3.gene.list

#there were a few significatnly enriched go terms
ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0003824      10970   3076    578     152     6526    228     1       false   1.1401627194158623E-9   1.0535103527402566E-6   0.0     "catalytic activity"
GO:0006082      10970   158     578     19      3214    118     4       false   2.8795731885141067E-6   0.0026607256261870346   7.223121488300162E-273  "organic acid metabolic process"
GO:0044710      10970   821     578     55      4142    156     2       false   3.1452519781148037E-6   0.0029062128277780788   0.0     "single-organism metabolic process"
GO:0043167      10970   1704    578     80      4351    141     1       false   1.3385238289602678E-5   0.012367960179592874    0.0     "ion binding"



#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/6_isoseq/idx4

ln -s /work/GIF4/usha/L3/Baum/01_camtech/09_isoseq_blast/idx4_isoforms_present.txt
ln -s ../ConsensusIsoformsNotRepeats
wget http://purl.obolibrary.org/obo/go.obo
#Performing GO enrichment
 less ConsensusIsoformsNotRepeats |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq |grep -f - idx4_isoforms_present.txt |awk '{print $2}' |grep -f - ../genome738sl.polished.consensus_isoforms.gff |bedtools intersect -wo -a - -b ../../12_functional/augustus_swissprot_iprscan.gff3 |cut -f 18- |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sed 's/;/\t/g' |grep -v "Parent" |cut -f 1 |sort|uniq >idx4.gene.list
ln -s ../../../57_secretome/ontologenizer/simpleformat.ids
ln -s ../../../57_secretome/ontologenizer/population
cp ../../../57_secretome/ontologenizer/Ontologizer.jar .

java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s idx4.gene.list

There were some significant GO terms
ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:1901265      11277   1566    1361    157     3206    238     2       false   2.3904488944275943E-8   3.234277354160535E-5    0.0     "nucleoside phosphate binding"
GO:0005575      11277   2016    1361    215     7387    593     1       false   3.9832563470057527E-7   5.389345837498783E-4    0.0     "cellular_component"
GO:0016875      11277   56      1361    17      175     21      1       false   1.1784027376654048E-6   0.0015943789040612927   3.5248560608876664E-47  "ligase activity, forming carbon-oxygen bonds"
GO:0005315      11277   8       1361    5       68      5       1       false   5.372151991993865E-6    0.0072685216451677      1.352812065169356E-10   "inorganic phosphate transmembrane transporter activity"
GO:0003824      11277   3076    1361    284     6526    506     1       false   1.5239154662360055E-5   0.020618576258173153    0.0     "catalytic activity"
GO:0016903      11277   30      1361    12      380     43      1       false   1.896083328243544E-5    0.025654007431135148    3.4729115893368666E-45  "oxidoreductase activity, acting on the aldehyde or oxo group of donors"
GO:0043167      11277   1704    1361    174     4351    352     1       false   2.9070266249236015E-5   0.03933207023521633     0.0     "ion binding"
GO:0097367      11277   1081    1361    120     4351    352     1       false   3.1299738691273826E-5   0.042348546449293485    0.0     "carbohydrate derivative binding"
GO:0045184      11277   106     1361    25      764     83      2       false   3.349753693696343E-5    0.04532216747571152     5.907339830536063E-133  "establishment of protein localization"
GO:0006082      11277   158     1361    28      3214    258     4       false   3.643631927858097E-5    0.04929833998392005     7.223121488300162E-273  "organic acid metabolic process"




#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/6_isoseq/idx5

ln -s /work/GIF4/usha/L3/Baum/01_camtech/09_isoseq_blast/idx5_isoforms_present.txt
ln -s ../ConsensusIsoformsNotRepeats
wget http://purl.obolibrary.org/obo/go.obo
#Performing GO enrichment
 less ConsensusIsoformsNotRepeats |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq |grep -f - idx5_isoforms_present.txt |awk '{print $2}' |grep -f - ../genome738sl.polished.consensus_isoforms.gff |bedtools intersect -wo -a - -b ../../12_functional/augustus_swissprot_iprscan.gff3 |cut -f 18- |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sed 's/;/\t/g' |grep -v "Parent" |cut -f 1 |sort|uniq >idx5.gene.list
ln -s ../../../57_secretome/ontologenizer/simpleformat.ids
ln -s ../../../57_secretome/ontologenizer/population
cp ../../../57_secretome/ontologenizer/Ontologizer.jar .

java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s idx5.gene.list     

ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0006414      11064   127     804     16      917     32      2       false   6.46769227868466E-7     7.289089198077612E-4    1.7167452646656717E-159 "translational elongation"
GO:1901265      11064   1566    804     90      3206    131     2       false   2.1594743110101664E-6   0.0024337275485084575   0.0     "nucleoside phosphate binding"
GO:0044267      11064   895     804     64      2393    111     2       false   7.316952913907738E-6    0.008246205933974021    0.0     "cellular protein metabolic process"
GO:0006082      11064   158     804     22      3214    167     4       false   1.3868917633764027E-5   0.01563027017325206     7.223121488300162E-273  "organic acid metabolic process"
GO:0044710      11064   821     804     66      4142    207     2       false   1.6705041627660293E-5   0.018826581914373152    0.0     "single-organism metabolic process"



#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/6_isoseq/idx6

ln -s /work/GIF4/usha/L3/Baum/01_camtech/09_isoseq_blast/idx6_isoforms_present.txt
ln -s ../ConsensusIsoformsNotRepeats
wget http://purl.obolibrary.org/obo/go.obo
#Performing GO enrichment
 less ConsensusIsoformsNotRepeats |cut -f 9 |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sort |uniq |grep -f - idx6_isoforms_present.txt |awk '{print $2}' |grep -f - ../genome738sl.polished.consensus_isoforms.gff |bedtools intersect -wo -a - -b ../../12_functional/augustus_swissprot_iprscan.gff3 |cut -f 18- |sed 's/ID=//g' |sed 's/\./\t/g' |cut -f 1 |sed 's/;/\t/g' |grep -v "Parent" |cut -f 1 |sort|uniq >idx6.gene.list
ln -s ../../../57_secretome/ontologenizer/simpleformat.ids
ln -s ../../../57_secretome/ontologenizer/population
cp ../../../57_secretome/ontologenizer/Ontologizer.jar .

java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s idx6.gene.list


ID      Pop.total       Pop.term        Study.total     Study.term      Pop.family      Study.family    nparents        is.trivial      p       p.adjusted      p.min   name
GO:0043603      10920   229     407     22      1435    33      1       false   4.442130189717575E-11   2.958458706351905E-8    9.249512420609186E-273  "cellular amide metabolic process"
GO:0006414      10920   127     407     11      917     13      2       false   1.4862420779511826E-8   9.898372239154875E-6    1.7167452646656717E-159 "translational elongation"
GO:0043604      10920   213     407     18      1069    25      2       false   1.6769689706077117E-8   1.116861334424736E-5    4.7911230994740565E-231 "amide biosynthetic process"
GO:0006518      10920   218     407     19      1658    48      2       false   2.9291402644045446E-6   0.0019508074160934268   1.9599105144441444E-279 "peptide metabolic process"
GO:1901566      10920   401     407     22      2429    52      2       false   7.615237827306782E-6    0.005071748392986317    0.0     "organonitrogen compound biosynthetic process"
GO:0006412      10920   127     407     11      1809    34      7       false   9.485998950685118E-6    0.006317675301156288    5.642801043472887E-199  "translation"
GO:0005737      10920   526     407     27      1165    33      1       false   1.336329531777886E-5    0.008899954681640721    0.0     "cytoplasm"
GO:0034660      10920   92      407     8       687     12      1       false   2.4479889629806258E-5   0.01630360649345097     7.361424346186463E-117  "ncRNA metabolic process"
GO:0030149      10920   4       407     3       93      3       3       false   3.082471525669198E-5    0.020529260360956855    3.4249683618547254E-7   "sphingolipid catabolic process"
GO:0006413      10920   127     407     11      2376    52      3       false   6.104089608645962E-5    0.040653236793582105    1.7201816062989902E-214 "translational initiation"
GO:0044763      10920   786     407     30      3559    70      2       false   7.017322995874808E-5    0.046735371152526226    0.0     "single-organism cellular process"
```
