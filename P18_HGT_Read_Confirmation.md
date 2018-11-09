# Need to confirm that the horizontal gene transfers are indeed integral to the contigs, and not just on the edges.  If so, how do they overlap with subread alignments?

### Fix subreads alignment with blasr to gff format
```
Matt had some doubts on the large number of HGT candidates (1092), and suggested to find pacbio reads that overlapped the HGT gene and an adjacent gene to be sure that we were not including contaminants in our assembly.


#/work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck/04_BlastPacBioReads
less ../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.out |awk '{if($6==0) {print $2,"Blasr","Subread","1",$7,"+",".",".",$1}else {print $2,"Blasr","Subread",$6,$7,"+",".",".",$1}' |tr " " "\t" >../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff

less ../../36_deduplicate122016/preadMapping/blasr/preads2genome.out |awk '{if($6==0) {print $2,"Blasr","pread","1",$7,"+",".",".",$1}else {print $2,"Blasr","pread",$6,$7,"+",".",".",$1}}' |tr " " "\t" >../../36_deduplicate122016/preadMapping/preads2genome.gff3

less ../../36_deduplicate122016/ccsMapping/blasrAlignment/ccs2genome.out |awk '{if($6==0) {print $2,"Blasr","ccsread","1",$7,"+",".",".",$1}else {print $2,"Blasr","ccsread",$6,$7,"+",".",".",$1}}' |tr " " "\t" >../../36_deduplicate122016/ccsMapping/blasrAlignment/ccs2genome.gff
```

### Check both adjacent genes to see if they are also HGT
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck

#Generate a gff that includes the gene 5' and gene 3' of the HGT event.
sort -V HGTMod.list |while read line; do grep --no-group-separator -A 1 -B 1 -w "$line" augustus.mod.gff3;done >HGTNSurroundingGenes.gff

paste <(awk -v x=1 'NR==x {print $0;x=x+3}'  HGTNSurroundingGenes.gff) <(awk -v x=3 'NR==x {print $0;x=x+3}'  HGTNSurroundingGenes.gff) |awk '{print $1,$2,$3,$4,$14,$15,$16,$17,$18}' |sed 's/g/g\t/2' |awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9$10-1}' |tr " " "\t" >HGTSurroundingGenes.gff

#how many genes are at the ends of scaffolds and generated an invalid gff line?
awk '$4>$5' HGTSurroundingGenes.gff >HGTEndsOfScaf.gff
69
#how many were not at the ends of scaffolds?
awk '$5>$4' HGTSurroundingGenes.gff >HGTNotEndsOfScaf.gff
1021

#how many of these 1021 have a subread overlapping the HGT gene and both adjacent genes at 100%
bedtools intersect -f 1 -wo -a HGTNotEndsOfScaf.gff  -b ../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |awk '{print $9}'|sort|uniq|wc
    118     118     769

"Only 118 HGT genes had a pacbio subread completely overlapping an HGT gene and both adjacent genes."

#how many of the 118 are in the high confidence list?
bedtools intersect -f 1 -wo -a HGTNotEndsOfScaf.gff  -b ../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |awk '{print $9}'|sort|uniq|grep -w -f - HighConfidence.list |wc
     20      20     131
"20 of the high confidence HGT events had a pacbio subreat completely overlapping an HGT gene and both adjacent genes"


#how many of the high confidence HGT events are at the ends of scaffolds
awk '{print $9}' HGTEndsOfScaf.gff |grep -w -f - HighConfidence.list |wc
      6       6      41

#how many of the HGT events that are novel to SCN are found at the ends of scaffolds?      
awk '{print $9}' HGTEndsOfScaf.gff |grep -w -f - NovelSCNHGTFormatted.list |wc
```

### 5' adjacent gene and HGT gene
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck/FivePrimeGeneHGTTest

#5' adjacent gene and HGT gene
sort -V ../HGTMod.list |while read line; do grep --no-group-separator -B 1 -w "$line" ../augustus.mod.gff3;done >HGTNSurroundingGenes.gff

paste <(awk -v x=1 'NR==x {print $0;x=x+2}'  HGTNSurroundingGenes.gff) <(awk -v x=2 'NR==x {print $0;x=x+2}'  HGTNSurroundingGenes.gff) |awk '{print $1,$2,$3,$4,$14,$15,$16,$17,$18}' |sed 's/g/g\t/2' |awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9$10}' |tr " " "\t" >HGTSurroundingGenes.gff

#removing genes at the ends of scaffolds
awk '$5>$4' HGTSurroundingGenes.gff >HGTNotEndsOfScaf.gff
1057

#How many HGT and their 5' adjacent gene are completely overlapped by at least one pacbio subread?
 bedtools intersect -f 1 -wo -a HGTNotEndsOfScaf.gff  -b ../../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |awk '{print $9}'|sort|uniq|wc
     159     159    1048
```

### 3' adjacent gene and HGT gene
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck/ThreePrimeGeneHGTTest

#3' adjacent gene and HGT gene
sort -V ../HGTMod.list |while read line; do grep --no-group-separator -A 1 -w "$line" ../augustus.mod.gff3;done >HGTNSurroundingGenes.gff

paste <(awk -v x=1 'NR==x {print $0;x=x+2}'  HGTNSurroundingGenes.gff) <(awk -v x=2 'NR==x {print $0;x=x+2}'  HGTNSurroundingGenes.gff) |awk '{print $1,$2,$3,$4,$14,$15,$16,$17,$18}' |sed 's/g/g\t/2' |awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9$10-1}' |tr " " "\t" >HGTSurroundingGenes.gff

#removing genes at the ends of scaffolds
awk '$5>$4' HGTSurroundingGenes.gff >HGTNotEndsOfScaf.gff
1056

bedtools intersect -f 1 -wo -a HGTNotEndsOfScaf.gff  -b ../../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |awk '{print $9}'|sort|uniq|wc
    153     153    1008
```

### Combination of both 5' and 3'
```
cat <(bedtools intersect -f 1 -wo -a HGTNotEndsOfScaf.gff  -b ../../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |awk '{print $9}'|sort|uniq) <(bedtools intersect -f 1 -wo -a ../ThreePrimeGeneHGTTest/HGTNotEndsOfScaf.gff  -b ../../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |awk '{print $9}'|sort|uniq)|sort|uniq|wc
    192     192    1274
cat <(bedtools intersect -f 1 -wo -a HGTNotEndsOfScaf.gff  -b ../../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |awk '{print $9}'|sort|uniq) <(bedtools intersect -f 1 -wo -a ../ThreePrimeGeneHGTTest/HGTNotEndsOfScaf.gff  -b ../../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |awk '{print $9}'|sort|uniq)|sort|uniq >EitherOverlap5prime3prime.list

How many of these genes have introns?
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck/01_IntronicHGT
cd FivePrimeGeneHGTTest/
 cat <(bedtools intersect -f 1 -wo -a HGTNotEndsOfScaf.gff  -b ../../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |awk '{print $9}'|sort|uniq) <(bedtools intersect -f 1 -wo -a ../ThreePrimeGeneHGTTest/HGTNotEndsOfScaf.gff  -b ../../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |awk '{print $9}'|sort|uniq)|sort|uniq >../01_IntronicHGT/OverlappedHGT.list
cd ../01_IntronicHGT/
grep -w -f OverlappedHGT.list augustus.exons.gff3 >OverlappedHGTExon.gff3
wc OverlappedHGT.list
 192  192 1274 OverlappedHGT.list
ess OverlappedHGTExon.gff3|grep -v "exon1"|awk '{print $9}' |sort|uniq|wc
    158     158    1051
"878/1092 and  158/192 total genes and read-confirmed genes had at least one intron"    

grep -w -f ../HGTMod.list augustus.exons.gff3 |grep -v "exon1"|awk '{print $9}' |sort|uniq|wc
    878     878    5660
"878/1092 and  158/192 total genes and read-confirmed genes had at least one intron"    
```

### What are the blast hits to these genes?
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck/02_HGTBLASTNR
grep -w -f HGTMod.list <(awk '$3=="mRNA"' ../../52_functional/augustus_swissprot_iprscan.gff3 |sed 's/\.t/\t/g' |sed 's/ID=//g') >HGTFunctional.gff
 cut -f 9- HGTFunctional.gff |grep -v "Caenorhabditis" |grep -v "Trypanosoma" |grep -v "Onchocerca" |grep -v "Haemonchus" |grep -v "Ascaris" |awk '{print $1}' |sort|uniq |wc
    854     854    5506


#854 of the 1092 HGT genes did not have a isoform top hit to another nematode, although Benjamin has another list of these.
```
### Apply these filters to the Nematoda genes that we are claiming are novel
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck/03_NovelHGT
grep -w -f NovelSCNHGTFormatted.list ../FivePrimeGeneHGTTest/EitherOverlap5prime3prime.list |wc
     19      19     130
This if I take the remaining 73 and remove those that have an intron how many are left.
grep -w -v -f ../01_IntronicHGT/OverlappedHGTExon.list   NovelSCNHGTFormatted.list |sort|uniq |grep -v -w -f ../FivePrimeGeneHGTTest/EitherOverlap5prime3prime.list - |wc
73
"None of these genes have an intron."

#how many had blast hits that were to nematodes in our dataset?
grep -w -v -f ../01_IntronicHGT/OverlappedHGTExon.list   NovelSCNHGTFormatted.list |sort|uniq|grep -v -w -f ../02_HGTBLASTNR/NotNematoda.list - |wc
     11      11      72
#how many does that leave     
grep -w -v -f ../01_IntronicHGT/OverlappedHGTExon.list   NovelSCNHGTFormatted.list |sort|uniq|grep -w -f ../02_HGTBLASTNR/NotNematoda.list - |wc
     81      81     521

None of the three datasets support 28/92 of the NovelSCN
grep -w -v -f ../01_IntronicHGT/OverlappedHGTExon.list   NovelSCNHGTFormatted.list |sort|uniq |grep -v -w -f ../FivePrimeGeneHGTTest/EitherOverlap5prime3prime.list - |grep -w -f ../02_HGTBLASTNR/NotNematoda.list |wc
     64      64     405
```

### Apply these filters to the Nematoda genes that we are claiming are novel
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck/03_NovelHGT
grep -w -f NovelSCNHGTFormatted.list ../FivePrimeGeneHGTTest/EitherOverlap5prime3prime.list |wc
     19      19     130
This if I take the remaining 73 and remove those that have an intron how many are left.
grep -w -v -f ../01_IntronicHGT/OverlappedHGTExon.list   NovelSCNHGTFormatted.list |sort|uniq |grep -v -w -f ../FivePrimeGeneHGTTest/EitherOverlap5prime3prime.list - |wc
73
"None of these genes have an intron."

#how many had blast hits that were to nematodes in our dataset?
grep -w -v -f ../01_IntronicHGT/OverlappedHGTExon.list   NovelSCNHGTFormatted.list |sort|uniq|grep -v -w -f ../02_HGTBLASTNR/NotNematoda.list - |wc
     11      11      72
#how many does that leave     
grep -w -v -f ../01_IntronicHGT/OverlappedHGTExon.list   NovelSCNHGTFormatted.list |sort|uniq|grep -w -f ../02_HGTBLASTNR/NotNematoda.list - |wc
     81      81     521

None of the three datasets support 28/92 of the NovelSCN
grep -w -v -f ../01_IntronicHGT/OverlappedHGTExon.list   NovelSCNHGTFormatted.list |sort|uniq |grep -v -w -f ../FivePrimeGeneHGTTest/EitherOverlap5prime3prime.list - |grep -w -f ../02_HGTBLASTNR/NotNematoda.list |wc
     64      64     405
```

### What are the reads that overlap these genes?

```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck/04_BlastPacBioReads
How many reads overlap these genes?
<sxh>


grep -w -f UnsupportedHGTs.list <(awk '$3=="gene"' ../augustus.mod.gff3 ) >UnsupportedHGTs.gff3

bedtools intersect  -wo -a UnsupportedHGTs.gff3 -b ../../36_deduplicate122016/subreadMapping/blasrAlignment/subreads2genome.gff |wc
      0       0       0
bedtools intersect  -wo -a UnsupportedHGTs.gff3 -b ../../36_deduplicate122016/preadMapping/preads2genome.gff3 |wc
      0       0       0
bedtools intersect  -wo -a UnsupportedHGTs.gff3 -b ../../36_deduplicate122016/ccsMapping/blasrAlignment/ccs2genome.gff |wc
      0       0       0
```
### Map genes to reads to find missing unmapped reads
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck/05_BLASTUnsupportedGenes
module load cdbfasta
cdbfasta gff2fasta.gene.fasta
ln -s ../04_BlastPacBioReads/UnsupportedHGTs.list  

less UnsupportedHGTs.list| cdbyank gff2fasta.gene.fasta.cidx >UnsupportedHGTs.fasta
grep -c ">" UnsupportedHGTs.fasta
64


#copy and modify below script to run on a large database.gmapl
cp ~/common_scripts/runGmap.sh .
vi runGmap.sh
sh gmap.sh
sh runGmap.sh Subreads /work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck/05_BLASTUnsupportedGenes/  SCN.all.subreads.sl.fasta UnsupportedHGTs.fasta

#all genes mapped to a read.
less Subreads.UnsupportedHGTs.psl |awk '{print $10}' |sort|uniq|wc
     64      64     405

#extract reads that genes mapped to.
cdbfasta SCN.all.subreads.sl.fasta
awk '{print $14}' Subreads.UnsupportedHGTs.psl |sort|uniq |cdbyank SCN.all.subreads.sl.fasta.cidx >UnsupportedHGTReads.fasta

#64 genes mapped to 201 reads.
grep -c ">" UnsupportedHGTReads.fasta
201

#modify nt location, and removed megablast and evalue threshold.
cp ~/common_analyses/runMegablast.sh .

sh runMegablast.sh UnsupportedHGTReads.fasta &
#so 8 reads have blast hits in nt greater than a score of 200
less UnsupportedHGTReads.vs.nt.cul5.1e5.megablast.out |awk '$3>200' |awk '{print $1}' |sort|uniq|wc
8

#all 8 of these blast hits are to Heterodera species.

#This leaves us with 56 possible HGT events that are unverified.
less UnsupportedHGTReads.vs.nt.cul5.1e5.megablast.out |awk '$3>200' |cut -f 1 |grep -w -f - Subreads.UnsupportedHGTs.psl|awk '{print $10 }' |grep -v -w -f - UnsupportedHGTs.list >UnsupportedHGTsMinusBLASTHits.list
wc UnsupportedHGTsMinusBLASTHits.list
 56  56 354

#How many were called as orthologs to another nematode in my analyses?
grep -w -f  UnsupportedHGTsMinusBLASTHits.list ../../62_totalOrthologues/SCN.all.orthologues.list |wc
     47      47     297

#How many does that leave that are unverified.
[remkv6@condo012 05_BLASTUnsupportedGenes]$ grep -v -w -f   ../../62_totalOrthologues/SCN.all.orthologues.list UnsupportedHGTsMinusBLASTHits.list |wc
      9       9      57

#If I directly blast these genes to nt, what is the best hit for each?
less UnsupportedHGTsMinusBLASTHitsNotOrtholog.list |cdbyank gff2fasta.gene.fasta.cidx >UnsupportedHGTsMinusBLASTHitsNotOrtholog.fasta
sh runMegablast.sh UnsupportedHGTsMinusBLASTHitsNotOrtholog.fasta
wc UnsupportedHGTsMinusBLASTHitsNotOrtholog.vs.nt.cul5.1e5.megablast.out
  1  24 185
#one hit to H. glycines cold shock gene, of about 118bp to a 700bp gene.  
#How many unverified genes left?
grep -v -w -f <(awk '{print $1}' UnsupportedHGTsMinusBLASTHitsNotOrtholog.vs.nt.cul5.1e5.megablast.out) UnsupportedHGTsMinusBLASTHitsNotOrtholog.list |sed '/^$/d' |wc
      8       8      50

grep -v -w -f <(awk '{print $1}' UnsupportedHGTsMinusBLASTHitsNotOrtholog.vs.nt.cul5.1e5.megablast.out) UnsupportedHGTsMinusBLASTHitsNotOrtholog.list |sed '/^$/d' >UnsupportedHGTsMinusBLASTHitsNotOrthologGeneBlast.list
less UnsupportedHGTsMinusBLASTHitsNotOrthologGeneBlast.list |awk '{print $1".t1"}' |cdbyank augustus.aa.cidx >UnsupportedHGTsMinusBLASTHitsNotOrthologGeneBlast.fasta


#running a tblastn on the amino acids to see if I can get a hit.
sh runMegablast.sh UnsupportedHGTsMinusBLASTHitsNotOrthologGeneBlast.fasta

#how many had a reasonable nematode hit?
less UnsupportedHGTsMinusBLASTHitsNotOrthologGeneBlast.vs.nt.cul5.1e5.megablast.out|awk '$3>100 {print $1}' |sort|uniq|wc
      3       3      28

#How many are left?
 less UnsupportedHGTsMinusBLASTHitsNotOrthologGeneBlast.vs.nt.cul5.1e5.megablast.out|awk '$3>100 {print $1}' |sed 's/\.t/\t/g' |awk '{print $1}' |sort|uniq|grep -w -v -f - UnsupportedHGTsMinusBLASTHitsNotOrthologGeneBlast.list|wc
      5       5      31
less UnsupportedHGTsMinusBLASTHitsNotOrthologGeneBlast.vs.nt.cul5.1e5.megablast.out|awk '$3>100 {print $1}' |sed 's/\.t/\t/g' |awk '{print $1}' |sort|uniq|grep -w -v -f - UnsupportedHGTsMinusBLASTHitsNotOrthologGeneBlast.list >UnsupportedHGTsMinusBLASTHitsNotOrthologGeneBlastProtBlast.list

#going to extract the reads that these genes map to and blastn, without a megablast

sh runMegablast.sh UnsupportedHGTsMinusBLASTHitsNotOrthologGeneBlastProtBlastReads.fasta

#none were able to show a reasonably good blast hit.  
```

### HGT Novel to nematodes vs other genomic strata
```
#make high confidence HGT gene list
#copied list from Benjamin's excel file.
less HGT.txt |tr " " "\n" |sed '/^$/d' >HGTOldNames.list
sed -i 's/augustus.aa/HGTOldNames.list/g' geneRenamer.sh
sh geneRenamer.sh &
cp HGTOldNames.list HGTNewNames.list

#how many of the 156 high confidence HGT are known effectors (gmap alignment)
sed 's/\./\t/2' HGTOldNames.list|awk '{print $1}' |cat <(awk '{print $2}' ../../58_Renamatorium/18_effectorRedo/121EffectorOldNewGeneNames.list) - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
      0       0       0
#how many are syntenic?
 sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/GenesInSynteny.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
     27      54     702

#how many HGT genes are found in tandem duplications
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/tandem.gene.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
     16      32     416

#How many are found in the effectorome?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/Intermixed.network - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
     10      30     381
#How many are secreted?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/Secretome.gene.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
      8      16     208
#how many are orthologous to a gene in a related species?
sed 's/\./\t/1' HGTOldNames.list|awk '{print $1}' |cat ../../62_totalOrthologues/SCN.all.orthologues.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
     77     154    1114
#How many are repeat-affected?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../72_RepeatNTandemFinalGFF/RepeatAffectedGenesNoSimple.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
     47      94    1222
#How many are affected by a family1265 LINE?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/Family1265.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
      1       2      26
#How many are affected by a family976 repeat?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/Family976.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
      0       0       0
How many are found in LTR retrotransposons?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat <( bedtools intersect -wo -a ../../58_Renamatorium/23_LTR_finder/NonRedundantLtrRetroelementMerge.gff -b ../../58_Renamatorium/1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq )  - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
      0       0       0
How many are found in DNA transposons?
  sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat <( bedtools intersect -wo -a ../../58_Renamatorium/24_IRF_DNATrans/SupportedIRFMergeClassified.gff -b ../../58_Renamatorium/1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq|wc )  - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
      0       0       0
```


### HGT 147 high confidence vs genomic strata
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/76_HGTCheck/07_HighConf156
mkdir 07_HighConf156
cd 07_HighConf156/
vi HGTHighOldNames.list
cp ../../58_Renamatorium/1_genomeNgff/geneRenamer.sh .
sed -i 's/augustus.aa/HGTHighOldNames.list/g' geneRenamer.sh
sh geneRenamer.sh &
 cp HGTHighOldNames.list HGTNewNames.list
 vi HGTHighOldNames.list
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |sort|uniq >fixed
mv fixed HGTNewNames.list

#how many of the 156 high confidence HGT are known effectors (gmap alignment)
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat <(awk '{print $2}' ../../58_Renamatorium/18_effectorRedo/121EffectorOldNewGeneNames.list) - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
      1       2      26
#which effector is it?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat <(awk '{print $2}' ../../58_Renamatorium/18_effectorRedo/121EffectorOldNewGeneNames.list) - |sort|uniq -c |sort -k1,1nr|awk '$1==2 {print $2}' |grep -w -f - ../../58_Renamatorium/18_effectorRedo/4SebastianEffvsGenes.list |less
Hetgly.G000016560 GLAND12 Pioneer


#how many are syntenic?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/GenesInSynteny.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
     39      78    1014

#how many HGT genes are found in tandem duplications
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/tandem.gene.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
     26      52     676

#How many are found in the effectorome?
 sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/Intermixed.network - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
     10      30     381

#How many are secreted?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/Secretome.gene.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
     14      28     364


#how many are orthologous to a gene in a related species?
sed 's/\./\t/1' HGTHighOldNames.list|awk '{print $1}' |sort|uniq |cat ../../62_totalOrthologues/SCN.all.orthologues.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
    120     240    1731


#How many are repeat-affected?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../72_RepeatNTandemFinalGFF/RepeatAffectedGenesNoSimple.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
     67     134    1742


#How many are affected by a family1265 LINE?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/Family1265.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
      1       2      26

#How many are affected by a family976 repeat?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat ../../58_Renamatorium/20_Expression/Family976.list - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
      0       0       0


How many are found in LTR retrotransposons?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat <( bedtools intersect -wo -a ../../58_Renamatorium/23_LTR_finder/NonRedundantLtrRetroelementMerge.gff -b ../../58_Renamatorium/1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq )  - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
      0       0       0

How many are found in DNA transposons?
sed 's/\./\t/2' HGTNewNames.list|awk '{print $1}' |cat <( bedtools intersect -wo -a ../../58_Renamatorium/24_IRF_DNATrans/SupportedIRFMergeClassified.gff -b ../../58_Renamatorium/1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq|wc )  - |sort|uniq -c |sort -k1,1nr|awk '$1==2' |wc
      0       0       0

```
