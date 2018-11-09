#  Trying to find out if tandem duplication has played a part in the evolution of virulence in H. glycines
### Using redtandem to identify tandem regions
```


#/data013/GIF/remkv6/Baum/CamTechGenomeComparison/40_tandemdups
#Setting up fasta file


bioawk -c fastx '{print ">SCN"$name" 1-"length($seq)" \n"$seq}' ../../18_mitochondria/swapMitoScaffolds/genome738sl.polished.mitoFixed.fa | sed 's/|quiver//g' >renamed.scn.fa
sed 's/>/>>scn/g' renamed.scn.fa |bioawk -c fastx '{print $name"_1-"length($seq)"\n"$seq}' |fold -79 >renamed.scn.reformatted.fa
sed -e 's/a/A/g' -e 's/c/C/g' -e 's/g/G/g' -e 's/t/T/g' -e 's/sCn/scn/g' renamed.scn.reformatted.fa >caps.renamed.scn.reformatted.fa
bioawk -c fastx '{print $name,length($seq)}' caps.renamed.scn.reformatted.fa |awk '$2>40000' |cdbyank caps.renamed.scn.reformatted.fa.cidx|fold -79 >subset.caps.renamed.scn.reformatted.fa

#redtandem pbs script that would not run


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=4:00:00
#PBS -N redtand738
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
module load redtandem
perl /shared/software/GIF/programs/redtandem/ReDtandem/ReDtandem.pl --species scn --noclean --dnafile ./subset.caps.renamed.scn.reformatted.fa

ssh condo "qstat -f ${PBS_JOBID} |head"

Tried to run redtandem on the whole genome, but some of the contigs were lagging and would not finish. SO I had to split the fasta


mkdir split
cd split
ln -s ../caps.renamed.scn.reformatted.fa .
~/common_scripts/fasta-splitter.pl --n-parts 736 caps.renamed.scn.reformatted.fa
for line in *.fa; do echo "perl /shared/software/GIF/programs/redtandem/ReDtandem/ReDtandem.pl --species scn --dnafile "$line;done >redtandem.sh

redtandem.pbs


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N redtand738
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
module load redtandem
module load parallel

parallel  --jobs 16 --sshloginfile $PBS_NODEFILE --joblog redtandem.log --workdir $PWD < redtandem.sh

ssh condo "qstat -f ${PBS_JOBID} |head"
Still had problems with 7 contigs which I ran in an interactive node and had to kill a tandem call once or twice for each.


awk '{print $14}' *.log |grep -w -v -f - <(ls *.fa) |sed 's/^/cat /g' |sed 's/$/ >>notfound.fa/g' >not.found.sh
grep ">" not.found.fa
I am pretty sure that these are giant tandem repeats.


000387 -- 1:>scn397_1-120455
#no real pattern here in jbrowse for this one.  not very many unique transcripts.
000731 -- 2010:>scn507_1-433642
# not very many unique transcripts. very few ccs reads spanning the length, last 40kb have super high ccs read depth
Total:  74,000
Count:  8
Mean:   9,250
Median: 7,500
Min:    6,000
Max:    25,000
000899 -- 9239:>scn527_1-904088
#This scaffold has two effectors:Esophageal and lcl|est_hggfha5A08.  with lots of repeats and multiple mapping transcripts.
Total:  157,000
Count:  15
Mean:   10,466
Median: 10,000
Min:    5,000
Max:    16,000

000221 -- 24309:>scn555_1-414772
#nothing really interesting herewith lots of repeats and multiple mapping transcripts..  
Total:  141,000
Count:  8
Mean:   17,625
Median: 16,000
Min:    7,000
Max:    32,000
000708 -- 31223:>scn559_1-471123
#same as above
000045 -- 39077:>scn676_1-375857
#two esophageal like proteins here, lower pacbio read depth.
Total:  205,000
Count:  20
Mean:   10,250
Median: 8,000
Min:    5,000
Max:    35,000

17syntmer -- 45343:>scn697_1-294611
#maybe 10 unique mapping transcripts here, low ccs read depth, lots mult transcripts
Total:  170,000
Count:  15
Mean:   11,333
Median: 8,000
Min:    6,000
Max:    30,000
Also looked at 000138, which is littered with transposons, and in the middle of them a relatively small Esophageal effector is encoded.


Total:  441,000
Count:  28
Mean:   15,750
Median: 15,000
Min:    5,000
Max:    32,000


qsub -I
cd $PBS_O_WORKDIR
module load redtandem
perl /shared/software/GIF/programs/redtandem/ReDtandem/ReDtandem.pl --species scn --dnafile notfound.fa

#concatenating all redtandem results
cat *.out* |sort -k6,6nr |awk '{print $0}' |sort -k1,1nr >final
#attempting to convert redtandem contig names to the original scn contig names.
while read line; do echo "grep -w "$line final;done <../renamingindex.reverse |awk '{print $1,$2,$3,$5,"|sed "$3,$4}' |sed "s|sed |sed 's/|g"|awk '{print $1,$2,$3,$4,$5,$6"/"$7"/g"}' |sed "s|/g|/g'|g" |awk '{
print $0" >> final.renamed.out" }' |sed 's/-w /-w "/g' |sed 's/ final /" final /g'|sed "s|sed 's/|awk '\$1==|g"| sed "s|/| {print |1" |sed 's/print /print "/g' |sed "s|/g'|- \$0}'|g"  |sed 's/-/",/2' >renamer.sh
#this applies renamer.sh to every redtandem.out* file.
for f in redtandem.out*;do echo "sed 's/final/"$f"/g' necc_Files/renamer.sh >$f.sh"; done >test
#this is a massive sh file that will run all of the renaming in parallel
for f in  redtandem.out*.sh;do echo "sh "$f; done >renameall

Now I need to make a comparison of these “tandemly duplicated regions”, because there are tons of them, and many of them are really large.


mkdir alignments
cd alignments
###moved all red tandem output to alignments

#Split tandem repeats into individual files.
for f in *.out; do split -l 1 $f $f ;done
for f in *.renamed.out*; do awk '{print $6-$5,$0}' $f |sort -k1,1Vr |awk '{print $2,$9}' |sed 's/,/\n/g' |awk 'NF' |awk 'BEGIN{i=0;}{if ($2==""){print i,$1;}else{print $0;i=$1}}' |sed 's/\.\./-/g' |sed 's/>//g' >$f.cdbyank; done
sed -i "s/^/cdbyank genome738sl.polished.mitoFixed.fixed.fa.cidx  -R -a '/g" *cdbyank

#get the cdbfasta database ready
sed 's/|quiver//g' ../../../18_mitochondria/swapMitoScaffolds/genome738sl.polished.mitoFixed.fa >genome738sl.polished.mitoFixed.fixed.fa
module load cdbfasta
sed -i "s/'//1" *cdbyank
for f in *cdbyank; do sh $f >$f.fa;done &

#tried to align with muscle, but only partially successful.
module load muscle
for f in *cdbyank.fa; do echo "muscle -in "$f" -out "$f".aln";done >muscle.sh

#Next attempt was to align these regions with blast, although the data is uninterpretable because of a naming issue earlier on.
for f in *cdbyank.fa; do echo "makeblastdb -in "$f" -dbtype nucl -out $f.blast.db; blastn -query $f -db $f.blast.db -out $f.blast.out -outfmt 6 ";done >blast.sh

This analysis can be finished and the tandem repeats verified if I blast the correct “ideal repeat unit” which is specified in column 4 and 5. With a proper blast, I can see how all of these tandem repeats are related and try to come up with a mechanism for how SCN evolves so quickly. From this analysis I can say that most of the tandem duplications identifed with redtandem are likely to be ancestral, although they are fewer. Most of the tandems are not the same size and are only duplications of parts of a consensus sequence.
what are the size of tandem repeats

because of the massive amounts of data that that redtandem generated,and since redtandem was focused at larger tandems repeats in the genome we focused on only those tandem repeats larger than 3kb.\ There are 67 tandem duplications that are larger than 10kb, and 330 that are larger than 3kb. The


#total count of tandem duplications in RedTandem
wc final
  2193  15351 348020 final
#count greater than 1kb
less final |awk '{ if(($5-$4)>1000) print $0}' |wc
    856    5992   81444
#count greater than 3kb
less final |awk '{ if(($5-$4)>3000) print $0}' |wc
    361    2527   36519
#count greater than 5kb
less final |awk '{ if(($5-$4)>5000) print $0}' |wc
    223    1561   23156



#this just adds the real scaffold name at column 1
#not sure how I made renamer.sh, but probably with a while loop.
<sxh>
sed -i 's/final/sorted.final/g' renamer.sh
sort -k1,1nr final >sorted.final
sh renamer.sh
wc sorted.final.renamed.out
  2193  17544 367107 sorted.final.renamed.out

Duplication size have a peak at pread length

Used RedTandem data only. #distribution of tandem repeats larger than 3kb. Hopefully this will reduce transposable elements, and focus on large scale variation. Is their a peak at our pread length? andrew says it was 12kb. what was the distribution of these after merging with cap3? did this have an effect?


less larger3kb.final.out |awk '{print $6-$5, $0}' |sort -k1,1nr|awk '{print $1/1000,$0}' |sed 's/\./\t/1' >sizeattached.final.renamed.txt
less larger3kb.final.out |awk '{print $6-$5, $0}' |sort -k1,1nr|awk '{print $1/1000,$0}' |sed 's/\./\t/1' |less
less larger3kb.final.out |awk '{print $6-$5, $0}' |sort -k1,1nr|awk '{print $1/1000,$0}' |sed 's/\./\t/1' |awk '{print $1}' |sort |uniq -c |sort -k1,1nr|less
doesnt appear to be a real peak.  it seems normally distributed.

Click to display ⇲

My thought is that by comparing the tandem duplications to the syntenic regions in the genome, I would figure out whether these tandem repeats are real or artifacts of the genome. If syntenic, then they are more likly to be artifacts of assembly. It seems strange, that many of these tandem duplicates have remarkably high identity (perc id blast), yet the gene models predictions are vastly different. Perhaps the poor mapping of rna-seq in these regions is the source of the the high frequency of “proteins of unknown function”. We've identified how to characterize this problem in heterozygous genome assemblies, and at the population level. If it is an artifact of assembly, then this could actually represent the diversity in the population, aka pan-genome.
Pangenome sizes can be quite large, so it makes sense that we are assembling this huge diversity from a population. Perhaps I should label these as duplicate artifact/pangenome tracks. If they are real duplications, the (variable gene calls and high perc iD???) then it may represent an individual in the population that was different, and simply had a large tandem duplication. Individual reads with tandem calls would be more likely to be included in the assembly.

#I wonder if this is where our duplicated buscos are hiding?
look pull out the gene names from the busco list and grep them out of merged gene gff. bedtools intersect
```
### subread overlap and annotation overlap counts
Individual TD investigation This 3kb repeat on scaffold 000693, has multiple tandems homologous to the first 25kb of scaffold 000745K
![Tandem duplicate with high read depth](assets/745k131bptandem.png)
```
submitted the longest ccs read to tandemrepeat finder online
The first 25kb of 000745K is comprised of a 131bp tandem repeat.

1--25462    131 194.5   131 98  0   49183   29  18  18  34  1.95
Looked at the largest red tandem duplicate located on

>53syntmer 734     94732   126820  95740   99765   3       94732..100772,107254..113294,120780..126820
This appears to cover 3 genic regions g28879, g28881,g28885.
The first two are highly identical, the last one is identical for the first 500bp and then appears to be divergent with interspersed blocks of homology in the the remaining ~4kb.
```
###  Comparing redtandem duplicates with subreads to see which overlap.
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/40_tandemdups/split\\
bedtools intersect -wo -a Repeats4overlap.gff -b ../../36_deduplicate122016/subreadMapping/subreads2genome.gff3 |awk '$4>($12*.9) &&(.9*$5)<$13 {print $1,$2,$3,$4,$5}' |sort |uniq|awk '{print $5-$4}' |summary.sh
Total:  17,378,922
Count:  17,568
Mean:   989
Median: 148
Min:    30
Max:    31,662

how many TD are there total?
wc Repeats4overlap.gff
  20009  180081 1285341 Repeats4overlap.gff


So TDs were at leaswt 90% verified by subread overlap for 17,568/20009 TDs. so 87.8% were verifiable, wihch includes the largest TD of 31kb.
<sxh>
What kind of annotations exist for the tandem duplicates?

bedtools intersect -wo -a ../40_tandemdups/split/Repeats4overlap.gff -b ../52_functional/augustusFunctionalAnnotation.gff3 |awk  '$12=="mRNA"' |cut -f 18 |grep -v "Not
e=NA" |grep "Note" |cut -d " " -f 2- |sed 's/;/\t/g' |cut -f 1 |sort |uniq -c |sort -k1,1nr |less
    214 Note=Similar to MUC19: Mucin-19 (Homo sapiens)
    186 Note=Similar to MATN1: Cartilage matrix protein (Gallus gallus)
    118 Note=Similar to MIMI_R196: Collagen-like protein 2 (Acanthamoeba polyphaga mimivirus)
    106 Note=Similar to fp-1: Foot protein 1 variant 1 (Perna viridis)
     89 Note=Similar to gpa-1: Guanine nucleotide-binding protein alpha-1 subunit (Caenorhabditis elegans)
     81 Note=Similar to BTBD2: BTB/POZ domain-containing protein 2 (Homo sapiens)
     79 Note=Similar to GBP: Glycophorin-binding protein (Plasmodium falciparum (isolate 3D7))
     78 Note=Similar to DDB_G0269086: Uncharacterized abhydrolase domain-containing protein DDB_G0269086 (Dictyostelium discoideum)
     66 Note=Similar to PF14_0175: Protein PF14_0175 (Plasmodium falciparum (isolate 3D7))
     55 Note=Similar to TbgDal_IV3690: Flagellar attachment zone protein 1 (Trypanosoma brucei gambiense (strain MHOM/CI/86/DAL972))
     54 Note=Similar to Syt4: Synaptotagmin-4 (Mus musculus)
     48 Note=Similar to Zan: Zonadhesin (Mus musculus)
     46 Note=Similar to Zfp26: Zinc finger protein 26 (Mus musculus)
     42 Note=Similar to C02C2.4: Uncharacterized transporter C02C2.4 (Caenorhabditis elegans)
     41 Note=Similar to K02A2.6: Uncharacterized protein K02A2.6 (Caenorhabditis elegans)
     40 Note=Similar to NEFH: Neurofilament heavy polypeptide (Homo sapiens)
     39 Note=Similar to nas-10: Zinc metalloproteinase nas-10 (Caenorhabditis elegans)
     39 Note=Similar to TCNA: Sialidase (Trypanosoma cruzi)
     34 Note=Similar to Slc16a14: Monocarboxylate transporter 14 (Mus musculus)
     33 Note=Similar to metS: Probable methionine--tRNA ligase%2C cytoplasmic (Dictyostelium discoideum)
     32 Note=Similar to C36A4.4: Probable UDP-N-acetylglucosamine pyrophosphorylase (Caenorhabditis elegans)
     32 Note=Similar to TRO: Trophinin (Homo sapiens)
     31 Note=Similar to MAGI3: Membrane-associated guanylate kinase%2C WW and PDZ domain-containing protein 3 (Homo sapiens)
     31 Note=Similar to sls: Titin (Drosophila melanogaster)
     29 Note=Similar to inx-3: Innexin-3 (Caenorhabditis elegans)
     27 Note=Similar to RPS6KA1: Ribosomal protein S6 kinase alpha-1 (Homo sapiens)
     27 Note=Similar to six1b: Homeobox protein six1b (Danio rerio)
     27 Note=Similar to ttn-1: Titin homolog (Caenorhabditis elegans)
     26 Note=Similar to Ppp1ca: Serine/threonine-protein phosphatase PP1-alpha catalytic subunit (Rattus norvegicus)
     25 Note=Similar to gcy-14: Receptor-type guanylate cyclase gcy-14 (Caenorhabditis elegans)
     25 Note=Similar to hmu: Halomucin (Haloquadratum walsbyi (strain DSM 16790 / HBSQ001))
     25 Note=Similar to PCNT: Pericentrin (Homo sapiens)
     25 Note=Similar to pif1: ATP-dependent DNA helicase PIF1 (Danio rerio)
     24 Note=Similar to clec-180: C-type lectin domain-containing protein 180 (Caenorhabditis elegans)
     24 Note=Similar to hlh-8: Twist-related protein (Caenorhabditis elegans)
     23 Note=Similar to F54H12.2: Uncharacterized protein F54H12.2 (Caenorhabditis elegans)
     22 Note=Similar to infB: Translation initiation factor IF-2 (Clostridium phytofermentans (strain ATCC 700394 / DSM 18823 / ISDg))
     22 Note=Similar to RRBP1: Ribosome-binding protein 1 (Canis lupus familiaris)
     21 Note=Similar to Pol-RFamide neuropeptides (Polyorchis penicillatus)
     20 Note=Similar to ANKS1A: Ankyrin repeat and SAM domain-containing protein 1A (Homo sapiens)
```

### synteny overlap
```
#How many syntenic regions do the duplicates overlap?
 bedtools intersect -wo -a Repeats4overlap.gff -b ../../46_SCN_syntenic_tracks/all.synteny.gff |awk '{print $18}' |sort|uniq|wc
    586     586   13371


#How many regions do these syntenic regions map to?
bedtools intersect -wo -a Repeats4overlap.gff -b ../../46_SCN_syntenic_tracks/all.synteny.gff |awk '{print $10,$13,$14}' |sort|uniq|wc
    558    1674   11536
```
### Effector overlap
```
#how many known effectors total in TD
cat TotalSectionsGene.list ../18_effectorRedo/121EffectorOnlyNewGeneNames.list |sort|uniq -c |awk '$1==2' |wc
     38      76     988
#how many effectors just in reiterated?     
[remkv6@condo031 27_TandemRedo]$ cat ReiteratedSectionsGene.list ../18_effectorRedo/121EffectorOnlyNewGeneNames.list |sort|uniq -c |awk '$1==2' |wc
     30      60     780
How many effectors in consensus?     
[remkv6@condo031 27_TandemRedo]$ cat ConsensusSectionsGene.list ../18_effectorRedo/121EffectorOnlyNewGeneNames.list |sort|uniq -c |awk '$1==2' |wc
     12      24     312
#remember, this is because of redtandems weird calling of consensus, so they do not always add up.


#how many of the predicted effector total?
cat ConsensusSectionsGene.list ../26_ExpressionSets/AllPredictedEffectors.list |sort|uniq -c |awk '$1==2' |wc
     36      72     936
[remkv6@condo031 27_TandemRedo]$ cat TotalSectionsGene.list ../26_ExpressionSets/AllPredictedEffectors.list |sort|uniq -c |awk '$1==2' |wc
    108     216    2808
[remkv6@condo031 27_TandemRedo]$ cat ReiteratedSectionsGene.list ../26_ExpressionSets/AllPredictedEffectors.list |sort|uniq -c |awk '$1==2' |wc
     91     182    2366

so out of 350 total predicted effectors, 30.8% of effector genes are contained within 14.7% of the genome.
```

### Genome without duplicates
```
#kept the consensus repeat unmasked in redtandem. #/work/GIF/remkv6/Baum/CamTechGenomeComparison/40_tandemdups/split/1tandemcopyonly.masked.genome738

Going to have to make 2 gff files.  One with all of the repeats, and one with the consensus repeat. Use bedtools to find in the repeats file which coordinates DO NOT OVERLAP the consensus.  Those can plug in below.
<sxh>
less sorted.final.renamed.out |awk '{print $1,$8}'| tr "," "\n" | awk 'BEGIN{i=0;}{if ( NF == 1 ){print i,$1;}else{print $0;i=$1;}}' |sed '/^$/d' |sed 's/\.\./\t/g'| awk '{print $1,"redtandem","gene",$2,$3,".","+",".","ID="$1$2$3";"}' |sed 's/>//g'|tr " " "\t" >Repeats4overlap.gff
#how many repeats in total
wc Repeats4overlap.gff
 20009  180081 1285341 Repeats4overlap.gff

awk 'print $1,"redtandem","gene",$5,$6,".","+",".","ID="$1$5$6"consensus;"}' sorted.final.renamed.out|sed 's/>//g'|tr " " "\t" >Repeats.consensus.4overlap.gff
awk '{if($4=="0") print $1,$2,$3,"1",$5,$6,$7,$8,$9; else print $0}' Repeats.consensus.4overlap.gff |tr " " "\t" >Repeats.consensus.4overlap1.gff
#how many repeats is this?
wc Repeats.consensus.4overlap1.gff
 2193  19737 159806 Repeats.consensus.4overlap1.gff

bedtools subtract -b Repeats.consensus.4overlap1.gff -a Repeats4overlap.gff >CoordinatesToMask.gff
#how many in the reiterated portion?
wc CoordinatesToMask.gff
20,654 156582 1118578 CoordinatesToMask.gff


###DO THE GENOME SIZES MATCH UP?
[remkv6@condo003 split]$ awk '{print $5-$4}' Repeats4overlap.gff |summary.sh
Total:  18,265,983
Count:  20,009
Mean:   912
Median: 134
Min:    30
Max:    31,662
[remkv6@condo003 split]$ awk '{print $5-$4}' Repeats.consensus.4overlap1.gff |summary.sh
Total:  3,962,596
Count:  2,193
Mean:   1,806
Median: 716
Min:    31
Max:    31,873
[remkv6@condo003 split]$ awk '{print $5-$4}' CoordinatesToMask.gff |summary.sh
Total:  14,561,727
Count:  20,711
Mean:   703
Median: 131
Min:    0
Max:    31,662

No, but they should not since the a reiteration frequently overlaps with a consensus, so if they overlapped in the gff, I removed them.  This is an improvement from the above, but still may not be perfect.

#####This breaks out the above consensus, removing 200kb from the total, and splitting some of the iterations into two because the consensus was in the middle of one reiteration, thereby splitting it into two.  This is just the way the red tandem software works.


#This masked those coordinates in the genome.
bedtools maskfasta -fi ../../18_mitochondria/swapMitoScaffolds/genome738sl.polished.mitoFixed.fa -bed CoordinatesToMask.gff  -fo 1tandemcopyonly.masked.genome738

#How much of the genome is masked?
~/common_scripts/new_Assemblathon.pl 1tandemcopyonly.masked.genome738
11.64% == 123860871 * 0.1164 = 14,417,405


How many genes does that take away?
bedtools intersect -wo -a CoordinatesToMask.gff -b ../../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3|awk '$12=="exon"' |cut -f 18 |sed 's/exon/\t/g' |cut -f 1 |sort|uniq |wc
  5917    5917   80618

#how many genes are left after tandem repeat expulsion
bedtools subtract -N -a ../../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 -b CoordinatesToMask.gff |awk '$3=="gene"' |wc
 23918  215262 1284963


#final gff after tandem removal
bedtools subtract -N -a ../../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 -b CoordinatesToMask.gff >1tandemcopyonly.masked.genome738.gff3

#if 85% of a feature
bedtools subtract -A -f .85 -a 1tandemcopyonly.masked.genome738.gff3 -b ../../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '$3=="gene"' |wc
 21882  196938 1175634

If 30% of a the feature is overlapped with the repeats
bedtools subtract -A -f .4 -a 1tandemcopyonly.masked.genome738.gff3 -b ../../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff |awk '$3=="gene"' |wc
 19974  179766 1072848

What is the distribution of tandem repeats in the 1 tandem masked genome?
less CoordinatesToMask.gff |awk '{print $5-$4}' |awk '{if($1<1000){print 1000}else print $1/1000}' |sed 's/\./\t/g'|awk '{if(NF>1) {print $1*1000} else {print $1}}' |sort |uniq -c |sort -k2,2 |less
```

### BUSCO overlap
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/47_busco_nematode/Busco_Prot/run_H.glycines.pep.fasta.out

less full_table_H.glycines.pep.fasta.out.tsv |awk '$2=="Duplicated"' |sed 's/\./\t/1' |grep -w "t1"  - |awk '{print $1}' |awk '{print $1}' |sort|uniq -c |awk '$1>1{print $2}' |grep -w -f - <(less full_table_H.glycines.pep.fasta.out.tsv |awk '$2=="Duplicated"') |awk '{print $3}' |sed 's/\./\t/1' |cut -f 1 |sort|uniq >DuplicatedBuscosGene.list


#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/26_ExpressionSets
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/47_busco_nematode/Busco_Prot/run_H.glycines.pep.fasta.out/DuplicatedBuscosGene.list .
cp ../1_genomeNgff/geneRenamer.sh .
sed -i 's/augustus.aa/DuplicatedBuscosGene.list/g' geneRenamer.sh

#How many duplicated buscos overlap with redtandem duplicates     
cat DuplicatedBuscosGene.list ../27_TandemRedo/TotalSectionsGene.list |sort|uniq -c |awk '$1==2' |wc
    79     158    2054

#what about just the reiterated portion?
cat DuplicatedBuscosGene.list ../27_TandemRedo/ReiteratedSectionsGene.list |sort|uniq -c |awk '$1==2' |wc
    71     142    1846
```

Ortholog overlap
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/40_tandemdups/split
#How many orthologues are found in the ~18MB of tandem duplications
bedtools intersect -wo -a Repeats4overlap.gff -b ../../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/ID=//g' |sed 's/;//g' |cat - ../../62_totalOrthologues/SCN.all.orthologues.list |sort|uniq -c |awk '$1==2' |wc
   2498    4996   36167

##############################################   
The below is outdated, and probably unnecessary to redo.  leaving as is.
#How many are in the 12.5MB reiterated portion of the TD?
bedtools intersect -wo -a CoordinatesToMask_sorted.gff -b ../../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/ID=//g' |sed 's/;//g' |cat - ../../62_totalOrthologues/SCN.all.orthologues.list |sort|uniq -c |awk '$1==2' |wc
   1767    3534   25547
To which species do these orthologues belong?
#b xylophilus
bedtools intersect -wo -a CoordinatesToMask_sorted.gff -b ../../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/ID=//g' |sed 's/;//g' |cat - ../../62_totalOrthologues/SCN.all.orthologues.list |sort|uniq -c |awk '$1==2 {print $2}' |cat - ../../62_totalOrthologues/Bxyl.ortho.list |sort|uniq -c| awk '$1==2'|wc
    559    1118    8069
#g ellingtonae
 bedtools intersect -wo -a CoordinatesToMask_sorted.gff -b ../../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/ID=//g' |sed 's/;//g' |cat - ../../62_totalOrthologues/SCN.all.orthologues.list |sort|uniq -c |awk '$1==2 {print $2}' |cat - ../../62_totalOrthologues/Gell.ortho.list |sort|uniq -c| awk '$1==2'|wc
   1150    2300   16593
#G. pallida
 bedtools intersect -wo -a CoordinatesToMask_sorted.gff -b ../../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/ID=//g' |sed 's/;//g' |cat - ../../62_totalOrthologues/SCN.all.orthologues.list |sort|uniq -c |awk '$1==2 {print $2}' |cat - ../../62_totalOrthologues/Gpal.ortho.list |sort|uniq -c| awk '$1==2'|wc
    848    1696   12246
#G. rostochiensis
bedtools intersect -wo -a CoordinatesToMask_sorted.gff -b ../../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/ID=//g' |sed 's/;//g' |cat - ../../62_totalOrthologues/SCN.all.orthologues.list |sort|uniq -c |awk '$1==2 {print $2}' |cat - ../../62_totalOrthologues/Gros.ortho.list |sort|uniq -c| awk '$1==2'|wc
   1082    2164   15629
#M. hapla
bedtools intersect -wo -a CoordinatesToMask_sorted.gff -b ../../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/ID=//g' |sed 's/;//g' |cat - ../../62_totalOrthologues/SCN.all.orthologues.list |sort|uniq -c |awk '$1==2 {print $2}' |cat - ../../62_totalOrthologues/Mhap.ortho.list |sort|uniq -c| awk '$1==2'|wc
    592    1184    8558   
#M. incognita
bedtools intersect -wo -a CoordinatesToMask_sorted.gff -b ../../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/ID=//g' |sed 's/;//g' |cat - ../../62_totalOrthologues/SCN.all.orthologues.list |sort|uniq -c |awk '$1==2 {print $2}' |cat - ../../62_totalOrthologues/Minc.ortho.list |sort|uniq -c| awk '$1==2'|wc
    535    1070    7733

###############################################################################  
```
### New Updated analysis of Tandem duplicates
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/27_TandemRedo

cd 27_TandemRedo/
cp ../../40_tandemdups/split/Repeats.consensus.4overlap1.gff  ConsensusSections.gff
cp ../../40_tandemdups/split/CoordinatesToMask.gff ReiteratedSections.gff
cp ../../40_tandemdups/split/Repeats4overlap.gff  TotalSections.gff
sed 's/>//g' scaffold.renamer.sh | sed 's/genome738sl.polished.mitoFixed.fa/ConsensusSections.gff/g' >scaffoldrenamer1
sed 's/>//g' scaffold.renamer.sh | sed 's/genome738sl.polished.mitoFixed.fa/ReiteratedSections.gff/g' >scaffoldrenamer2
sed 's/>//g' scaffold.renamer.sh | sed 's/genome738sl.polished.mitoFixed.fa/TotalSections.gff/g' >scaffoldrenamer3
sh scaffoldrenamer1 &
sh scaffoldrenamer2 &
sh scaffoldrenamer3 &
#how many genes overlap these TD
bedtools intersect -wo -a ReiteratedSections.gff -b ../1_genomeNgff/fixed.augustus.gff3|awk '$12=="exon"' |cut -f 18 |sed 's/exon/\t/g' |cut -f 1 |sort|uniq |wc                                               5917    5917  142008
bedtools intersect -wo -a ConsensusSections.gff -b ../1_genomeNgff/fixed.augustus.gff3|awk '$12=="exon"' |cut -f 18 |sed 's/exon/\t/g' |cut -f 1 |sort|uniq |wc
   1947    1947   46728
bedtools intersect -wo -a TotalSections.gff -b ../1_genomeNgff/fixed.augustus.gff3|awk '$12=="exon"' |cut -f 18 |sed 's/exon/\t/g' |cut -f 1 |sort|uniq |wc
   6767    6767  162408
Getting Gene lists for above for overlap comparisons
bedtools intersect -wo -a ReiteratedSections.gff -b ../1_genomeNgff/fixed.augustus.gff3|awk '$12=="exon"' |cut -f 18 |sed 's/exon/\t/g' |cut -f 1 |sort|uniq |sed 's/\./\t/2' |sed 's/ID=//g' |cut -f 1 |sed 's/T/G/g'|sort|uniq>ReiteratedSectionsGene.list
[remkv6@condo031 27_TandemRedo]$ bedtools intersect -wo -a ConsensusSections.gff -b ../1_genomeNgff/fixed.augustus.gff3|awk '$12=="exon"' |cut -f 18 |sed 's/exon/\t/g' |cut -f 1 |sort|uniq |sed 's/\./\t/2' |sed 's/ID=//g' |cut -f 1 |sed 's/T/G/g'|sort|uniq>ConsensusSectionsGene.list
bedtools intersect -wo -a TotalSections.gff -b ../1_genomeNgff/fixed.augustus.gff3|awk '$12=="exon"' |cut -f 18 |sed 's/exon/\t/g' |cut -f 1 |sort|uniq |sed 's/\./\t/2' |sed 's/ID=//g' |cut -f 1 |sed 's/T/G/g'|sort|uniq>TotalSectionsGene.list   

ONTOLOGIZER
bedtools intersect -wo -a ../TotalSections.gff -b ../../1_genomeNgff/fixed.augustus.gff3|awk '$12=="exon"' |cut -f 18 |sed 's/exon/\t/g' |cut -f 1 |sort|uniq |sed 's/\./\t/2' |sed 's/ID=//g' |cut -f 1 |sed 's/T/G/g'|sort|uniq>TotalSectionsGene.list
ln -s ../../10_tandemDups/ontologenizer/go.obo
ln -s ../../10_tandemDups/ontologenizer/Ontologizer.jar
ln -s ../../10_tandemDups/ontologenizer/simpleformat.ids
less ../../1_genomeNgff/fixed.augustus.gff3|awk '$3=="exon"' |cut -f 9 |sed 's/exon/\t/g' |cut -f 1 |sort|uniq |sed 's/\./\t/2' |sed 's/ID=//g' |cut -f 1 |sed 's/T/G/g'|sort|uniq>PopulationGene.list


module load java
java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p PopulationGene.list -s TotalSectionsGene.list
################################
no GO terms were enriched here.
```
