# How does the 738 genome compare to the previously existing SCN genome?

```
#/data021/GIF/remkv6/Baum/CamTechGenomeComparison/16_jgi_738_synteny/opscan

ln -s ../../15_transcriptsTo738/Braker/braker/final.repeat.scaffs.738/augustus.aa final.repeat.scaffs.738.aa
ln -s ../../15_transcriptsTo738/Braker/braker/final.repeat.scaffs.738/augustus.gff3 final.repeat.scaffs.738.gff3
ln -s ../../15_transcriptsTo738/Braker/braker/genome.738/augustus.aa genome.738.aa
ln -s ../../15_transcriptsTo738/Braker/braker/genome.738/augustus.gff3 genome.738.gff3
ln -s /data021/GIF//genome/sequences/Hetgly/JGI/Hetgly_JGI.gff3
ln -s /data021/GIF//genome/sequences/Hetgly/JGI/Hetgly_JGI.pep.fasta
cp ../../06_iadhore_allscaf/opscan/opscan.pbs

sed -e 's/ID=//g' genome.738.gff3 -e 's/;/\t/g' -e 's/NODE_//g'|grep gene -|awk '{print $3 "\t" $9 "\t" substr($1, 0, 8) "\t" $4 "\t" $5 "\t" $7 "\t" "t"}' |sed 's/scaffold_//g' |sort -k 5 -V|awk '{print ">"$2 "\t" "1" "\t" $4 "\t" $5 "\t" $6 "\t" $3}'|sort -k 1 -V |sed -e 's/F+//g' -e 's/F-//g' -e 's/F//g' >augustus.738.gff3.reformat1
tr "\n" "\t" < genome.738.aa |sed 's/>/\n>/g'|sed 's/\t//g'|grep ".t1"|sed  's/.t1/\t/g'>7384join2
join -1 1 -2 1 augustus.738.gff3.reformat1 7384join2 >joined.738.3
sed 's/ /\n/6' joined.738.3 >operon.db

#some modifications were needed for the gff3 from JGI.  the # typo was removed
awk 'FS="\t", $1=$1"JGI"{print}' Hetgly_JGI.gff3 >Hetgly_JGI.gff3.mod
sed -e 's/ID=//g' Hetgly_JGI.gff3.mod -e 's/;/\t/g' -e 's/NODE_//g'|grep gene -|awk '{print $3 "\t" $9 "\t"$1 "\t" $4 "\t" $5 "\t" $7 "\t" "t"}'|sed 's/scaffold_//g' |sort -k 5 -V|awk '{print ">"$2 "\t" "1" "\t" $4 "\t" $5 "\t" $6 "\t" $3}'|sort -k 1 -V  >JGI.gff3.reformat1
 tr "\n" "\t" <Hetgly_JGI.pep.fasta |sed 's/>/\n>/g'|sed 's/\t//g'|grep ".t1"|sed  's/.t1/\t/g'>JGI4join2
join -1 1 -2 1 JGI.gff3.reformat1 JGI4join2 >joined.JGI.3
#says they aren't sorted, but the correct number of genes is in the output.
sed 's/ /\n/6' joined.JGI.3 >genome.db

mkdir ../Nematode
cd ../Nematode

sed -e 's/ID=//g' ../opscan/genome.738.gff3 -e 's/;/\t/g' -e 's/NODE_//g'|grep gene -|awk '{print $3 "\t" $9 "\t" substr($1, 0, 8) "\t" $4 "\t" $5 "\t" $7 "\t" "t"}' |sed 's/scaffold_//g' |sort -k 5 -V|awk '{print ">"$2 "\t" "1" "\t"  4 "\t" $5 "\t" $6 "\t" $3}'|sort -k 1 -V |sed -e 's/F+//g' -e 's/F-//g' -e 's/F//g' >augustus.738.gff3.reformat1
tr "\n" "\t" < ../opscan/genome.738.aa |sed 's/>/\n>/g'|sed 's/\t//g'|grep ".t1"|sed  's/.t1/\t/g'>7384join2
join -1 1 -2 1 augustus.738.gff3.reformat1 7384join2 >joined.738.3
sed 's/ /\n/6' joined.738.3 >Nematodeoperon.db

### this step is no longer needed because of a slight alteration in the opscan.pbs script.
###cp OPSCAN.o111450.condo iADHoRe_PLAZA.table
###grep "CL"  ../opscan/OPSCAN.o111450.condo |awk '{print $3 "\t" $6}'  >iADHoRe_PLAZA.table
cp iADHoRe_PLAZA.table ../Nematode/
cd ../Nematode/
mkdir subject
mkdir query
cd query/
#FYI the next few steps can take a few minutes if there are lots of scaffolds to index
tr "\n" "\t" <../operon.db |sed 's/>/\n>/g'|awk '{print $1$5 " " $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini

cd ../subject/

tr "\n" "\t" < ../genome.db |sed 's/>/\n>/g'|awk '{print $1$5 " " $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini

cd ../
cat query/query.ini subject/subject.ini >tmp
tr "\t" " "<tmp >Nematode.ini

#this is what is put in the Nematode.ini file
genome=query
#sequence stuff
#need this empty newline
genome=subject
#sequence stuff
prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=true
blast_table= iADHoRe_PLAZA.table
output_path= output
alignment_method=gg2
gap_size=15
cluster_gap=20
level_2_only=false
q_value=.9


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=1:00:00
#PBS -N iadhore
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -q debug
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module load iAdHoRe/3.0.01
i-adhore Nematode.ini

##in case you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"
```

### Post processing

```
#number of contigs in JGI assembly with synteny to 738 genome
 cat <(awk '$2!=$4{print $3}' multiplicons.txt|grep JGI -) <(awk '$2!=$4{print $5}' multiplicons.txt|grep JGI -)|sort |uniq|sed -e 's/JGI//g' -e 's/^/scaffold_/g'|sort|uniq|wc
679 scaffolds

#number of scaffolds in the 738.genome.fa with synteny to JGI assembly
cat <(awk '$2!=$4 {print $3}' multiplicons.txt|grep -v JGI -) <(awk '$2!=$4 {print $5}' multiplicons.txt|grep -v JGI -)|sort |uniq|sed -e 's/JGI//g' -e 's/^/scaffold_/g'|sort|uniq|wc
573 scaffolds

#total length of scaffolds in JGIthat have synteny in 738, (this is not syntenic length)
join  -1 1 -2 1 <(sort Hetgly_JGI_length.txt) <(cat <(awk '$2!=$4{print $3}' multiplicons.txt|grep JGI -) <(awk '$2!=$4{print $5}' multiplicons.txt|grep JGI -)|sort |uniq|sed -e 's/JGI//g' -e 's/^/scaffold_/g'|sort)|awk '{print $2}'| awk -F, '{sum+=$1} END{print sum}'
122188557
#total length of scaffolds in 738 that have synteny in JGI, (this is not syntenic length)
join  -1 1 -2 1 <(sort genome.738.chr.len) <(cat <(awk '$2!=$4{print $3}' multiplicons.txt|grep -v JGI -) <(awk '$2!=$4{print $5}' multiplicons.txt|grep -v JGI -)|sort |uniq)|awk '{print $2}'|awk -F, '{sum+=$1} END{print sum}'
98881958
```

### circos
```
 #backtrack This is the initial synteny that used 1:1 rbbh instead of families

 cp -rf ../26_GloboderaSynteny/circos .
ln -s ../../27_pyscafGlobodera/pyScaf/SCNJGI/Hetgly_JGI.fasta
ln -s ../../27_pyscafGlobodera/pyScaf/SCNJGI/genome738.5.2.tsv
module load cdbfasta
module load bioawk
cdbfasta Hetgly_JGI.fasta
sed 's/ /</g' genome738.5.2.tsv|awk 'NF>2 {print $4}'|sed 's/</ /g'|awk 'NF>5' |sed 's/|quiver//g'|sed 's/</\n/g'|sed 's/>//g' |cdbyank genome738sl.polished.mitoFixed.noquiver.fa.cidx |bioawk -c fastx '{print "chr - "$name,$name,"0",length($seq),"blue"}' >738.mergeset.withoutjgi.kary
awk '{print $3}' 738.mergeset.withoutjgi.kary |grep -f - ../Nematode/output/multiplicons.txt |awk '{if ($2=="subject") print $3; else if($4=="subject")print $5}'|sort|uniq|sed 's/JGI//g'|sed 's/^/scaffold_/g'|cdbyank Hetgly_JGI.fasta.cidx|bioawk -c fastx '{print "chr - "$name,$name,"0",length($seq),"green"}' |cat - 738.mergeset.withoutjgi.kary >jgimerged.all.kary
ln -s ../opscan/Hetgly_JGI.gff3
ln -s ../../26_GloboderaSynteny/circos/738.polished.mitofixed.repmod.gff3
ln -s ../opscan/738.polished.mitofixed.repmod.gff3
ln -s ../../26_GloboderaSynteny/opscan/738.polished.mitofixed.repmod.gff3
grep gene Hetgly_JGI.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}' |sed 's/Note/\t/g'|awk '{print $1 "\t" $6"\t"$7}'|sed 's/\./\t/g' >tmp
cat <(grep gene 738.polished.mitofixed.repmod.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}') <(awk 'NF > 3 {print $1"\t"$3"\t"$4}' tmp) <(awk 'NF == 3 {print $0}' tmp ) >both.gene.coord
cat <(awk '$2=="query" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../Nematode/output/genes.txt|sort -k 1,1 -V ) <(awk '$2=="subject" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../Nematode/output/genes.txt|sort -k 1,1 -V) >reformat.genes.txt
 join -a 1  -1 1 -2 1 <(sort -k 1,1 -V reformat.genes.txt) <(sort -k 1,1 -V both.gene.coord) >gene.positions
 sed 's/ /_/2' gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >gene.positions.1
cat <(grep gene 738.polished.mitofixed.repmod.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}') <(awk 'NF > 3 {print $1"\t"$3"\t"$4}' tmp) <(awk 'NF == 3 {print $0}' tmp ) >both.gene.coord
grep gene Hetgly_JGI.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}' |sed 's/Note/\t/g'|awk '{print $1 "\t" $6"\t"$7}'|sed 's/\./\t/g' |awk 'NF >3 {print $1"\t"$3"\t"$4}' >TMP1
grep gene Hetgly_JGI.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}' |sed 's/Note/\t/g'|awk '{print $1 "\t" $6"\t"$7}'|sed 's/\./\t/g' |awk 'NF == 3 {print $0}' >TMP2
cat TMP1 TMP2 >JGI.gene.coord
#I think this is where I realized that I needed to process JGI and Camtech separately

awk '$2=="subject" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../Nematode/output/genes.txt|sort -k 1,1 -V >JGI.reformat.genes.txt
awk '$2=="query" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../Nematode/output/genes.txt|sort -k 1,1 -V >Camtech.reformat.genes.txt
join -a 1  -1 1 -2 1 <(sort -k 1,1 -V JGI.reformat.genes.txt) <(sort -k 1,1 -V JGI.gene.coord) >JGI.gene.positions
ln -s ../../19_elegans3SCN/opscan/genome738sl.polished.gff3
ln -s ../opscan/genome.738.gff3
grep gene genome.738.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}' >738.gene.coord
join -a 1  -1 1 -2 1 <(sort -k 1,1 -V Camtech.reformat.genes.txt) <(sort -k 1,1 -V 738.gene.coord) >Camtech.gene.positions
sed 's/ /_/2' Camtech.gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >Camtech.gene.positions.1
sed 's/ /_/2' JGI.gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >JGI.gene.positions.1
awk '{if ($2=="subject") print $1"\t"$4"\t"$5"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8"\t"$11"\t"$12"\t"$9"\t"$10"\t"$12; else print $0}'  ../Nematode/output/multiplicons.txt |awk '$4!=$2' >../Nematode/output/multiplicons.mod.txt
paste <(awk '{print $2"_"$3"_"$9}' ../Nematode/output/multiplicons.mod.txt |while read line; do grep -w $line Camtech.gene.positions.1 ;done |awk '{print $1,$2}' ) <(awk '{print
$2"_"$3"_"$10}' ../Nematode/output/multiplicons.mod.txt|while read line; do grep -w $line Camtech.gene.positions.1 ;done| awk '{print $1,$3}') <(awk '{print $4"_"$5"_"$11}' ../Nematode/output/multiplicons.mod.txt|while read line; do grep -w $line JGI.gene.positions.1 ;done | awk '{print $1,$2}' |uniq ) <(awk '{print $4"_"$5"_"$12}' ../Nematode/output/multiplicons.mod.txt |while read line; do grep -w $line JGI.gene.positions.1 ;done |awk '{print $1,$3}' |uniq)|sed 's/_/ /g'|awk '{print $2,$4,$8,$10,$12,$16}'>syntenic.ribbons.txt
awk '{print $3}' jgimerged.all.kary |grep -w -f - syntenic.ribbons.txt >syntenic.ribbons.subset.txt
sed 's/scaffold_//g' jgimerged.all.kary |awk '{if($7=="green") print $1,$2,$3"JGI",$4"JGI",$5,$6,$7; else print $0}' >improvedJGI.scn.kary

/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links syntenic.ribbons.subset.txt -karyotype improvedJGI.scn.kary -init_order 1055JGI,108JGI,1098JGI,109JGI,10JGI,112JGI,1154JGI,11JGI,129JGI,12JGI,133JGI,136JGI,138JGI,141JGI,145JGI,149JGI,151JGI,189JGI,1JGI,216JGI,230JGI,240JGI,24JGI,258JGI,266JGI,26JGI,270JGI,280JGI,28JGI,2JGI,320JGI,32JGI,336JGI,33JGI,352JGI,354JGI,381JGI,394JGI,417JGI,419JGI,423JGI,426JGI,43JGI,49JGI,50JGI,516JGI,522JGI,58JGI,610JGI,653JGI,68JGI,755JGI,85JGI,8JGI,9JGI -static_rx 1055JGI,108JGI,1098JGI,109JGI,10JGI,112JGI,1154JGI,11JGI,129JGI,12JGI,133JGI,136JGI,138JGI,141JGI,145JGI,149JGI,151JGI,189JGI,1JGI,216JGI,230JGI,240JGI,24JGI,258JGI,266JGI,26JGI,270JGI,280JGI,28JGI,2JGI,320JGI,32JGI,336JGI,33JGI,352JGI,354JGI,381JGI,394JGI,417JGI,419JGI,423JGI,426JGI,43JGI,49JGI,50JGI,516JGI,522JGI,58JGI,610JGI,653JGI,68JGI,755JGI,85JGI,8JGI,9JGI

sed -i 's/JGI//g' syntenic.ribbons.subset
awk '($3-$2)>30000' syntenic.ribbons.subset >syntenic.ribbons.subset.subset

cat <(awk '{print $1"\n"$4}' syntenic.ribbons.subset.subset|sed 's/tme/tmer/g' |cdbyank genome738sl.polished.mitoFixed.noquiver.fa.cidx| bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" green"}' ) <(awk '{print "scaffold_"$1"\nscaffold_"$4}' syntenic.ribbons.subset.subset|cdbyank Hetgly_JGI.fasta.cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" blue"}' ) |sort|uniq |sed 's/scaffold_//g' >kary.families
 /shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links syntenic.ribbons.subset.subset -karyotype kary.families -init_order 000001,000002,000011,000014,000015,000019,000021,000028,000029,000033,000037,000038,000041,000049,000058,000059,000066,000080K,000082,000089,000091,000099,000117,000118,000119,000120,000122,000123,000131,000132,000136,000138,000139,000167,000168,000171,000177,000188,000190,000195,000203,000212,000213,000217,000220,000221,000237,000251,000263,000286,000300,000308,000309,000324,000396,000404,000434K,000502,000523,000614K,000794K,000816K,000823K,000893K,001116K,001124,002293,002368,1syntmer,20syntmer,22syntmer,24syntmer,30syntmer,33syntmer,34syntmer,35syntmer,36syntmer,39syntmer,41syntmer,45syntmer,46syntmer,55syntmer -static_rx 000001,000002,000011,000014,000015,000019,000021,000028,000029,000033,000037,000038,000041,000049,000058,000059,000066,000080K,000082,000089,000091,000099,000117,000118,000119,000120,000122,000123,000131,000132,000136,000138,000139,000167,000168,000171,000177,000188,000190,000195,000203,000212,000213,000217,000220,000221,000237,000251,000263,000286,000300,000308,000309,000324,000396,000404,000434K,000502,000523,000614K,000794K,000816K,000823K,000893K,001116K,001124,002293,002368,1syntmer,20syntmer,22syntmer,24syntmer,30syntmer,33syntmer,34syntmer,35syntmer,36syntmer,39syntmer,41syntmer,45syntmer,46syntmer,55syntmer
awk '{print $1, $2,$5}' genome738sl.polished.mitoFixed.fa.out.gff|awk '{print $1 "\t" $3/10000}'  | sed 's/\./\t/g' | awk '{print $1 "\t" $2}'| sort |uniq -c|awk '{print $2 "\t" $3*10000 "\t" $3*10000+10000 "\t" $1}' |awk 'NR>3' >repeat.histogram
```

![circosPlot](assets/JGIsynteny.png)
