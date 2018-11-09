#  Need to check synteny with other nematodes to verify our assembly quality, and identify those more highly conserved genomic regions

### C. elegans
```
###/data021/GIF/remkv6/Baum/CamTechGenomeComparison/19_elegans3SCN/opscan What happens to the missing classes of genes in parasitic nematodes , like the classes discussed Ditylenchus destructor genome paper. What synteny exists between the species?
This is the polished 738 genome that still has the large 000114 mitochondrial scaffold, not the modified version

Download the c. elegans gff3, .aa file, and genome fasta


wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/genomic c_elegans.canonical_bioproject.current_development.genomic.fa.gz
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS255/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS255.annotations.gff3.gz
wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS255/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS255.protein.fa.gz

#This gets input ready for opscan

sed -e 's/ID=//g' genome738sl.polished.gff3 -e 's/;/\t/g' -e 's/NODE_//g'|grep gene -|awk '{print $3 "\t" $9 "\t" substr($1, 0,8) "\t" $4 "\t" $5 "\t" $7 "\t" "t"}' |sed 's/scaffold_//g' |sort -k 5 -V|awk '{print ">"$2 "\t" "1" "\t" $4 "\t" $5 "\t" $6 "\t" $3}'|sort -k 1 -V |sed -e 's/F+//g' -e 's/F-//g' -e 's/F//g' >augustus.738.gff3.reformat1
[2016-10-06 08:26:18] tr "\n" "\t" < genome738sl.polished.aa |sed 's/>/\n>/g'|sed 's/\t//g'|grep ".t1"|sed  's/.t1/\t/g'>7384join2
[2016-10-06 08:26:27] join -1 1 -2 1 augustus.738.gff3.reformat1 7384join2 >joined.738.3
[2016-'10-06 08:26:35] sed 's/ /\n/6' joined.738.3 >operon.db


#I've have not done the c. elegans genome conversion for opscan previously, so there were lots of things that changed. Wow C. elegans gff was a huge pain.

sed -e 's/ID=//g' c_elegans.PRJNA13758.WS255.annotations.gff3 -e 's/;/\t/g' -e 's/NODE_//g'|grep gene -|awk '{print $3 "\t" $9 "\t" substr($1, 0, 8) "\t" $4 "\t" $5 "\t" $7 "\t"
"t"}' |sed 's/scaffold_//g' |sort -k 5 -V|awk '{print ">"$2 "\t" "1" "\t" $4 "\t" $5 "\t" $6 "\t" $3}'|sort -k 1 -V |sed -e 's/F+//g' -e 's/F-//g' -e 's/F//g'|sed 's/Gene://g' >c.elegans.join1
awk '{print $1" "$3}' c_elegans.PRJNA13758.WS255.protein.fa|tr "\n" "\t"| sed 's/>/\n>/g'|cut -d " " -f 2-|sed 's/WBGene/>WBGene/g'|sort -k1,1V |awk '!seen[$1]++'  >c.elegans.join2
16616:16616  [2016-10-06 11:08:35] join -1 1 -2 1 <(sort c.elegans.join1) <(sort c.elegans.join2) >c.elegans.join32
join  -1 1 -2 1 <(sort c.elegans.join1) <(sort c.elegans.join2 ) >c.elegans.join3
sed 's/ /\n/6' c.elegans.join3|sed '/^>/! s/ //g' >genome.db


mkdir ../Nematode
cp operon.db ../Nematode .
cp genome.db ../Nematode .
cp ../opscan/iADHoRe_PLAZA.table.families .
cd ../Nematode
pwd
/data021/GIF/remkv6/Baum/CamTechGenomeComparison/19_elegans3SCN/Nematode
mkdir query
mkdir subject
cp ../../16_jgi_738_synteny/Nematode/iadhore.pbs .
cd query
#FYI the next few steps can take a few minutes if there are lots of scaffolds to index

```

Created a new script to find the optimal gap and gapcluster settings for iadhore.
```
i-adhore Nematode.ini
sed -e 's/gap_size=10/gap_size=20/g'  Nematode.ini -e 's/cluster_gap=20/cluster_gap=30/g'    -e 's/output_path= output/output_path= 30output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                        
sed -e 's/gap_size=20/gap_size=30/g'  Nematode.ini -e 's/cluster_gap=30/cluster_gap=40/g'    -e 's/output_path= 30output/output_path= 40output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=30/gap_size=40/g'  Nematode.ini -e 's/cluster_gap=40/cluster_gap=50/g'    -e 's/output_path= 40output/output_path= 50output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=40/gap_size=50/g'  Nematode.ini -e 's/cluster_gap=50/cluster_gap=60/g'    -e 's/output_path= 50output/output_path= 60output/g' > tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=50/gap_size=60/g'  Nematode.ini -e 's/cluster_gap=60/cluster_gap=70/g'    -e 's/output_path= 60output/output_path= 70output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=60/gap_size=70/g'  Nematode.ini -e 's/cluster_gap=70/cluster_gap=80/g'    -e 's/output_path= 70output/output_path= 80output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=70/gap_size=80/g'  Nematode.ini -e 's/cluster_gap=80/cluster_gap=90/g'    -e 's/output_path= 80output/output_path= 90output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=80/gap_size=90/g'  Nematode.ini -e 's/cluster_gap=90/cluster_gap=100/g'   -e 's/output_path= 90output/output_path= 100output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=90/gap_size=100/g' Nematode.ini -e 's/cluster_gap=100/cluster_gap=110/g'  -e 's/output_path= 100output/output_path= 110output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=100/gap_size=110/g' Nematode.ini -e 's/cluster_gap=110/cluster_gap=120/g' -e 's/output_path= 110output/output_path= 120output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=110/gap_size=120/g' Nematode.ini -e 's/cluster_gap=120/cluster_gap=130/g' -e 's/output_path= 120output/output_path= 130output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=120/gap_size=130/g' Nematode.ini -e 's/cluster_gap=130/cluster_gap=140/g' -e 's/output_path= 130output/output_path= 140output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=130/gap_size=140/g' Nematode.ini -e 's/cluster_gap=140/cluster_gap=150/g' -e 's/output_path= 140output/output_path= 150output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=140/gap_size=150/g' Nematode.ini -e 's/cluster_gap=150/cluster_gap=160/g' -e 's/output_path= 150output/output_path= 160output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=150/gap_size=160/g' Nematode.ini -e 's/cluster_gap=160/cluster_gap=170/g' -e 's/output_path= 160output/output_path= 170output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=160/gap_size=170/g' Nematode.ini -e 's/cluster_gap=170/cluster_gap=180/g' -e 's/output_path= 170output/output_path= 180output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=170/gap_size=180/g' Nematode.ini -e 's/cluster_gap=180/cluster_gap=190/g' -e 's/output_path= 180output/output_path= 190output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=180/gap_size=190/g' Nematode.ini -e 's/cluster_gap=190/cluster_gap=200/g' -e 's/output_path= 190output/output_path= 200output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=190/gap_size=200/g' Nematode.ini -e 's/cluster_gap=200/cluster_gap=210/g' -e 's/output_path= 200output/output_path= 210output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=200/gap_size=210/g' Nematode.ini -e 's/cluster_gap=210/cluster_gap=220/g' -e 's/output_path= 210output/output_path= 220output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=210/gap_size=220/g' Nematode.ini -e 's/cluster_gap=220/cluster_gap=230/g' -e 's/output_path= 220output/output_path= 230output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=220/gap_size=230/g' Nematode.ini -e 's/cluster_gap=230/cluster_gap=240/g' -e 's/output_path= 230output/output_path= 240output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=230/gap_size=240/g' Nematode.ini -e 's/cluster_gap=240/cluster_gap=250/g' -e 's/output_path= 240output/output_path= 250output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=240/gap_size=250/g' Nematode.ini -e 's/cluster_gap=250/cluster_gap=260/g' -e 's/output_path= 250output/output_path= 260output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=250/gap_size=260/g' Nematode.ini -e 's/cluster_gap=260/cluster_gap=270/g' -e 's/output_path= 260output/output_path= 270output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=260/gap_size=270/g' Nematode.ini -e 's/cluster_gap=270/cluster_gap=280/g' -e 's/output_path= 270output/output_path= 280output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=270/gap_size=280/g' Nematode.ini -e 's/cluster_gap=280/cluster_gap=290/g' -e 's/output_path= 280output/output_path= 290output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=280/gap_size=290/g' Nematode.ini -e 's/cluster_gap=290/cluster_gap=300/g' -e 's/output_path= 290output/output_path= 300output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=290/gap_size=300/g' Nematode.ini -e 's/cluster_gap=300/cluster_gap=310/g' -e 's/output_path= 300output/output_path= 310output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=300/gap_size=310/g' Nematode.ini -e 's/cluster_gap=310/cluster_gap=320/g' -e 's/output_path= 310output/output_path= 320output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=310/gap_size=320/g' Nematode.ini -e 's/cluster_gap=320/cluster_gap=330/g' -e 's/output_path= 320output/output_path= 330output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=320/gap_size=330/g' Nematode.ini -e 's/cluster_gap=330/cluster_gap=340/g' -e 's/output_path= 330output/output_path= 340output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=330/gap_size=340/g' Nematode.ini -e 's/cluster_gap=340/cluster_gap=350/g' -e 's/output_path= 340output/output_path= 350output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=340/gap_size=350/g' Nematode.ini -e 's/cluster_gap=350/cluster_gap=360/g' -e 's/output_path= 350output/output_path= 360output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=350/gap_size=360/g' Nematode.ini -e 's/cluster_gap=360/cluster_gap=370/g' -e 's/output_path= 360output/output_path= 370output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=360/gap_size=370/g' Nematode.ini -e 's/cluster_gap=370/cluster_gap=380/g' -e 's/output_path= 370output/output_path= 380output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=370/gap_size=380/g' Nematode.ini -e 's/cluster_gap=380/cluster_gap=390/g' -e 's/output_path= 380output/output_path= 390output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=380/gap_size=390/g' Nematode.ini -e 's/cluster_gap=390/cluster_gap=400/g' -e 's/output_path= 390output/output_path= 400output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=390/gap_size=400/g' Nematode.ini -e 's/cluster_gap=400/cluster_gap=410/g' -e 's/output_path= 400output/output_path= 410output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=400/gap_size=410/g' Nematode.ini -e 's/cluster_gap=410/cluster_gap=420/g' -e 's/output_path= 410output/output_path= 420output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=410/gap_size=420/g' Nematode.ini -e 's/cluster_gap=420/cluster_gap=430/g' -e 's/output_path= 420output/output_path= 430output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=420/gap_size=430/g' Nematode.ini -e 's/cluster_gap=430/cluster_gap=440/g' -e 's/output_path= 430output/output_path= 440output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=430/gap_size=440/g' Nematode.ini -e 's/cluster_gap=440/cluster_gap=450/g' -e 's/output_path= 440output/output_path= 450output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=440/gap_size=450/g' Nematode.ini -e 's/cluster_gap=450/cluster_gap=460/g' -e 's/output_path= 450output/output_path= 460output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=450/gap_size=460/g' Nematode.ini -e 's/cluster_gap=460/cluster_gap=470/g' -e 's/output_path= 460output/output_path= 470output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=460/gap_size=470/g' Nematode.ini -e 's/cluster_gap=470/cluster_gap=480/g' -e 's/output_path= 470output/output_path= 480output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=470/gap_size=480/g' Nematode.ini -e 's/cluster_gap=480/cluster_gap=490/g' -e 's/output_path= 480output/output_path= 490output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini                                                                                           
sed -e 's/gap_size=480/gap_size=490/g' Nematode.ini -e 's/cluster_gap=490/cluster_gap=500/g' -e 's/output_path= 490output/output_path= 500output/g'> tmp; mv tmp Nematode.ini
i-adhore Nematode.ini            
cut -f 2- all.multiplicons.txt |sort -k2,2V |uniq|awk '{print NR"\t"$0}' >all.uniq.multiplicons   
```
Circos
```

cp ../../26_GloboderaSynteny/circos/bands.conf .
cp ../../26_GloboderaSynteny/circos/circos.test.delete .
cp ../../26_GloboderaSynteny/circos/karyotype.both.again.subset.roundup .
cp ../../26_GloboderaSynteny/circos/ticks.conf .
cp ../../26_GloboderaSynteny/circos/ideogram.conf .

	cat <(grep gene ../opscan/genome738sl.polished.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}') <(awk '$3=="gene"' ../opscan/c_elegans.PRJNA13758.WS255.annotations.gff3 |sed 's/;/\t/g'|awk '{print $9"\t"$4"\t"$5}' |sed 's/ID=//g' |grep "Gene" -|sed 's/Gene://g') >both.gene.coord
33807  [2016-11-17 12:58:55] cat <(awk '$2=="H.glycines" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../Nematode/output/genes.txt|sort -k 1,1 -V ) <(awk '$2=="C.elegans" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../Nematode/output/genes.txt|sort -k 1,1 -V) >reformat.genes.txt
join -a 1  -1 1 -2 1 <(sort -k 1,1 -V reformat.genes.txt) <(sort -k 1,1 -V both.gene.coord) >gene.positions
sed 's/ /_/2' gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >gene.positions.1
paste <(awk '{print $2"_"$3"_"$9}' ../Nematode/all.uniq.multiplicons |while read line; do grep -w $line gene.positions.1 ;done |awk '{print $1,$2}' ) <(awk '{print $2"_"$3"_"$10}
' ../Nematode/all.uniq.multiplicons|while read line; do grep -w $line gene.positions.1 ;done| awk '{print $1,$3}') <(awk '{print $4"_"$5"_"$11}' ../Nematode/all.uniq.multiplicons|while read line; do grep -w
$line gene.positions.1 ;done | awk '{print $1,$2}' ) <(awk '{print $4"_"$5"_"$12}' ../Nematode/all.uniq.multiplicons |while read line; do grep -w $line gene.positions.1 ;done |awk '{print $1,$3}' )|sed 's/_/ /g'|awk '{print $2,$4,$8,$10,$12,$16}'|sort|uniq>syntenic.ribbons.txt

karyotype file update for relevant scaffolds


awk '{print $1}' syntenic.ribbons.txt |uniq|cdbyank genome738sl.polished.mitoFixed.noquiver.fa.cidx |bioawk -c fastx '{print $name, length($seq)}'
bioawk -c fastx '{print $name,length($seq)}' c_elegans.canonical_bioproject.current_development.genomic.fa
Final output found a few syntenic regions, but some are probably artifactual due to the large gap size allowed.The ribbons were huge here, so Im only showing the threads, which do not represent the actual size of each syntenic region

modified circos plot to find the highest confidence multiplicons that are unique


less all.uniq.multiplicons |awk '{print $9,$10,$11,$12}' |sort|uniq|awk '{print $0,$2-$1,$4-$3}' |sort -k1,1n >tmp
awk '{print $1"\t"$2"\t"$3"\t"$4}' tmp |grep -f - all.uniq.multiplicons >multiplicon.retry
paste <(awk '{print $2"_"$3"_"$9}' ../Nematode/multiplicon.retry |while read line; do grep -w $line gene.positions.1 ;done |awk '{print $1,$2}' ) <(awk '{print $2"_"$3"_"$10}' ../Nematode/multiplicon.retry|while read line; do grep -w $line gene.positions.1 ;done| awk '{print $1,$3}') <(awk '{print $4"_"$5"_"$11}' ../Nematode/multiplicon.retry|while read line; do grep -w $line gene.positions.1 ;done | awk '{print $1,$2}' ) <(awk '{print $4"_"$5"_"$12}' ../Nematode/multiplicon.retry |while read line; do grep -w $line gene.positions.1 ;done |awk '{print $1,$3}' )|sed 's/_/ /g'|awk '{print $2,$4,$8,$10,$12,$16}'|sort|uniq>syntenic.ribbons.txt
awk '{print $1}' syntenic.ribbons.txt |uniq|cdbyank genome738sl.polished.mitoFixed.noquiver.fa.cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" green"}' >kary.redo
cat kary.redo karyotype.both.again.subset.roundup >kary.both
```
  ![Circos C. elegans](assets/CircosCelegans.png)

```
When C. elegans was updated with familial orthologous calculations, no difference was found. This multiplicon below, was the only synteny found. This was less than B. xylophilus, but could represent the high contiguity of the C. elegans genome.
000457 64162 85295 I 8606404 8616420
8 SCN 000457 CEL I 2 3 7 15 21 1107 1111 0

Ha ha, between 7 genes in scn, and 5 genes in C elegans, with three anchors. very small region.

```
### Globodera pallida
```


###/data021/GIF/remkv6/Baum/CamTechGenomeComparison/21_globoderapSCNsynteny
This is the closest related nematode to H. glycines, so if synteny exists, I should find synteny here


[2016-10-07 09:48:04] wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS7/species/globodera_pallida/PRJEB123/globodera_pallida.PRJEB123.WBPS7.genomic_softmasked.fa.gz
[2016-10-07 09:48:31] wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS7/species/globodera_pallida/PRJEB123/globodera_pallida.PRJEB123.WBPS7.annotations.gff3.gz
[2016-10-07 09:48:42] wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS7/species/globodera_pallida/PRJEB123/globodera_pallida.PRJEB123.WBPS7.protein.fa.gz
[2016-10-07 09:48:53] gunzip *
##create opscan input files

t" $7 "\t" "t"}' |sed 's/pathogens_Gpal_scaffold_//g' |sort -k 5 -V|awk '{print ">"$2 "\t" "1" "\t" $4 "\t" $5 "\t" $6 "\t" $3}'|sort -k 1 -V |sed 's/gene://g' >gpallida.join1
 tr "\n" "\t" <globodera_pallida.PRJEB123.WBPS7.protein.fa |sed 's/>/\n>/g'|tr " " "\t" | cut  -f 1,4-|sed 's/\t/%/1'|sed 's/\t//g'|sed 's/%/\t/g' >gpallida.join2
join -1 1 -2 1 <(sort gpallida.join1) <(sort gpallida.join2) >gpallida.join3
sed 's/ /\n/6' gpallida.join3 >genome.db
##files for the 738 genome already existed and were used.
##move to nematode subdirectory
###/data021/GIF/remkv6/Baum/CamTechGenomeComparison/21_globoderapSCNsynteny/Nematode
cd ../Nematode
cp ../opscan/iADHoRe_PLAZA.table .
cp ../opscan/genome.db .
cp ../opscan/operon.db .
cd query
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini

cd ../subject/

<note important>#not sure if the tab is correct.  It could be a space that is needed between $5 and $6.important</note>
tr "\n" "\t" < ../genome.db |sed 's/>/\n>/g'|awk '{print $1$5 "\t" $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini

cd ../
cat query/query.ini subject/subject.ini >tmp
tr "\t" " "<tmp >Nematode.ini

###modified the Nematode.ini file as follows
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

Results of Globodera pallida synteny

Yes, there appears to be synteny conserved between these two species. It is comprised of 341 multiplicons over ~1600 genes. AlignmentMultiplicon3.svg contains a small inversion.
Further investigation is needed…
Andrew has requested to know if globodera can stitch Heterodera scaffolds together.
Redo of Synteny with Globodera using re-brakered masked genome

###/data021/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny


cp -rf ../21_globoderapSCNsynteny/opscan/ .
cp -rf ../21_globoderapSCNsynteny/Nematode/ . &
cd opscan
#cleaned out old 738 genome files so that the new ones can be set up.
ln -s ../../22_RepeatModeler/braker/braker/genome738sl.polished.mitoFixed.noquiver.masked/augustus.aa 738.polished.mitofixed.repmod.aa
ln -s ../../22_RepeatModeler/braker/braker/genome738sl.polished.mitoFixed.noquiver.masked/augustus.gff3 738.polished.mitofixed.repmod.gff3

 sed -e 's/ID=//g' 738.polished.mitofixed.repmod.gff3 -e 's/;/\t/g' -e 's/NODE_//g'|grep gene -|awk '{print $3 "\t" $9 "\t" substr($1, 0,8) "\t" $4 "\t" $5 "\t" $7 "\t" "t"}' |sed 's/scaffold_//g' |sort -k 5 -V|awk '{print ">"$2 "\t" "1" "\t" $4 "\t" $5 "\t" $6 "\t" $3}'|sort -k 1 -V |sed -e 's/F+//g' -e 's/F-//g' -e 's/F//g' >augustus.738.gff3.reformat1
 tr "\n" "\t" < 738.polished.mitofixed.repmod.aa |sed 's/>/\n>/g'|sed 's/\t//g'|grep ".t1"|sed  's/.t1/\t/g'>7384join2
join -1 1 -2 1 augustus.738.gff3.reformat1 7384join2 >joined.738.3
sed 's/ /\n/6' joined.738.3 >operon.db

#Globodera was already set up, so did not have to do anything.
qsub opscan.pbs

cd ../Nematode
#cleaned out old files from Nematode directory, keeping only iadhore.pbs, query and output folders.
cp ../opscan/operon.db .
cp ../opscan/genome.db .
cp ../opscan/iADHoRe_PLAZA.table .
cd query
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini

cd ../subject/
tr "\n" "\t" < ../genome.db |sed 's/>/\n>/g'|awk '{print $1$5 "\t" $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini

cd ../
cat query/query.ini subject/subject.ini >tmp
tr "\t" " "<tmp >Nematode.ini
###modified the Nematode.ini file as follows
genome=Heterodera
#sequence stuff
#need this empty newline
genome=Globodera
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

```
CIRCOS VISUALIZATION
```
#/data021/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos
Some files needed for visualization are generic. Be sure to copy them to the circos directory.
This includes: bands.conf, ideogram.conf, and ticks.conf.
#generate karyotype file

Gets scaffold names and lengths for karyotype file.
awk '{print $1}' trait* |sort|uniq|cdbyank ../genome738sl.polished.mitoFixed.noquiver.fa.cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" green"}' >H.glycines.kary
awk '{print $1"\n"$4}' syntenic.ribbons.txt |sort|uniq|cdbyank globodera_pallida.PRJEB123.WBPS7.genomic_softmasked.renamed.fa.cidx |bioawk -c fastx '{print "chr - "$name" "$name"
0 " length($seq)" green"}' >globo.kary1
cat hetero.kary1 globo.kary1 >karyotype.both.txt
karyotype.both.txt – just a subset to show the organization needed.

cat <(grep gene 738.polished.mitofixed.repmod.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}') <(grep gene globodera_pallida.PRJEB123.WBPS7.annotations.gff3 |sed 's/ID=gene://g'|sed 's/;Name.*//g'|grep -v "#" -| grep -v ID -|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}')  >both.gene.coord
cat <(awk '$2=="Heterodera" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../Nematode/output/genes.txt|sort -k 1,1 -V ) <(awk '$2=="GLOBODERA" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../Nematode/output/genes.txt|sort -k 1,1 -V) >reformat.genes.txt

join -a 1  -1 1 -2 1 <(sort -k 1,1 -V reformat.genes.txt) <(sort -k 1,1 -V both.gene.coord) >gene.positions
#note if there is a "_" in the scaffold or gene name, this does not work.  
sed 's/ /_/2' gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >gene.positions.1

paste <(awk '{print $2"_"$3"_"$9}' ../Nematode/output/multiplicons.txt |while read line; do grep -w $line gene.positions.1 ;done |awk '{print $1,$2}' ) <(awk '{print $2"_"$3"_"$10}' ../Nematode/output/multiplicons.txt|while read line; do grep -w $line gene.positions.1 ;done| awk '{print $1,$3}') <(awk '{print $4"_"$5"_"$11}' ../Nematode/output/multiplicons.txt|while read line; do grep -w $line gene.positions.1 ;done | awk '{print $1,$2}' ) <(awk '{print $4"_"$5"_"$12}' ../Nematode/output/multiplicons.txt |while read line; do grep -w $line gene.positions.1 ;done |awk '{print $1,$3}' )|sed 's/_/ /g'|awk '{print $2,$4,$8,$10,$12,$16}'>syntenic.ribbons.txt
sort -k 1,1 syntenic.ribbons.txt >test
cat <(awk '{print $1}' test) <(awk '{print $4}' test)|sort|uniq -c |sort -k 1,1 -nr|awk '$1!=1 {print $2}'|grep -w -f - test >syntenic.ribbons.subset
Create track of SL transcripts to align with synteny on circos plot.


mkdir splicedLeaders
cd splicedLeaders
awk '{print $3}' ../karyotype.both.again.subset.txt |cdbyank ../genome738sl.polished.mitoFixed.noquiver.fa.cidx|makeblastdb -in - -title test -dbtype nucl -out 738.subset.blastdb
blastn -db 738.subset.blastdb -query bothtestall.cdhit.fasta -outfmt '6 stitle sstart send' -out SLto738subset.blast.out
#convert to 100kb bin histogram for plot
awk '{print $1 "\t" $2/100000}' *.blast.out    | sed 's/\./\t/g' | awk '{print $1 "\t" $2}'| sort |uniq -c|awk '{print $2 "\t" $3*100000 "\t" $3*100000+100000 "\t" $1}' >SLto738subset.histogram
#Dont forget to add a plot to the circos.conf file to visualize.
Create a track for repetitive elements in the genome


#this is a gff generated by the repeatmasker program.
#converted to 1000bp bins
awk '{print $1,$4}' genome738sl.polished.mitoFixed.fa.out.gff|awk '{print $1 "\t" $2/1000}'  |sed 's/\./\t/g' |awk '{print $1,$2}'| sort |uniq -c|awk '{print $2 "\t" $3*1000 "\t" $3*1000+1000 "\t" $1}'>repeatmasker.histogram
#dont forget to add a plot to circos.conf

Circos conf file. Make sure to open edit to get rid of all the gibberish.
RESULT –> there are two large peaks on the plot that have lots of SL transcripts. –One of the peaks mapped to scaffold 000041, at 40-41kb. This hit sa gene of about 300bp in the unmasked braker run. This is not an SLRNA.
The only homology I could find was to a single G.max gene. Glyma.09G230300, which is annotated as an LRR gene. This gene is inverted in synteny between G. max and castor bean.
http://chibba.agtec.uga.edu/duplication/usr/locus/0af7dcbb2679dd948674071bdbbe42d7.html
The gene is anotated as a Clavata1 gene, which is involved in flower and meristem development. Arabidopsis mutants form extra floral organs. https://www.ncbi.nlm.nih.gov/pubmed/8287795
I bet that this gene is involved in the development of the syncytium.
This analysis was left in the spliced leader directory.
/data021/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/splicedLeaders

Cripes those ribbons are everywhere, lets straighten them out with a handy dandy circos tool.
?
1
2
3
4

#figure out which chromosomes/scaffolds that you want to keep in order.
awk '{print $3}' hetero.kary1|tr "\n" ","|less
#there were some problems match contigs, as iadhore truncated syntmer to syntme, and the original fasta still contains "|quiver".  Make sure everything matches across the synteny and karyotype file.
/shared/software/GIF/programs/circos/0.69.2/../circos-tools-0.22/tools/orderchr/bin/orderchr -links syntenic.ribbons.subset -karyotype karyotype.both.again.subset.txt -init_order  000001,000002,000011,000014,000015,000019,000021,000028,000037,000041,000049,000058,000059,000066,000117,000118,000138,000161,000171,000177,000190,000195,000212,000217,000220,000251,000263,000286,000300,000309,000348,000358K,000404,000502,000540,000614K,000640K,000928,000975K,001124,002293,002368,20syntmer,22syntmer,34syntmer,35syntmer,55syntmer -static_rx 000001,000002,000011,000014,000015,000019,000021,000028,000037,000041,000049,000058,000059,000066,000117,000118,000138,000161,000171,000177,000190,000195,000212,000217,000220,000251,000263,000286,000300,000309,000348,000358K,000404,000502,000540,000614K,000640K,000928,000975K,001124,002293,002368,20syntmer,22syntmer,34syntmer,35syntmer,55syntmer

awk '($6-$5)>70000 {print $1"\n"$4}' syntenic.ribbons.txt |sort|uniq|cdbyank genome738sl.polished.mitoFixed.noquiver.fa.cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" green"}' >H.glycines.kary
awk '($6-$5)>70000 {print $1"\n"$4}' syntenic.ribbons.txt |sort|uniq|cdbyank globodera_pallida.PRJEB123.WBPS7.genomic_softmasked.renamed.fa.cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" green"}' >globo.kary1
cat hetero.kary1 globo.kary1 >karyotype.both.txt
```
  ![Circos C. elegans](assets/CircosGpallida.png)

### Globodera rostochiensis
```
Setting up opscan


wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/079/975/GCA_900079975.1_nGr/GCA_900079975.1_nGr_genomic.fna.gz
gunzip *
perl ~/common_scripts/gff2fasta.pl GCA_900079975.1_nGr_genomic.fna G.ros.gff3 G.ros
[2016-12-29 07:01:08] awk '$3=="gene"' G.ros.gff3 |sed 's/ID=/>/g'|sed 's/;.*//g' |awk '{print $9,"1",$4,$5,$7,$1}' |paste - <(bioawk -c fastx '{print $name, $seq}' G.ros.pep.fasta) |awk '{print $1,$2,$3,$4,$5,$6,$8}' |sed 's/ /\n/6' >genome.db

#copy the operon.db from pallida folder
#copy opscan.pbs from pallida folder
qsub opscan.pbs

Setup for iadhore


cat <(grep "CL" iADHoRe_PLAZA.pretable|awk '{print $3,"Family"$2"\n"$6,"Family"$2}') <(grep "not" iADHoRe_PLAZA.pretable|awk '{print $2}') <(grep "SG" iADHoRe_PLAZA.pretable|awk '{print $2,"Family"NR+10000"\n"$5,"Family"NR+10000}' ) <(grep ">" genome.db|grep -v -f <(grep "GROS_g" iADHoRe_PLAZA.pretable |awk '{print $6}' |grep "frz3" - )|sed 's/>//g'|awk '{print $1}' )|sort -k 1,1V |uniq| awk '{if(NF>1) print $0; else print $1,"Family"NR+20000}'|tr " " "\t" > iADHoRe_PLAZA.table.families

mkdir subject
mkdir query
cd query
tr "\n" "\t" <../operon.db |sed 's/>/\n>/g'|awk '{print $1$5 "\t" $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini
cd ../subject/
tr "\n" "\t" < ../genome.db |sed 's/>/\n>/g'|awk '{print $1$5 "\t" $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
cd ../
cat query/query.ini subject/subject.ini >tmp
tr "\t" " "<tmp >Nematode.ini

#Parameters for Nematode.ini


prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=true
output_path= output
alignment_method=gg2
gap_size=15
cluster_gap=20
level_2_only=true
q_value=.9
blast_table=iADHoRe_PLAZA.table.families
table_type=family

#setting up circos


cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/bands.conf            .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/ticks.conf            .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/circos.test.delete    .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/housekeeping.conf     .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/ideogram.conf         .

cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/genome738sl.polished.mitoFixed.noquiver.fa.cidx .
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/genome738sl.polished.mitoFixed.noquiver.fa .
module load cdbfasta
module load bioawk

awk '$2=="SCN" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../output/genes.txt|sort -k 1,1 -V >Hgly.reformat.genes.txt
awk '$2!="SCN" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../output/genes.txt|sort -k 1,1 -V >subj.reformat.genes.txt
awk '$3=="gene"' ../../../circos/augustus.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}' >Hgly.gene.coord
awk '$3=="gene"' ../G.ros.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}' >subj.gene.coord
join -a 1  -1 1 -2 1 <(sort -k 1,1 -V Hgly.reformat.genes.txt) <(sort -k 1,1 -V Hgly.gene.coord) >Hgly.gene.positions
join -a 1  -1 1 -2 1 <(sort -k 1,1 -V subj.reformat.genes.txt) <(sort -k 1,1 -V subj.gene.coord) >subj.gene.positions
sed 's/ /_/2' Hgly.gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >Hgly.gene.positions.1
sed 's/ /</2' subj.gene.positions |sed 's/ /</2'|cut -d " " -f 1- |sed 's/</ /1' | sed 's/</_/g' |awk '{print $2"_"$3,$4,$5}' >subj.gene.positions.1

awk '{if ($2!="SCN") print $1"\t"$4"\t"$5"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8"\t"$11"\t"$12"\t"$9"\t"$10"\t"$12; else print $0}'  ../output/multiplicons.txt |awk '$4!=$2' >multiplicons.mod.txt

paste <(awk '{print $2"_"$3"_"$9}' multiplicons.mod.txt |while read line; do grep -w $line Hgly.gene.positions.1 ;done |awk '{print $1,$2}' ) <(awk '{print $2"_"$3"_"$10}' multiplicons.mod.txt|while read line; do grep -w $line Hgly.gene.positions.1 ;done| awk '{print $1,$3}') <(awk '{print $4"_"$5"_"$11}' multiplicons.mod.txt|while read line; do grep -w $line subj.gene.positions.1 ;done | awk '{print $1,$2}' |uniq ) <(awk '{print $4"_"$5"_"$12}' multiplicons.mod.txt |while read line; do grep -w $line subj.gene.positions.1 ;done |awk '{print $1,$3}' |uniq)|sed 's/_/ /g'|awk '{print $2,$4,$8,$10"_"$11,$13,$18}' |sed 's/syntme/syntmer/g' |sed 's/syntmerr/syntmer/g' >syntenic.ribbons.txt

awk '{print $1"\n"$4}' syntenic.ribbons.txt |sort |uniq|cdbyank ../*cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" blue"}' >subj.kary
awk '{print $1"\n"$4}' syntenic.ribbons.txt |sort |uniq|sed 's/syntme/syntmer/g' |sed 's/syntmerr/syntmer/g' |cdbyank genome738sl.polished.mitoFixed.noquiver.fa.cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" green"}' >H.glycines.kary
cat H.glycines.kary subj.kary >karyotype.both.txt
```
  ![Circos G. rostochiensis](assets/CircosGrostochiensis.png)

### Meloidogyne hapla
```
#opscan setup

wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/meloidogyne_hapla/PRJNA29083/meloidogyne_hapla.PRJNA29083.WBPS8.genomic.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/meloidogyne_hapla/PRJNA29083/meloidogyne_hapla.PRJNA29083.WBPS8.annotations.gff3.gz
gunzip *
perl ~/common_scripts/gff2fasta.pl meloidogyne_hapla.PRJNA29083.WBPS8.genomic.fa meloidogyne_hapla.PRJNA29083.WBPS8.annotations.gff3 meloidogyne_hapla

awk '$3=="gene"' meloidogyne_hapla.PRJNA29083.WBPS8.annotations.gff3 |sed 's/;/\t/g' |sed 's/ID=gene:/>/g' |awk '{print $9,"1",$4,$5,$7,$1}' |paste - <(bioawk -c fastx '{print $name, $seq}' meloidogyne_hapla.pep.fasta) |awk '{print $1,$2,$3,$4,$5,$6,$8}' |sed 's/ /\n/6' >genome.db
cp ../../26_GloboderaSynteny/otherGloboderaGenomes/rostochiensis/operon.db .
#copy opscan.pbs to current dir
qsub opscan.pbs

#iadhore setup

cd query
tr "\n" "\t" <../operon.db |sed 's/>/\n>/g'|awk '{print $1$5 "\t" $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini
cd ../subject/
<note important>#not sure if the tab is correct.  It could be a space that is needed between $5 and $6.important</note>
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
cd ../
cat query/query.ini subject/subject.ini >tmp
tr "\t" " "<tmp >Nematode.ini

grep "MhA1_Contig" iADHoRe_PLAZA.pretable |awk '{print $6}' |grep "MhA1_Contig" -  > remove.GPLIN.list
44949  [2016-12-29 10:58:06] cat <(grep "CL" iADHoRe_PLAZA.pretable|awk '{print $3,"Family"$2"\n"$6,"Family"$2}') <(grep "not" iADHoRe_PLAZA.pretable|awk '{print $2}') <(grep "SG" iADHoRe_PLAZA.pretable|awk '{print $2,"Family"NR+10000"\n"$5,"Family"NR+10000}' ) <(grep ">" genome.db|grep -w -v -f remove.GPLIN.list - |sed 's/>//g'|awk '{print $1}' )|sort -k 1,1V |uniq| awk '{if(NF>1) print $0; else print $1,"Family"NR+20000}'|tr " " "\t"  >iADHoRe_PLAZA.table.families

#Nematode.ini parameters

prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=true
output_path= output
alignment_method=gg2
gap_size=15
cluster_gap=20
level_2_only=true
q_value=.9
blast_table=iADHoRe_PLAZA.table.families
table_type=family

#Circos Setup

cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/bands.conf            .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/ticks.conf            .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/circos.test.delete    .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/housekeeping.conf     .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/ideogram.conf         .

cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/genome738sl.polished.mitoFixed.noquiver.fa.cidx .
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/genome738sl.polished.mitoFixed.noquiver.fa .
cdbfasta ../meloidogyne_hapla.PRJNA29083.WBPS8.genomic.fa
module load cdbfasta
module load bioawk

awk '$3=="gene"'  /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/augustus.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}' >Hgly.gene.coord
awk '$3=="gene"' ../meloidogyne_hapla.PRJNA29083.WBPS8.annotations.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|sed 's/Name/\t/g' |sed 's/gene://g' |awk '{print $1 "\t" $6 "\t" $7}' >subj.gene.coord
join -a 1  -1 1 -2 1 <(sort -k 1,1 -V Hgly.reformat.genes.txt) <(sort -k 1,1 -V Hgly.gene.coord) >Hgly.gene.positions
join -a 1  -1 1 -2 1 <(sort -k 1,1 -V subj.reformat.genes.txt) <(sort -k 1,1 -V subj.gene.coord) >subj.gene.positions
sed 's/ /_/2' Hgly.gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >Hgly.gene.positions.1
sed 's/ /_/2' subj.gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >subj.gene.positions.1


awk '{if ($2!="SCN") print $1"\t"$4"\t"$5"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8"\t"$11"\t"$12"\t"$9"\t"$10"\t"$12; else print $0}'  ../output/multiplicons.txt |awk '$4!=$2' >multiplicons.mod.txt

paste <(awk '{print $2"_"$3"_"$9}' multiplicons.mod.txt |while read line; do grep -w $line Hgly.gene.positions.1 ;done |awk '{print $1,$2}' ) <(awk '{print $2"_"$3"_"$10}' multiplicons.mod.txt|while read line; do grep -w $line Hgly.gene.positions.1 ;done| awk '{print $1,$3}') <(awk '{print $4"_"$5"_"$11}' multiplicons.mod.txt|while read line; do grep -w $line subj.gene.positions.1 ;done | awk '{print $1,$2}' |uniq ) <(awk '{print $4"_"$5"_"$12}' multiplicons.mod.txt |while read line; do grep -w $line subj.gene.positions.1 ;done |awk '{print $1,$3}' |uniq) |sed -e 's/_/ /1' -e 's/_/ /1' -e 's/_/ /1' -e 's/_/ /1'  -e 's/_/ /1' -e 's/_/ /1' -e 's/_/ /1' -e 's/_/ /1'|sed 's/ /_/8' |sed 's/_/ /3' | awk '{print $2,$4,$8,$10,$12,$16 }' >syntenic.ribbons.txt

awk '{print $1"\n"$4}' syntenic.ribbons.txt |sort |uniq|cdbyank ../*cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" blue"}' >subj.kary
awk '{print $1"\n"$4}' syntenic.ribbons.txt |sort |uniq|sed 's/syntme/syntmer/g' |sed 's/syntmerr/syntmer/g' |cdbyank genome738sl.polished.mitoFixed.noquiver.fa.cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" green"}' >H.glycines.kary
cat H.glycines.kary subj.kary > karyotype.both.txt
```  
  ![Circos M. hapla](assets/CircosMhapla.png)


### Meloidogyne incognita

```
opscan setup


44497  [2016-12-29 05:01:08] wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/meloidogyne_incognita/PRJEA28837/meloidogyne_incognita.PRJEA28837.WBPS8.genomic.fa.gz
44498  [2016-12-29 05:01:17] wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/meloidogyne_incognita/PRJEA28837/meloidogyne_incognita.PRJEA28837.WBPS8.annotations.gff3.gz

awk '$3=="gene"' meloidogyne_incognita.PRJEA28837.WBPS8.annotations.gff3 |sed 's/;/\t/g' |sed 's/ID=gene:/>/g' |awk '{print $9,"1",$4,$5,$7,$1}' |paste - <(bioawk -c fastx '{print $name, $seq}' meloidogyne_incognita.PRJEA28837.WBPS8.protein.fa) |awk '{print $1,$2,$3,$4,$5,$6,$8}' |sed 's/ /\n/6' >genome.db
#copy operon.db to here
#copy opscan.pbs to here
qsub opscan.pbs
iadhore setup


cat <(grep "CL" iADHoRe_PLAZA.pretable|awk '{print $3,"Family"$2"\n"$6,"Family"$2}') <(grep "not" iADHoRe_PLAZA.pretable|awk '{print $2}') <(grep "SG" iADHoRe_PLAZA.pretable|awk
'{print $2,"Family"NR+10000"\n"$5,"Family"NR+10000}' ) <(grep ">" genome.db|grep -v -f <(grep "Minc" iADHoRe_PLAZA.pretable |awk '{print $6}' |grep "Minc" - )|sed 's/>//g'|awk '{print $1}' )|sort -k 1,1V |uniq| awk '{if(NF>1) print $0; else print $1,"Family"NR+20000}'|tr " " "\t" > iADHoRe_PLAZA.table.families

cd query
tr "\n" "\t" <../operon.db |sed 's/>/\n>/g'|awk '{print $1$5 "\t" $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini
cd ../subject/
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
cd ../
cat query/query.ini subject/subject.ini >tmp
tr "\t" " "<tmp >Nematode.ini
Nematode.ini setup


prob_cutoff=0.001
anchor_points=3
number_of_threads=16
visualizeAlignment=true
output_path= output
alignment_method=gg2
gap_size=15
cluster_gap=20
level_2_only=true
q_value=.9
blast_table=iADHoRe_PLAZA.table.families
table_type=family
#circos setup

cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/bands.conf            .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/ticks.conf            .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/circos.test.delete    .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/housekeeping.conf     .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/ideogram.conf         .

cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/genome738sl.polished.mitoFixed.noquiver.fa.cidx .
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/genome738sl.polished.mitoFixed.noquiver.fa .
cdbfasta ../meloidogyne_incognita.PRJEA28837.WBPS8.genomic.fa
module load cdbfasta
module load bioawk

less ../output/multiplicons.txt |awk '{print $3"\n"$6}' |sort|uniq | cdbyank genome738sl.polished.mitoFixed.noquiver.fa.cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" green"}' >H.glycines.kary
less ../output/multiplicons.txt |awk '{print $3"\n"$6}' |sort|uniq | cdbyank ../*cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 " length($seq)" blue"}' >subj.kary
cat H.glycines.kary subj.kary >karyotype.both.txt
awk '$2=="SCN" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../output/genes.txt|sort -k 1,1 -V >Hgly.reformat.genes.txt
awk '$2!="SCN" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../output/genes.txt|sort -k 1,1 -V >subj.reformat.genes.txt
awk '$3=="gene"'  /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/augustus.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\
t" $5 "\t" $6}' >Hgly.gene.coord
awk '$3=="gene"' ../meloidogyne_incognita.PRJEA28837.WBPS8.annotations.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|less |awk '{print $1 "\t" $5 "\t" $6}' |sed 's/biotype=protein_coding//g'| sed 's/gene.*Name=//g' >subj.gene.coord
join -a 1  -1 1 -2 1 <(sort -k 1,1 -V Hgly.reformat.genes.txt) <(sort -k 1,1 -V Hgly.gene.coord) >Hgly.gene.positions
join -a 1  -1 1 -2 1 <(sort -k 1,1 -V subj.reformat.genes.txt) <(sort -k 1,1 -V subj.gene.coord) >subj.gene.positions
sed 's/ /_/2' Hgly.gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >Hgly.gene.positions.1
sed 's/ /_/2' subj.gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >subj.gene.positions.1

awk '{if ($2!="SCN") print $1"\t"$4"\t"$5"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8"\t"$11"\t"$12"\t"$9"\t"$10"\t"$12; else print $0}'  ../output/multiplicons.txt |awk '$4!=$2' >multiplicons.mod.txt
paste <(awk '{print $2"_"$3"_"$9}' multiplicons.mod.txt |while read line; do grep -w $line Hgly.gene.positions.1 ;done |awk '{print $1,$2}' ) <(awk '{print $2"_"$3"_"$10}' multiplicons.mod.txt|while read line; do grep -w $line Hgly.gene.positions.1 ;done| awk '{print $1,$3}') <(awk '{print $4"_"$5"_"$11}' multiplicons.mod.txt|while read line; do grep -w $line subj.gene.positions.1 ;done | awk '{print $1,$2}' |uniq ) <(awk '{print $4"_"$5"_"$12}' multiplicons.mod.txt |while read line; do grep -w $line subj.gene.positions.1 ;done |awk '{print $1,$3}' |uniq) |sed 's/_/ /g' |awk '{print $2,$4,$8,$10,$12,$16}' >syntenic.ribbons.txt

awk '{print $1}' syntenic.ribbons.txt |sort|uniq |cdbyank genome738sl.polished.mitoFixed.noquiver.fa.cidx |bioawk -c fastx '{print "chr - "$name" "$name" 0 "length($seq)" green"}' |cat karyotype.both.txt - >tmp
mv tmp karyotype.both.txt
```

  ![Circos M. incognita](assets/CircosMincognita.png)

### Globodera ellingotnae

```
Globodera ellingotnae

Have to run braker, as there is not a gff3 available. #Braker setup

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/723/225/GCA_001723225.1_ASM172322v1/GCA_001723225.1_ASM172322v1_genomic.fna.gz
fastq-dump --split-files -O .  SRR3162514
gunzip *
cp ~/common_scripts/runBraker.sh .
#runBraker.pl

#!/bin/bash
# works on L3
# change this to condo is you need to run it there
# needs rnaseq reads as well as the genome to be annotated
# if you have multiple RNA seq libraries merge them together (all R1's and all R2's seperately)
if [ $# -lt 3 ] ; then
        echo "usage: runBraker.sh <RNAseq_R1> <RNAseq_R2> <genome.fasta>"
        echo ""
        echo "To align RNAseq reads to genome and run Braker gene prediction program"
        echo ""
exit 0
fi
module use /shared/software/GIF/modules
module load hisat2
module load braker/1.9

cp $GENEMARK_PATH/gm_key ~/.gm_key

R1="$1"
R2="$2"
GENOME="$3"
BASE=$(basename ${GENOME%.*})

hisat2-build ${GENOME} ${GENOME%.*}
hisat2 -p 15 -x ${GENOME%.*} -1 ${R1} -2 ${R2} | samtools view -bS - > ${BASE}_rnaseq.bam
samtools sort -m 5G ${BASE}_rnaseq.bam > ${BASE}_sorted_rnaseq.bam
braker.pl --cores=16 --overwrite --species=${BASE} --genome=${GENOME} --bam=${BASE}_sorted_rnaseq.bam --gff3

#PBS SCRIPT


#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=96:00:00
#PBS -N BrakerG.ell
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module purge
sh ./runBraker.sh ./SRR3162514_1.fastq ./SRR3162514_2.fastq ./GCA_001723225.1_ASM172322v1_genomic.fna

ase you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"
#preparing files for opscan


awk '$3=="gene"' augustus.gff3 |sed 's/;/\t/g' |sed 's/ID=gene://g' |awk '{print $9,"1",$4,$5,$7,$1}' |paste - <(bioawk -c fastx '{print $name, $seq}' augustus.aa) |awk '{print ">transcript:"$1,$2,$3,$4,$5,$6,$8}' |sed 's/ /\n/6' >genome.db
sed -i 's/>/>k/g' genome.db
#copy operon.db to here
#copy opscan.pbs to here
qsub opscan.pbs
#setup for iadhore


cat <(grep "CL" iADHoRe_PLAZA.pretable|awk '{print $3,"Family"$2"\n"$6,"Family"$2}') <(grep "not" iADHoRe_PLAZA.pretable|awk '{print $2}') <(grep "SG" iADHoRe_PLAZA.pretable|awk '{print $2,"Family"NR+10000"\n"$5,"Family"NR+10000}' ) <(grep ">" genome.db|grep -v -f <(grep "Minc" iADHoRe_PLAZA.pretable |awk '{print $6}' |grep "Minc" - )|sed 's/>//g'|awk '{print $1}' )|sort -k 1,1V |uniq| awk '{if(NF>1) print $0; else print $1,"Family"NR+20000}'|tr " " "\t" > iADHoRe_PLAZA.table.families

cd query
#FYI the next few steps can take a few minutes if there are lots of scaffolds to index
tr "\n" "\t" <../operon.db |sed 's/>/\n>/g'|awk '{print $1$5 "\t" $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "query/"$1}' input.txt)>query.ini
cd ../subject/
<note important>#not sure if the tab is correct.  It could be a space that is needed between $5 and $6.important</note>
tr "\n" "\t" < ../genome.db |sed 's/>/\n>/g'|awk '{print $1$5 "\t" $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls *lst >input.txt
paste <(cut -f 1 -d "." input.txt) <(awk '{print "subject/"$1}' input.txt)>subject.ini
cd ../
cat query/query.ini subject/subject.ini >tmp
tr "\t" " "<tmp >Nematode.ini
#setup for circos


cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/bands.conf            .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/ticks.conf            .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/circos.test.delete    .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/housekeeping.conf     .
cp /data013/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/ideogram.conf         .
module load cdbfasta
module load bioawk
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/genome738sl.polished.mitoFixed.noquiver.fa.cidx .
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos/genome738sl.polished.mitoFixed.noquiver.fa .
cdbfasta ../../../GCA_001723225.1_ASM172322v1_genomic.fna
../../../GCA_001723225.1_ASM172322v1_genomic.fna were indexed in file ../../../GCA_001723225.1_ASM172322v1_genomic.fna.cidx
ln -s ../../../GCA_001723225.1_ASM172322v1_genomic.fna
ln -s ../../../GCA_001723225.1_ASM172322v1_genomic.fna.cidx

awk '$2=="SCN" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../output/genes.txt|sort -k 1,1 -V >Hgly.reformat.genes.txt
awk '$2!="SCN" {print $1 "\t" $2 "\t" $3 "\t" $4}' ../output/genes.txt|sort -k 1,1 -V | awk 'NR>1' >subj.reformat.genes.txt
awk '$3=="gene"' ../../../../../circos/augustus.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}' >Hgly.gene.coord
awk '$3=="gene"' ../augustus.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}' |sed 's/g/kg/g'| awk 'NR>1'  >subj.gene.coord

join -a 1  -1 1 -2 1 <(sort -k 1,1 -V Hgly.reformat.genes.txt) <(sort -k 1,1 -V Hgly.gene.coord) >Hgly.gene.positions
join -a 1  -1 1 -2 1 <(sort -k 1,1 -V subj.reformat.genes.txt) <(sort -k 1,1 -V subj.gene.coord) >subj.gene.positions
sed 's/ /_/2' Hgly.gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >Hgly.gene.positions.1
sed 's/ /_/2' subj.gene.positions |sed 's/ /_/2'|cut -d " " -f 2- >subj.gene.positions.1

awk '{if ($2!="SCN") print $1"\t"$4"\t"$5"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8"\t"$11"\t"$12"\t"$9"\t"$10"\t"$12; else print $0}'  ../output/multiplicons.txt |awk '$4!=$2' >multiplicons.mod.txt
paste <(awk '{print $2"_"$3"_"$9}' multiplicons.mod.txt |while read line; do grep -w $line Hgly.gene.positions.1 ;done |awk '{print $1,$2}' ) <(awk '{print $2"_"$3"_"$10}' multiplicons.mod.txt|while read line; do grep -w $line Hgly.gene.positions.1 ;done| awk '{print $1,$3}') <(awk '{print $4"_"$5"_"$11}' multiplicons.mod.txt|while read line; do grep -w $line subj.gene.positions.1 ;done | awk '{print $1,$2}' |uniq ) <(awk '{print $4"_"$5"_"$12}' multiplicons.mod.txt |while read line; do grep -w $line subj.gene.positions.1 ;done |awk '{print $1,$3}' |uniq)|sed 's/_/ /g'|awk '{print $2,$4,$8,$10,$12,$16}' |sed 's/syntme/syntmer/g' |sed 's/syntmerr/syntmer/g' >syntenic.ribbons.txt
```
  ![Circos G. ellingtonae](assets/CircosGellingtonae.png)

### Bursephelenchus xylophilus
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/51_B.xylophilus


sed 's/;/\t/g' B.xylophilus.gff3 |awk '$3=="CDS"' |sed 's/ID=gene://g' |awk '{print $9,"1",$4,$5,$7,$1}' | paste - <(bioawk -c fastx '{print $name, $seq}' B.xylophilus.pep.fasta ) |awk '{print $1,$2,$3,$4,$5,$6,$8}' |sed 's/ /\n/6' >genome.db


1SynChro/Opscan/bin/opscan -b -Q -C -O operon.db genome.db >iADHoRe_PLAZA.pretable
cat <(grep "CL" iADHoRe_PLAZA.pretable|awk '{print $3,"Family"$2"\n"$6,"Family"$2}') <(grep "not" iADHoRe_PLAZA.pretable|awk '{print $2}') <(grep "SG" iADHoRe_PLAZA.pretable|awk '{print $2,"Family"NR+10000"\n"$5,"Family"NR+10000}' ) <(grep ">" genome.db|grep -v -f <(grep "BXY" iADHoRe_PLAZA.pretable |awk '{print $6}' |grep "BXY" - )|sed 's/>//g'|awk '{print $1}' )|sort -k 1,1V |uniq| awk '{if(NF>1) print $0; else print $1,"Family"NR+20000}'|tr " " "\t" > iADHoRe_PLAZA.table.families
```
No multiplicons between H. glycines and B.xylophilus

##  Stats post processing
```


Need to compile general stats for orthologues, syntenic region counts, sizes with summary.sh.
G.pallida

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/Nematode/output


# multiplicons without internal synteny to H.glycines
less multiplicons.txt |awk '$2!=$4 && NR!=1' |wc
    388    5044   21887

#the largest multiplicon was 85 genes wide and between scaffold 17 in G.pallida and 000263 in H.glycines.
less multiplicons.txt |awk 'NR!=1' |sort -k8,8nr |awk ' NR==1'
19      GLOBODERA       17              Heterodera      000263  2       11      85      23      52      21      104     0
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/circos


grep "000263" syntenic.ribbons.txt
19:17 196482 351274 000263 146745 567547
41:48 35609 176255 000263 1042093 1349159
58:40 68883 129409 000263 10057 82556
382:48 210650 260985 000263 796869 842705

#size of largest syntenic region in each
grep "000263" syntenic.ribbons.txt |awk 'NR==1 {print ($3-$2),($6-$5)}'
154792 420802
(size more than twice as wide in H. glycines.)

#a difference of bp
grep "000263" syntenic.ribbons.txt |awk 'NR==1 {print ($6-$5)-($3-$2)}'
266010
#how much more distance there is in H.glycines.
grep "000263" syntenic.ribbons.txt |awk 'NR==1 {print ($6-$5)/($3-$2)}'
2.7185

#total syntenic distance in G. pallida
less syntenic.ribbons.txt |awk '$4!=$1' |awk '{if(length($4)>5) print $4,$5,$6,$1,$2,$3; else print $0}' |awk '{if($6>$5) print $6-$5; else print $5-$6}' |~/common_scripts/summary.sh
Total:  13,609,472
Count:  400
Mean:   34,023
Median: 25,252
Min:    2,669
Max:    251,295


#total syntenic distance in H. glycines
[remkv6@condo2017 circos]$ less syntenic.ribbons.txt |awk '$4!=$1' |awk '{if(length($4)>5) print $4,$5,$6,$1,$2,$3; else print $0}' |awk '{if($3>$2) print $3-$2; else print $2-$3}' |~/common_scripts/summary.sh
Total:  16,477,006
Count:  400
Mean:   41,192
Median: 29,595
Min:    2,672
Max:    420,802

G. rostochiensis

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/otherGloboderaGenomes/rostochiensis


#of multiplicons
awk 'NR!=1' multiplicons.txt |wc
    454    5902   23607
 #count after removing self synteny
less multiplicons.txt |awk '$2!=$4 && NR!=1' |wc
    439    5707   22880
less multiplicons.txt |awk 'NR!=1' |sort -k8,8nr |awk ' NR==1'
1       Gros    GROS_00001              SCN     000001  2       98      170     33      184     11      159     0
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/otherGloboderaGenomes/rostochiensis/circos

#larger in scn
grep "000001" syntenic.ribbons.txt |awk 'NR==1' |awk 'NR==1 {print ($3-$2),($6-$5)}'
634642 568977
#~66kb larger
[remkv6@condo2017 circos]$ grep "000001" syntenic.ribbons.txt |awk 'NR==1' |awk 'NR==1 {print ($6-$5)-($3-$2)}'
-65665
#rostochiensis has 89% of H. glycines distance
[remkv6@condo2017 circos]$ grep "000001" syntenic.ribbons.txt |awk 'NR==1' |awk 'NR==1 {print ($3-$2)/($6-$5)}'
1.11541

G. rostochiensis syntenic distance
less syntenic.ribbons.txt |awk '$4!=$1' |awk '{if(length($4<7)) print $0; else print $4,$5,$6,$1,$2,$3}'|awk '{if($6>$5) print $6-$5; else print $5-$6}' |~/common_scripts/summary.sh
Total:  20,631,011
Count:  439
Mean:   46,995
Median: 29,258
Min:    0
Max:    568,977

H. glycines syntenic distance
less syntenic.ribbons.txt |awk '$4!=$1' |awk '{if(length($4<7)) print $0; else print $4,$5,$6,$1,$2,$3}'|awk '{if($3>$2) print $3-$2; else print $2-$3}' |~/common_scripts/summary.sh
Total:  23,210,928
Count:  439
Mean:   52,872
Median: 31,466
Min:    3,098
Max:    634,642

G. ellingtonae

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/otherGloboderaGenomes/ellingtonae/braker/GCA_001723225.1_ASM172322v1_genomic/output


#number of syntenic regions - self homology
less multiplicons.txt |awk '$2!=$4 && NR!=1' |wc
    362    4706   22650
less multiplicons.txt |awk 'NR!=1' |sort -k8,8nr |awk ' NR==1'
1       ellingotnae     MEIZ01000025            SCN     000502  2       90      175     4       158     53      207     0

H.glycines and G.ellingtonae length
grep "000502" syntenic.ribbons.txt |awk 'NR==1 {print ($3-$2),($6-$5)}'
601928 461543
#length diff
grep "000502" syntenic.ribbons.txt |awk 'NR==1 {print ($6-$5)-($3-$2)}'
-140385
#how many times larger
grep "000502" syntenic.ribbons.txt |awk 'NR==1 {print ($3-$2)/($6-$5)}'
1.30416

#syntenic stats for G. ellingtonae
less syntenic.ribbons.txt |awk '$4!=$1' |awk '{if(length($4<7)) print $0; else print $4,$5,$6,$1,$2,$3}'|awk '{if($6>$5) print $6-$5; else print $5-$6}' |~/common_scripts/summary.sh
Total:  24,121,193
Count:  362
Mean:   66,633
Median: 37,660
Min:    3,787
Max:    558,196

#syntenic stats for H. glycines
less syntenic.ribbons.txt |awk '$4!=$1' |awk '{if(length($4<7)) print $0; else print $4,$5,$6,$1,$2,$3}'|awk '{if($3>$2) print $3-$2; else print $2-$3}' |~/common_scripts/summary.sh
Total:  31,255,836
Count:  362
Mean:   86,342
Median: 48,260
Min:    4,117
Max:    751,156

M.hapla

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/39_Meloidogyne/hapla/output


awk '$2!=$4 && NR!=1' |wc
    112    1456    6285
less multiplicons.txt |awk 'NR!=1' |sort -k8,8nr |awk ' NR==1'
10      hapla   MhA1_Contig990          SCN     000011  2       6       24      2       20      78      95      0
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/39_Meloidogyne/hapla/circos


#largest distances, H. glycines and M. hapla
grep "000011" syntenic.ribbons.txt |awk 'NR==1 {print ($3-$2),($6-$5)}'
71846 40160
31kb shorter in M.hapla
grep "000011" syntenic.ribbons.txt |awk 'NR==1 {print ($6-$5)-($3-$2)}'
-31686
How much larger is H. glycines
grep "000011" syntenic.ribbons.txt |awk 'NR==1' |awk 'NR==1 {print ($3-$2)/($6-$5)}'
1.78899


H. glycines syntenic regions
 less syntenic.ribbons.txt |awk '$4!=$1' |awk '{if(length($4>5)) print $4,$5,$6,$1,$2,$3; else print $0}' |awk '{if($6>$5) print $6-$5; else print $5-$6}' |~/common_scripts/summary.sh
Total:  2,650,178
Count:  112
Mean:   23,662
Median: 17,085
Min:    2,933
Max:    105,074

M.hapla Syntenic regions
[remkv6@condo2017 circos]$ less syntenic.ribbons.txt |awk '$4!=$1' |awk '{if(length($4>5)) print $4,$5,$6,$1,$2,$3; else print $0}' |awk '{if($3>$2) print $3-$2; else print $2-$3}' |~/common_scripts/summary.sh
Total:  2,218,827
Count:  112
Mean:   19,810
Median: 18,014
Min:    2,887
Max:    64,572

M. incognita

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/39_Meloidogyne/incognita/output


less multiplicons.txt |awk '$2!=$4 && NR!=1' |wc
     15     195     828
less multiplicons.txt |awk '$2!=$4 && NR!=1' |sort -k8,8nr |awk ' NR==1'
2       incognita       MiV1ctg45               SCN     35syntme        2       8       36      1       25      34      67      0
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/39_Meloidogyne/incognita/circos


#length of longest synt region (profile length) Heterodera then incognita
grep "35syntme" syntenic.ribbons.txt |awk 'NR==1 {print ($3-$2),($6-$5)}'
109680 124997
#how much larger incognita
grep "35syntme" syntenic.ribbons.txt |awk 'NR==1 {print ($6-$5)-($3-$2)}'
15317
#Incognita is greater than 1.13 times larger
[remkv6@condo2017 circos]$ grep "35syntme" syntenic.ribbons.txt |awk 'NR==1 {print ($6-$5)/($3-$2)}'
1.13965


#H.glycines syntenic length
less syntenic.ribbons.txt |awk '$4!=$1' |awk '{if(length($4>5)) print $4,$5,$6,$1,$2,$3; else print $0}' |awk '{if($6>$5) print $6-$5; else print $5-$6}' |~/common_scripts/summary.sh
Total:  381,190
Count:  15
Mean:   25,412
Median: 16,872
Min:    5,341
Max:    109,680

#M. incognita syntenic length
less syntenic.ribbons.txt |awk '$4!=$1' |awk '{if(length($4>5)) print $4,$5,$6,$1,$2,$3; else print $0}' |awk '{if($3>$2) print $3-$2; else print $2-$3}' |~/common_scripts/summary.sh
Total:  511,585
Count:  15
Mean:   34,105
Median: 28,978
Min:    5,199
Max:    124,997
Looks like only M. incognita has a larger syntenic than H. glycines. These are most likely intergenic regions of interspersed repeat expansion and contraction etc. H. glycines is highly repetitive compared to the others it looks like.

C. elegans

A single syntenic region.
000457 64162 85295 I 8606404 8616420
8 SCN 000457 CEL I 2 3 7 15 21 1107 1111 0

#these are all hand calculations since there was just one.


H glycines
Total: 21,133
Count: 1
genes:7
C. elegans
Total: 10,016
Count: 1
genes: 5
#How much larger was the syntenic region in H. glycines?
11,117.

#How many times larger was the syntenic region in H. glycines
2.11

Orthologous Gene details

Analyzing the common orthologous genes between all of the species.


#how many families are there between SCN and GROS
less iADHoRe_PLAZA.table.families |awk '{print $2}' |sort| uniq -c |awk '$1!=1' |wc
   7556   15112  148533


grep "GROS" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|wc
grep "GROS" iADHoRe_PLAZA.table.families |awk '{print $1}' |sort|uniq|wc
grep "GROS" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|grep -w -f - iADHoRe_PLAZA.table.families |grep -v "GROS" |awk '{print $1}' |sort|uniq|wc

grep "MhA" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|wc
grep "MhA" iADHoRe_PLAZA.table.families |awk '{print $1}' |sort|uniq|wc
grep "MhA" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|grep -w -f - iADHoRe_PLAZA.table.families |grep -v "MhA" |awk '{print $1}' |sort|uniq|wc

grep "Minc" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|wc
grep "Minc" iADHoRe_PLAZA.table.families |awk '{print $1}' |sort|uniq|wc
grep "Minc" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|grep -w -f - iADHoRe_PLAZA.table.families |grep -v "Minc" |awk '{print $1}' |sort|uniq|wc

grep "GPLIN" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|wc
grep "GPLIN" iADHoRe_PLAZA.table.families |awk '{print $1}' |sort|uniq|wc
grep "GPLIN" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|grep -w -f - iADHoRe_PLAZA.table.families |grep -v "GPLIN" |awk '{print $1}' |sort|uniq|wc

grep "kg" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|wc
grep "kg" iADHoRe_PLAZA.table.families |awk '{print $1}' |sort|uniq|wc
grep "kg" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|grep -w -f - iADHoRe_PLAZA.table.families |grep -v "kg" |awk '{print $1}' |sort|uniq|wc


grep "BXY" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|wc
20708
grep "BXY" iADHoRe_PLAZA.table.families |awk '{print $1}' |sort|uniq|wc
17704
grep "BXY" iADHoRe_PLAZA.table.families |awk '{print $2}' |sort|uniq|grep -w -f - iADHoRe_PLAZA.table.families |grep -v "BXY" |awk '{print $1}' |sort|uniq|wc
 4061
All combined calculations were done here #/work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny


#How many orthologous genes are there across Heterodera and all of Globodera?
cat opscan/iADHoRe_PLAZA.pretable otherGloboderaGenomes/rostochiensis/iADHoRe_PLAZA.pretable otherGloboderaGenomes/ellingtonae/braker/GCA_001723225.1_ASM172322v1_genomic/iADHoRe_PLAZA.pretable |grep "SG" |awk '{print $2}' |sort |uniq -c |sort -k1,1nr|awk '$1==3' |wc
   2723    5446   39624
#How many orthologous genes are there cross Heterodera and all of Meloidogyne?
cat  ../39_Meloidogyne/hapla/iADHoRe_PLAZA.pretable ../39_Meloidogyne/incognita/iADHoRe_PLAZA.pretable |grep "SG" |awk '{print $2}' |sort |uniq -c |sort -k1,1nr|awk '$1==2' |wc
   1900    3800   27623

#How many are conserved between Heterodera, Globodera, and Meloidogyne?   
cat opscan/iADHoRe_PLAZA.pretable otherGloboderaGenomes/rostochiensis/iADHoRe_PLAZA.pretable otherGloboderaGenomes/ellingtonae/braker/GCA_001723225.1_ASM172322v1_genomic/iADHoRe_PLAZA.pretable ../39_Meloidogyne/hapla/iADHoRe_PLAZA.pretable ../39_Meloidogyne/incognita/iADHoRe_PLAZA.pretable |grep "SG" |awk '{print $2}' |sort |uniq -c |sort -k1,1nr|awk '$1==5' |wc
    867    1734   12633

#How many are conserved across Heterodara, Globodera, Meloidogyne, and B. xylophilus
cat opscan/iADHoRe_PLAZA.pretable otherGloboderaGenomes/rostochiensis/iADHoRe_PLAZA.pretable otherGloboderaGenomes/ellingtonae/braker/GCA_001723225.1_ASM172322v1_genomic/iADHoRe_PLAZA.pretable ../39_Meloidogyne/hapla/iADHoRe_PLAZA.pretable ../39_Meloidogyne/incognita/iADHoRe_PLAZA.pretable ../51_B.xylophilus/iADHoRe_PLAZA.pretable|grep "SG" |awk '{print $2}' |sort |uniq -c |sort -k1,1nr|awk '$1==6' |wc
    602    1204    8763

#what about Heterodera, Globodera, Meloidogyne, and C. elegans
     Needed to rerun opscan with family settings
     #/work/GIF/remkv6/Baum/CamTechGenomeComparison/19_elegans3SCN/opscan
     ../../51_B.xylophilus/1SynChro/Opscan/bin/opscan -b -Q -C -O operon.db genome.db >iADHoRe_PLAZA.pretable &
cat opscan/iADHoRe_PLAZA.pretable otherGloboderaGenomes/rostochiensis/iADHoRe_PLAZA.pretable otherGloboderaGenomes/ellingtonae/braker/GCA_001723225.1_ASM172322v1_genomic/iADHoRe_PLAZA.pretable ../39_Meloidogyne/hapla/iADHoRe_PLAZA.pretable ../39_Meloidogyne/incognita/iADHoRe_PLAZA.pretable ../19_elegans3SCN/opscan/iADHoRe_PLAZA.pretable |grep "SG" |awk '{print $2}' |sort |uniq -c |sort -k1,1nr|awk '$1==6' |wc
     69     138    1005
#common orthologues B.xylophilus
cat ../51_B.xylophilus/iADHoRe_PLAZA.pretable|grep "SG" |awk '{print $2}' |sort |uniq -c |sort -k1,1nr|wc
3004

#common orthologues C. elegans
cat iADHoRe_PLAZA.pretable|grep "SG" |awk '{print $2}' |sort |uniq -c |sort -k1,1nr|wc
   2632    5264   38594
There seem to be very few. orthologues and lots of duplication in all of the other species… (mult genes in a family with 1 H.glycines gene)
```
  ![Synteny and Orthologs](assets/finalMitochondrialAnnotation.png)
  ![Synteny and Orthologs](assets/SpeciesGroupsOrthology.png)
  ![Synteny on phylogeny](assets/SyntenyPhylogeny.png)

### syntenic table stats
```
#orthologs in SCN when compared to G. pallida
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep -v "GPLIN" iADHoRe_PLAZA.table.families) |awk '{print $1}' |sort|uniq|wc
   6506    6506   42002
#orthologs in G. pallida when compared to SCN
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep -v "g" iADHoRe_PLAZA.table.families) |awk '{print $1}' |sort|uniq|wc
   4588    4588   73408
#Total orthologous families
 less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |awk '$1>1' |wc
   5845   11690  114533

#orthologs in SCN when compared to G. rostochiensis
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep -v "GROS" iADHoRe_PLAZA.table.families) |awk '{p
rint $1}' |sort|uniq|wc
   8180    8180   52732

#orthologs in G. rostochiensis when compared to SCN
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep "GROS" iADHoRe_PLAZA.table.families) |awk '{print $1}' |sort|uniq|wc
   6079    6079   72948

#Total orthologous families
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |awk '$1>1' |wc
   7556   15112  148533

#orthologs in SCN when compared to G. ellingtonae
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep -v "k" iADHoRe_PLAZA.table.families) |awk '{print $1}' |sort|uniq|wc
   8523    8523   54969


#orthologs in G. ellingtonae when compared to SCN
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep  "k" iADHoRe_PLAZA.table.families) |awk '{print $1}' |sort|uniq|wc
   6190    6190   46241

#Total orthologous families
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |awk '$1>1' |wc
   7771   15542  152729

#orthologs in SCN when compared to M. hapla
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep -v "M" iADHoRe_PLAZA.table.families) |awk '{print $1}' |sort|uniq|wc
   4553    4553   29352
#orthologs in M. hapla when compared to SCN
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep  "M" iADHoRe_PLAZA.table.families) |awk '{print $1}' |sort|uniq|wc
   3291    3291   88044

#Total orthologous families
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |awk '$1>1' |wc
   4117    8234   80578

#orthologs in SCN when compared to M. incognita
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep -v "M" iADHoRe_PLAZA.table.families) |awk '{print $1}' |sort|uniq|wc
   4148    4148   26781

#orthologs in M. incognita when compared to SCN
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep  "M" iADHoRe_PLAZA.table.families) |awk '{print $1}' |sort|uniq|wc
   3014    3014   30140
#Total orthologous families
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |awk '$1>1' |wc
   3775    7550   73868

#orthologs in SCN when compared to B. xylophilus
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep -v "B" iADHoRe_PLAZA.table.families) |awk '{print $1}' |sort|uniq|wc
   4061    4061   26180
#orthologs in B. xylophilus when compared to SCN
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep  "B" iADHoRe_PLAZA.table.families) |awk '{print $1}' |sort|uniq|wc
   3004    3004   36048
#Total orthologous families
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |awk '$1>1' |wc
   3727    7454   72984

#orthologs in SCN when compared to C. elegans
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep  -v "WB" iADHoRe_PLAZA.table.families) |awk '{print $1}
' |sort|uniq|wc
   3785    3785   25025

#orthologs in C. elegans when compared to SCN
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |sort -k 1,1nr |awk '$1>1' |awk '{print $2}' |grep -w  -f -  <(grep  "WB" iADHoRe_PLAZA.table.families) |awk '{print $1}' |
sort|uniq|wc
   2632    2632   39480

#Total orthologous families
less iADHoRe_PLAZA.table.families |cut -f 2 |sort|uniq -c |awk '$1>1' |wc
   3338    6676   65238
===GO enrichment in orthologous genes and non orthologous genes===
<sxh>
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/62_totalOrthologues/Ontologizer
ln -s ../../57_secretome/ontologenizer/simpleformat.ids
ln -s ../../57_secretome/ontologenizer/population
ln -s ../../57_secretome/ontologenizer/Ontologizer.jar
ln -s ../SCN.all.orthologues.list
wget http://purl.obolibrary.org/obo/go.obo
 module load java
java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s SCN.all.orthologues.list
java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s antiOrthologue.list

No significant enrichment was found.
```

###  Effector overlap
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/46_SCN_syntenic_tracks


#these are duplicate in the various species
less all.synteny.gff |awk '{print $5-$4}' |summary.sh
Total:  74,234,204
Count:  1,331
Mean:   55,773
Median: 31,670
Min:    2,672
Max:    751,156

#merged the coordinates from above
merge -i sorted.all.synteny.gff |awk '{print $1,"synteny","region",$2,$3,".","+","."}' |tr " " "\t" >Merged.All.gff

#Had to remove zeros
 vi Merged.All.gff

#how large is syteny
bedtools merge -i sorted.all.synteny.gff |awk '{print $3-$2}' |summary.sh
Total:  37,886,455
Count:  404
Mean:   93,778
Median: 53,520
Min:    4,578
Max:    932,890

How many genes are syntenic?
bedtools intersect -wo -a Merged.All.gff -b ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3|awk '$11=="gene" {print $17}' |wc
  10282   10282  109202

How many orthologs are syntenic?
bedtools intersect -wo -a Merged.All.gff -b ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3|awk '$11=="gene" {print $17}' |sed 's/ID=//g' |sed 's/;//g' |cat - ../62_totalOrthologues/SCN.all.orthologues.list |sort|uniq -c |awk '$1==2' |wc
3420

bedtools intersect -wo -a Merged.All.gff -b ../29_effectorMapping/effectors.gmapped.gff3 |awk '$11=="gene"' |sort -k9,9 -k12,13 |wc
     30     540    5623

30-7= 23 distinct effector loci, 21 of which are distinct effector types.
(GLAND 11,GLAND12,GLAND14, GLAND15, GLAND 17, GLAND 3, (GLAND 9 x2), 10C02,16A01,16H02,(17G06 x2), 19B10,19C07, 21E12 x3, 30C02,30D08 x2,3H07, 45D07, 4E02,4G06,8A07)

#which effectors are syntenic?
bedtools intersect -wo -a Merged.All.gff -b ../29_effectorMapping/effectors.gmapped.gff3 |awk '$11=="gene"' |sort -k9,9 -k12,13 |awk '{print $17}'|sed 's/\./\t/g' |awk '{print $1}'|sort|uniq >AllSyntenicEffectors.list

#which effectors are orthologous?
 #/work/GIF/remkv6/Baum/CamTechGenomeComparison/29_effectorMapping
bedtools intersect -wo -a effectors.gmapped.gff3 -b ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |awk '$3=="exon" && $15=="exon"' |awk '{print $21}' |sed 's/\./\t/g'|sed 's/ID=//g' |awk '{print $1}' |sort|uniq |cat - ../62_totalOrthologues/SCN.all.orthologues.list  |sort|uniq -c |awk '$1==2' |awk '{print $2}' |grep -f - <(bedtools intersect -wo -a effectors.gmapped.gff3 -b ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |awk '$3=="exon" && $15=="exon"') |awk '{print $9}' |sed 's/;/\t/g'|sed 's/\./\t/g' |awk '{print $1}' |sort|uniq >../62_totalOrthologues/OrthologousEffectors.list

#which effectors are both orthologous and syntenic?
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/62_totalOrthologues
cat OrthologousEffectors.list ../46_SCN_syntenic_tracks/AllSyntenicEffectors.list |sort|uniq -c |sort -k1,1nr|less
 2 ID=GLAND12|Pioneer
      2 ID=GLAND15|Pioneer
      2 ID=GLAND3|Pioneer,12H04family
      2 ID=GLAND9|Pioneer
      2 ID=lcl|flhggfha16A01|Pioneer,30D08/21E12family
      2 ID=lcl|flhggfha17G06|Pioneer
      2 ID=lcl|flhggfha21E12|Pioneer,30D08/16A01family
      2 ID=lcl|flhggfha30C02|Pioneer
      2 ID=lcl|flhggfha30D08|Pioneer,16A01/21E12family
      2 ID=lcl|flhggfha3H07|Ubiquitinextension
      2 ID=lcl|flhggfha45D07|ChorismateMutase
      2 ID=lcl|flhggfha4G06|Ubiquitinextensionprotein
      2 ID=lcl|flhggfha8A07|Pioneer
```
