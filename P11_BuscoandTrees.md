#  Need to run BUSCO 2 on genome to see how complete it is

```
Busco/2.0 was released with a nematode specific set of genes, so wanted to test it out.


module load busco/2.0
python3 /shared/software/GIF/programs/busco/2.0/BUSCO.py -i genome738sl.polished.mitoFixed.fa -l /shared/software/GIF/programs/busco/2.0/nematoda_odb9/ -o test1 -c 16 -m geno
python3 /shared/software/GIF/programs/busco/2.0/BUSCO.py -i Trinity.fasta -l /shared/software/GIF/programs/busco/2.0/nematoda_odb9/ -o testgland -c 16 -m tran
python3 /shared/software/GIF/programs/busco/2.0/BUSCO.py -i test.pep.fasta -l /shared/software/GIF/programs/busco/2.0/nematoda_odb9/ -o testpep -c 16 -m prot
Genome busco results.


#genome
Results:
C:54.4%[S:45.7%,D:8.7%],F:10.3%,M:35.3%,n:982
534 Complete BUSCOs (C)
449 Complete and single-copy BUSCOs (S)
85 Complete and duplicated BUSCOs (D)
101 Fragmented BUSCOs (F)
347 Missing BUSCOs (M)
982 Total BUSCO groups

#gene models cdna
INFO    Results:
INFO    C:70.5%[S:49.7%,D:20.8%],F:7.3%,M:22.2%,n:982
INFO    692 Complete BUSCOs (C)
INFO    488 Complete and single-copy BUSCOs (S)
INFO    204 Complete and duplicated BUSCOs (D)
INFO    72 Fragmented BUSCOs (F)
INFO    218 Missing BUSCOs (M)
INFO    982 Total BUSCO groups searched

#gene models amino acids
Results:
INFO    C:71.9%[S:50.3%,D:21.6%],F:7.9%,M:20.2%,n:982
INFO    706 Complete BUSCOs (C)
INFO    494 Complete and single-copy BUSCOs (S)
INFO    212 Complete and duplicated BUSCOs (D)
INFO    78 Fragmented BUSCOs (F)
INFO    198 Missing BUSCOs (M)
INFO    982 Total BUSCO groups
It appears that Busco cannot deal with the funky splicing that SCN has, but if I allow braker to figure out the gene models first, then it works better.

Are these same buscos missing in globodera and meloidogyne?

No matter, if this splicing problem is conserved between Globodera or Meloidogyne, the same buscos should be difficult to find. I will see how much they overlap.


ln -s /work/GIF/remkv6/Baum/CamTechGenomeComparison/27_pyscafGlobodera/pyScaf/2692JGI/pyramidAll/G.ellingtonae.fa
ln -s /work/GIF/remkv6/Baum/CamTechGenomeComparison/27_pyscafGlobodera/pyScaf/2692JGI/pyramidAll/G.pallida.fa
ln -s /work/GIF/remkv6/Baum/CamTechGenomeComparison/27_pyscafGlobodera/pyScaf/2692JGI/pyramidAll/G.rostochiensis.fa
ln -s /work/GIF/remkv6/Baum/CamTechGenomeComparison/27_pyscafGlobodera/pyScaf/2692JGI/pyramidAll/M.hapla.fa
ln -s /work/GIF/remkv6/Baum/CamTechGenomeComparison/27_pyscafGlobodera/pyScaf/2692JGI/pyramidAll/M.incognita.fa
#this sets up all the genomes to run busco with the --long parameter which allows augustus to self train
ls -1 *.fa|while read line; do echo "python3 /shared/software/GIF/programs/busco/2.0/BUSCO.py -i "$line" -l /shared/software/GIF/programs/busco/2.0/nematoda_odb9/ -o "$line"_out -c 16 -m geno --long"; done >busco.sh
Ran this with JobR_condo.sh
when all is done, run this.
cat */missing*|grep -v “#” |sort|uniq -c |awk '$1<3' > nonunique.missing.tsv

BUSCO all nematodes for phylogenetic tree

So we switched to Condo2017 and some things changed in the analysis while programs were down. I ran busco on all of the related nematodes and found that all genomes had poorer scores than what was reported. Busco cannot predict the splicing patterns in the nematodes, and this means that an rna seq guided approach to identify the genes is necessary.
Approach now, extract gene models from all species, run BUSCO, make phylogenetic tree, ka/ks analysis.

Downloads Extractions


#softlinked all of the fa files and ran with --long parameters
for f in *fa; do echo "python3 /work/GIF/software/programs/busco/2.0/BUSCO.py -f -i "$f" -l /work/GIF/software/programs/busco/2.0/nematoda_odb9/ -o "$f".out -t "$f".tmp -c 16 -m geno --long"; done >busco.sh

Click to display ⇲

runbusco.sh

Click to display ⇲

Busco scores were really poor and errors found in Busco gene prediction, as seen earlier protein scores are much higher for SCN. So will plan to do all busco in protein


#/work/GIF/remkv6/Baum/CamTechGenomeComparison/47_busco_nematode/Busco_Prot
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/meloidogyne_hapla/PRJNA29083/meloidogyne_hapla.PRJNA29083.WBPS8.protein.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/globodera_pallida/PRJEB123/globodera_pallida.PRJEB123.WBPS8.protein.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/globodera_pallida/PRJEB123/globodera_pallida.PRJEB123.WBPS8.annotations.gff3.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/meloidogyne_hapla/PRJNA29083/meloidogyne_hapla.PRJNA29083.WBPS8.annotations.gff3.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/meloidogyne_incognita/PRJEA28837/meloidogyne_incognita.PRJEA28837.WBPS8.annotations.gff3.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/079/975/GCA_900079975.1_nGr/GCA_900079975.1_nGr_genomic.gff.gz
mv GCA_900079975.1_nGr_genomic.gff.gz G.rostochiensis.aa
mv G.rostochiensis.aa G.rostochiensis.gff
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/040/885/GCA_001040885.1_S_ratti_ED321/GCA_001040885.1_S_ratti_ED321_genomic.gff.gz
mv G.rostochiensis.gff G.rostochiensis.gff.gz

mv GCA_001040885.1_S_ratti_ED321_genomic.gff.gz S.ratti.gff.gz
mv globodera_pallida.PRJEB123.WBPS8.annotations.gff3.gz G.pallida.gff3.gz
mv meloidogyne_incognita.PRJEA28837.WBPS8.annotations.gff3.gz M.incognita.gff3.gz
mv meloidogyne_hapla.PRJNA29083.WBPS8.annotations.gff3.gz M.hapla.gff.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/meloidogyne_incognita/PRJEA28837/meloidogyne_incognita.PRJEA28837.WBPS8.genomic.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/meloidogyne_hapla/PRJNA29083/meloidogyne_hapla.PRJNA29083.WBPS8.genomic.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/globodera_pallida/PRJEB123/globodera_pallida.PRJEB123.WBPS8.genomic.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/040/885/GCA_001040885.1_S_ratti_ED321/GCA_001040885.1_S_ratti_ED321_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/079/975/GCA_900079975.1_nGr/GCA_900079975.1_nGr_genomic.fna.gz
gunzip *
cp ../../32_genePredictionComp/UnmaskedBraker/braker/genome738sl.polished.mitoFixed/sorted_unmasked_738mitofixed.gff3.gz .
ln -s ../../32_genePredictionComp/UnmaskedBraker/genome738sl.polished.mitoFixed.fa
ln -s ../../26_GloboderaSynteny/otherGloboderaGenomes/ellingtonae/GCA_001723225.1_ASM172322v1_genomic.fna G.ellingtonae.fa
ln -s ../../26_GloboderaSynteny/otherGloboderaGenomes/ellingtonae/braker/GCA_001723225.1_ASM172322v1_genomic/augustus.gff3 G.ellingtonae.gff3
gunzip sorted_unmasked_738mitofixed.gff3.gz
mv GCA_001040885.1_S_ratti_ED321_genomic.fna S.ratti.fa
mv GCA_900079975.1_nGr_genomic.fna G.rostochiensis.fa
mv globodera_pallida.PRJEB123.WBPS8.genomic.fa G.pallida.fa
mv meloidogyne_hapla.PRJNA29083.WBPS8.genomic.fa M.hapla.fa
mv meloidogyne_incognita.PRJEA28837.WBPS8.genomic.fa M.incognita.fa
<sxh>
All files downloaded and named appropriately
Converting to protein or downloading seqs.
<sxh>
paste <(ls -1 *fa) <(ls -1 *.gf*) |while read line; do echo "perl ~/common_scripts/gff2fasta.pl "$line" "$line; done >pepmaker.sh
sh pepmaker.sh
pepmaker.sh <hidden>


perl ~/common_scripts/gff2fasta.pl G.ellingtonae.fa G.ellingtonae.gff3 G.ellingtonae
perl ~/common_scripts/gff2fasta.pl G.pallida.fa G.pallida.gff3 G.pallida
perl ~/common_scripts/gff2fasta.pl G.rostochiensis.fa G.rostochiensis.gff G.rostochiensis
perl ~/common_scripts/gff2fasta.pl M.hapla.fa M.hapla.gff M.hapla.fa M.hapla
perl ~/common_scripts/gff2fasta.pl M.incognita.fa M.incognita.gff3 M.incognita
perl ~/common_scripts/gff2fasta.pl S.ratti.fa S.ratti.gff S.ratti
perl ~/common_scripts/gff2fasta.pl genome738sl.polished.mitoFixed.fa sorted_unmasked_738mitofixed.gff3 genome738sl.polished.mitoFixed

CONT.


rm sorted_unmasked_738mitofixed.gff3
ln -s ../../32_genePredictionComp/UnmaskedBraker/braker/genome738sl.polished.mitoFixed/augustus.gff3 genome738sl.polished.mitoFixed_ALLDATA.gff3
perl ~/common_scripts/gff2fasta.pl genome738sl.polished.mitoFixed.fa genome738sl.polished.mitoFixed_ALLDATA.gff3 test
####This gff was wonky and would not work with gff2fasta.  Easier to download existing records.
###perl  ~/common_scripts/gff2fasta.pl S.ratti.fa S.ratti.gff S.ratti
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/040/885/GCA_001040885.1_S_ratti_ED321/GCA_001040885.1_S_ratti_ED321_protein.faa.gz
gunzip GCA_001040885.1_S_ratti_ED321_protein.faa.gz
mv GCA_001040885.1_S_ratti_ED321_protein.faa S.ratti.pep
ln -s ../../26_GloboderaSynteny/otherGloboderaGenomes/rostochiensis/G.ros.gff3
perl ~/common_scripts/gff2fasta.pl G.rostochiensis.fa G.ros.gff3 G.rostochiensis
ln -s ../../26_GloboderaSynteny/otherGloboderaGenomes/rostochiensis/G.ros.pep.fasta
ln -s ../../26_GloboderaSynteny/otherGloboderaGenomes/rostochiensis/G.ros.gene.fasta

wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/bursaphelenchus_xylophilus/PRJEA64437/bursaphelenchus_xylophilus.PRJEA64437.WBPS8.genomic.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/bursaphelenchus_xylophilus/PRJEA64437/bursaphelenchus_xylophilus.PRJEA64437.WBPS8.annotations.gff3.gz
gunzip *
ls
mv bursaphelenchus_xylophilus.PRJEA64437.WBPS8.annotations.gff3 B.xylophilus.gff3
mv bursaphelenchus_xylophilus.PRJEA64437.WBPS8.genomic.fa B.xylophilus.fa
###this didnt work with gff2fasta, so downloaded theirs
###perl ~/common_scripts/gff2fasta.pl B.xylophilus.fa B.xylophilus.gff3 B.xylophilus
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS8/species/bursaphelenchus_xylophilus/PRJEA64437/bursaphelenchus_xylophilus.PRJEA64437.WBPS8.protein.fa.gz
gunzip bursaphelenchus_xylophilus.PRJEA64437.WBPS8.protein.fa.gz
mv bursaphelenchus_xylophilus.PRJEA64437.WBPS8.protein.fa B.xylophilus.pep.fasta
perl ~/common_scripts/gff2fasta.pl genome738sl.polished.mitoFixed.fa genome738sl.polished.mitoFixed_ALLDATA.gff3 cp G.ros.pep.fasta G.rostochiensis.pep.fastaSCN
unlink G.ros.pep.fasta
mv M.hapla.fa.pep.fasta M.hapla.pep.fasta
mv SCN.pep.fasta H.glycines.pep.fasta
```


### BUSCO run
using the protein sequences
```
ls -1 *pep.fasta |while read line;do echo "python3 BUSCO.py -f -i "$line" -l /work/GIF/software/programs/GIF/busco/2.0/nematoda_odb9/ -o "$line".out -t "$line".tmp -c 16 -m prot";done >busco.sh

python makeSLURMs.py 1 busco.sh
ln -s ../../32_genePredictionComp/UnmaskedBraker/braker/genome738sl.polished.mitoFixed/augustus.gff3 genome738sl.polished.mitoFixed_ALLDATA.gff3


This is the actual busco run that ran. There were lots of things that changed when switching to condo2017. The python has to be pointing directly to python3, the path to busco and nematoda_odb9 changed.


/opt/rit/app/python/3.6.0/bin/python $BUSCO_HOME/BUSCO.py -f -i B.xylophilus.pep.fasta -l /work/GIF/software/programs/busco/2.0/nematoda_odb9/ -o B.xylophilus.pep.fasta.out -t B.xylophilus.pep.fasta.tmp -c 16 -m prot
/opt/rit/app/python/3.6.0/bin/python $BUSCO_HOME/BUSCO.py -f -i G.ellingtonae.pep.fasta -l /work/GIF/software/programs/busco/2.0/nematoda_odb9/ -o G.ellingtonae.pep.fasta.out -t G.ellingtonae.pep.fasta.tmp -c 16 -m prot
/opt/rit/app/python/3.6.0/bin/python $BUSCO_HOME/BUSCO.py -f -i G.pallida.pep.fasta -l /work/GIF/software/programs/busco/2.0/nematoda_odb9/ -o G.pallida.pep.fasta.out -t G.pallida.pep.fasta.tmp -c 16 -m prot
/opt/rit/app/python/3.6.0/bin/python $BUSCO_HOME/BUSCO.py -f -i G.rostochiensis.pep.fasta -l /work/GIF/software/programs/busco/2.0/nematoda_odb9/ -o G.rostochiensis.pep.fasta.out -t G.rostochiensis.pep.fasta.tmp -c 16 -m prot
/opt/rit/app/python/3.6.0/bin/python $BUSCO_HOME/BUSCO.py -f -i H.glycines.pep.fasta -l /work/GIF/software/programs/busco/2.0/nematoda_odb9/ -o H.glycines.pep.fasta.out -t H.glycines.pep.fasta.tmp -c 16 -m prot
/opt/rit/app/python/3.6.0/bin/python $BUSCO_HOME/BUSCO.py -f -i M.hapla.pep.fasta -l /work/GIF/software/programs/busco/2.0/nematoda_odb9/ -o M.hapla.pep.fasta.out -t M.hapla.pep.fasta.tmp -c 16 -m prot
/opt/rit/app/python/3.6.0/bin/python $BUSCO_HOME/BUSCO.py -f -i M.incognita.pep.fasta -l /work/GIF/software/programs/busco/2.0/nematoda_odb9/ -o M.incognita.pep.fasta.out -t M.incognita.pep.fasta.tmp -c 16 -m prot


BUSCO2.0    nematoda_odb9 B.xylophilus      C:79.7%[S:77.6%,D:2.1%],F:6.0%,M:14.3%,n:982            783 Complete BUSCOs (C)     762 Complete and single-copy BUSCOs (S)     21  Complete and duplicated BUSCOs (D)      59  Fragmented BUSCOs (F)       140 Missing BUSCOs (M)      982 Total BUSCO groups searched                     
BUSCO2.0    nematoda_odb9 H.glycines        C:71.9%[S:50.3%,D:21.6%],F:7.9%,M:20.2%,n:982           706 Complete BUSCOs (C)     494 Complete and single-copy BUSCOs (S)     212 Complete and duplicated BUSCOs (D)      78  Fragmented BUSCOs (F)       198 Missing BUSCOs (M)      982 Total BUSCO groups searched                     
BUSCO2.0    nematoda_odb9 G.ellingtonae     C:70.7%[S:61.5%,D:9.2%],F:10.2%,M:19.1%,n:982           694 Complete BUSCOs (C)     604 Complete and single-copy BUSCOs (S)     90  Complete and duplicated BUSCOs (D)      100 Fragmented BUSCOs (F)       188 Missing BUSCOs (M)      982 Total BUSCO groups searched                     
BUSCO2.0    nematoda_odb9 G.rostochiensis   C:70.7%[S:68.2%,D:2.5%],F:9.4%,M:19.9%,n:982            695 Complete BUSCOs (C)     670 Complete and single-copy BUSCOs (S)     25  Complete and duplicated BUSCOs (D)      92  Fragmented BUSCOs (F)       195 Missing BUSCOs (M)      982 Total BUSCO groups searched                     
BUSCO2.0    nematoda_odb9 M.hapla           C:59.0%[S:57.0%,D:2.0%],F:9.7%,M:31.3%,n:982            580 Complete BUSCOs (C)     560 Complete and single-copy BUSCOs (S)     20  Complete and duplicated BUSCOs (D)      95  Fragmented BUSCOs (F)       307 Missing BUSCOs (M)      982 Total BUSCO groups searched                     
BUSCO2.0    nematoda_odb9 M.incognita       C:50.8%[S:34.1%,D:16.7%],F:8.2%,M:41.0%,n:982           499 Complete BUSCOs (C)     335 Complete and single-copy BUSCOs (S)     164 Complete and duplicated BUSCOs (D)      81  Fragmented BUSCOs (F)       402 Missing BUSCOs (M)      982 Total BUSCO groups searched                     
BUSCO2.0    nematoda_odb9 G.pallida         C:50.7%[S:46.6%,D:4.1%],F:9.4%,M:39.9%,n:982            498 Complete BUSCOs (C)     458 Complete and single-copy BUSCOs (S)     40  Complete and duplicated BUSCOs (D)      92  Fragmented BUSCOs (F)       392 Missing BUSCOs (M)      982 Total BUSCO groups searched                         
```
### post analysis
```
Extract all of the complete busco genes


####cat run*/full* |grep -v "#" |awk '$2=="Complete"' |awk '{print $1}' |sort |uniq -c|sort  -k1,1nr |awk '$1==7 {print $2}' >single_copy_busco_IDs
####cat single_copy_busco_IDs |while read line; do echo "grep \""$line"\"  run*/full*  >"$line".list"; done >extractSingleCopy.sh


cdbfasta B.xylophilus.pep.fasta
cdbfasta G.pallida.pep.fasta
cdbfasta G.ellingtonae.pep.fasta
cdbfasta G.rostochiensis.pep.fasta
cdbfasta H.glycines.pep.fasta
cdbfasta M.hapla.pep.fasta
cdbfasta M.incognita.pep.fasta

####ls -1 *list |while read line; do echo "awk '{print \$3}' "$line" |cdbyank \$1 >>"$line".fasta"; done  >fasta.extraction.sh
####sh fasta.extraction.sh B.xylophilus.pep.fasta.cidx
####sh fasta.extraction.sh G.ellingtonae.pep.fasta.cidx
####sh fasta.extraction.sh G.pallida.pep.fasta.cidx
####sh fasta.extraction.sh G.rostochiensis.pep.fasta.cidx
####sh H.glycines.pep.fasta.cidx
####sh fasta.extraction.sh H.glycines.pep.fasta.cidx
####sh fasta.extraction.sh M.hapla.pep.fasta.cidx
####sh fasta.extraction.sh M.incognita.pep.fasta.cidx


Number of busco genes found in all 7 nematodes
ls -1 *list |wc
67

#get genes that are found in more than 3 species
cat run*/full* |grep -v "#" |awk '$2=="Complete"' |awk '{print $1}' |sort |uniq -c|sort  -k1,1nr |awk '$1>3 {print $2}' >single_copy_busco_IDs
cat single_copy_busco_IDs |while read line; do echo "grep \""$line"\"  run*/full*  >"$line".list"; done >extractSingleCopy.sh
sh extractSingleCopy.sh
#how many genes are found in more than 3 species.  This should be more powerful for the tree.
less extractSingleCopy.sh |wc
   651    2604   32550

ls -1 *list |while read line; do echo "awk '{print \$3}' "$line" |cdbyank \$1 >>"$line".fasta"; done  >fasta.extraction.sh
sh fasta.extraction.sh B.xylophilus.pep.fasta.cidx
sh fasta.extraction.sh G.ellingtonae.pep.fasta.cidx
sh fasta.extraction.sh G.pallida.pep.fasta.cidx
sh fasta.extraction.sh G.rostochiensis.pep.fasta.cidx
sh H.glycines.pep.fasta.cidx
sh fasta.extraction.sh H.glycines.pep.fasta.cidx
sh fasta.extraction.sh M.hapla.pep.fasta.cidx
sh fasta.extraction.sh M.incognita.pep.fasta.cidx
```
### guidance
```
All fasta files have been generated, now time for alignment.


#/work/GIF/remkv6/Baum/CamTechGenomeComparison/47_busco_nematode/Busco_Prot
sed -i 's/*//g' *.list.fasta
for file in *list.fasta; do echo "perl  /work/GIF/software/programs/guidence/2.02/www/Guidance/guidance.pl --seqFile "$file" --msaProgram PRANK --seqType aa --proc_num 16 --outDir guidance_"$file"_test"; done >guidance.cmds
#remove aruns email
#module load  GIF/guidence
python makeSLURMp.py 100 guidance.cmds
for f in *sub; do qsub $f; done

#Get the alignment files prepped for raxml


cp ~/common_scripts/runRAxML.sh
#modified to make it work, cut had to be set at -2. and fixed some bugs. pushed, so should be current.

#had to rename all of the prank alignment files so that the names were unique
for subdir in guidance_E*; do mv $subdir/MSA.PRANK.aln $subdir.MSA.PRANK.aln; done;
#generate commands
for f in  *test.MSA.PRANK.aln; do echo "sh ./runRAxML.sh $f AA";done >raxml.cmds
#split
python makeSLURMs.py 93 raxml.cmds
#submit
for f in *sub; do qsub $f; done
4
```
### Nematoda protein alignment using exonerate
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/71_BuscoExonerate
wget http://busco.ezlab.org/v2/datasets/nematoda_odb9.tar.gz
ln -s ../58_Renamatorium/1_genomeNgff/genome738sl.polished.mitoFixed.fa
awk 'NR<17 {print NR}' exonerate.out |awk -v interation=1 '{print "exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid "$1" --querychunktotal 16 --ryo \">%ti (%tab - %tae)\\n%tas\\n\" >>exonerate.out &"}' >exonerate.sh

exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 1 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 2 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 3 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 4 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 5 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 6 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 7 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 8 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 9 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 10 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 11 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 12 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 13 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 14 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 15 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
exonerate --model protein2genome nematoda_odb9/ancestral genome738sl.polished.mitoFixed.fa --querychunkid 16 --querychunktotal 16 --ryo ">%ti (%tab - %tae)\n%tas\n" >>exonerate.out &
```
### Duplicated BUSCO gene overlaps
```

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/47_busco_nematode/Busco_Prot/run_H.glycines.pep.fasta.out

#how many duplicated busco groups are there, now that we know busco was calling isoforms of the same gene, as duplicates?
less full_table_H.glycines.pep.fasta.out.tsv |awk '$2=="Duplicated"' |sed 's/\./\t/1' |grep -w "t1"  - |awk '{print $1}' |awk '{print $1}' |sort|uniq -c |awk '$1>1' |wc
   158     316    3160
#what is the new duplicated percentage?
158/982=16.1%
#this used to be 22% and 212 genes

#How many genes represent these busco gene groups?
less full_table_H.glycines.pep.fasta.out.tsv |awk '$2=="Duplicated"' |sed 's/\./\t/1' |grep -w "t1"  - |awk '{print $1}' |awk '{print $1}' |sort|uniq -c |awk '$1>1{print $2}' |grep -w -f - <(less full_table_H.glycines.pep.fasta.out.tsv |awk '$2=="Duplicated"') |awk '{print $3}' |sed 's/\./\t/1' |cut -f 1 |sort|uniq|wc
   343     343    2295

#How many busco gene groups get added to the total?
212-158=54 genes to add.
#old single copy count was 494
494+54=548 busco gene groups
#what is the new percentage
548/982 = 55.8%

less full_table_H.glycines.pep.fasta.out.tsv |awk '$2=="Duplicated"' |sed 's/\./\t/1' |grep -w "t1"  - |awk '{print $1}' |awk '{print $1}' |sort|uniq -c |awk '$1>1{print $2}' |grep -w -f - <(less full_table_H.glycines.pep.fasta.out.tsv |awk '$2=="Duplicated"') |awk '{print $3}' |sed 's/\./\t/1' |cut -f 1 |sort|uniq >DuplicatedBuscosGene.list



#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/26_ExpressionSets
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/47_busco_nematode/Busco_Prot/run_H.glycines.pep.fasta.out/DuplicatedBuscosGene.list .
cp ../1_genomeNgff/geneRenamer.sh .
sed -i 's/augustus.aa/DuplicatedBuscosGene.list/g' geneRenamer.sh


cat <(cat DuplicatedBuscosGene.list SupportedIRFMergeClassified.list |sort|uniq -c |awk '$1==2{print $2}') <(cat DuplicatedBuscosGene.list  NonRedundantLtrRetroelementMerge.list|sort|uniq -c |awk '$1==2{print $2}') <(cat DuplicatedBuscosGene.list ../27_TandemRedo/TotalSectionsGene.list |sort|uniq -c |awk '$1==2{print $2}') |sort|uniq|wc
   119     119    2142
```
### Removing wrongly duplicated buscos for G. ellingtonae
```

How many busco gene groups are there that are duplicated?
less full_table_G.ellingtonae.pep.fasta.out.tsv |awk '$2=="Duplicated"' |sed 's/\./\t/1' |grep -w "t1"  - |awk '{print $1}' |awk '{print $1}' |sort|uniq -c |awk '$1>1' |wc
    35      70     700
What is the duplicated percentage?
35/982 = 3.6%
How many are considered complete now?
#old count was 604 single and 90 duplicated
Now 604  + 55 = 659
659/982 = 67.1%
```
