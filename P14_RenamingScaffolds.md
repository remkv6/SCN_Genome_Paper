###  Need to rename scaffolds so that they are all similarly named
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium
Renaming scaffold names

#renames scaffolds in the gff
cp ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 .
bioawk -c fastx '{print $name,length($seq)}' genome738.mitofixed.fa |sort -k2,2nr|awk '{print "sed s/"$1"/scaffold_"NR"/g"}' |sed "s/sed  s/sed 's/g" |sed "s/$/' fixed.augustus.gff3/g" |sed "s/sed s/sed -i 's/g" |sort -k2,2Vr >gff.renamer.sh
sh gff.renamer.sh &

Renaming gene names in gff\\
less augustus.gff3 |sed 's/\.exon.*;/;/g'|sed 's/\.CDS.*;/;/g' |sort -k1,1 -k4,5n |sed 's/Name=;//g' >fixed.augustus.gff3
maker_map_ids --prefix Hetgly. --justify 9 fixed.augustus.gff3 --abrv_tran T --abrv_gene G --iterate 1 --suffix '.' >renaming.index
map_gff_ids renaming.index fixed.augustus.gff3
#there were some errors based on line 117 of map_gff_ids.  All genes appear to be changed.  The alias is not getting placed in the col9
#he instructions to this need to be revised to show it wants `Parent=gene;ID=gene;Alias=gene` in column 9 of the gff

#renames the fasta file scaffold names
sed  's/augustus.gff3/genome738sl.polished.mitoFixed.fa/g' gff.renamer.sh >scaffold.renamer.sh
sed -i 's/\^/\^>/g' scaffold.renamer.sh
sh scaffold.renamer.sh &
```

#### RepeatMasker/repeatmodeler
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/repeats

mkdir repeats
cd repeats/
cp ../../22_RepeatModeler/genome738sl.polished.mitoFixed.fa.out.gff .
cp ../1_genomeNgff/gff.renamer.sh .
sed -i 's/augustus.gff3/genome738sl.polished.mitoFixed.fa.out.gff/g' gff.renamer.sh
sh gff.renamer.sh &
```
#### Renaming scaffolds from synteny analysis
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/3_synteny
I think three files need to be changed in each of these. syntenic ribbons.txt, karyotype.both.txt, and circos.test.delete


#
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/3_synteny
cp -rf ../../26_GloboderaSynteny/circos pallida
cp -rf ../../39_Meloidogyne/hapla/circos/ hapla
cp -rf ../../39_Meloidogyne/incognita/circos/ incognita/

cp ../1_genomeNgff/gff.renamer.sh karyotype.renamer.sh

sed  's|*/karyotype.both.txt|*/syntenic.ribbons.txt|g' karyotype.renamer1.sh >ribbon.renamer.sh
sh ribbon.renamer.sh &
sed  's|*/karyotype.both.txt|*/circos.test.delete|g' karyotype.renamer1.sh >conf.renamer.sh
sh conf.renamer.sh &
sed  's|augustus.gff3|*/karyotype.both.txt|g' karyotype.renamer.sh |sed 's|/|\t|2' |awk '{print $1,$2,$4,$5}' |sed "s/ / \'/2" >karyotype.renamer1.sh

ELLINGTONAE
cp -rf ../../26_GloboderaSynteny/otherGloboderaGenomes/ellingtonae/braker/GCA_001723225.1_ASM172322v1_genomic/circos/ ellingtonae


#karyotype.both.txt
cp  ../karyotype.renamer1.sh .
#looks like this karyotype file was a simple addition.
sed -i 's|*/karyotype.both.txt| H.glycines.kary|g' karyotype.renamer1.sh
sh karyotype.renamer1.sh
cat H.glycines.kary subj.kary >karyotype.both.txt


#syntenic.ribbons.txt
cut -d " " -f 4- syntenic.ribbons.txt >elling.syntenicribbons.4paste
cut -d " " -f -3 syntenic.ribbons.txt >hgly.syntenicribbons.4paste
cp ../karyotype.renamer1.sh .
sed -i 's|*/karyotype.both.txt|hgly.syntenicribbons.4paste|g' karyotype.renamer1.sh
sh karyotype.renamer1.sh
paste hgly.syntenicribbons.4paste elling.syntenicribbons.4paste |tr "\t" " " >syntenic.ribbons.txt


#circos.test.delete
#COPY THE ellingtonae scaffolds
sed  's/hgly.syntenicribbons.4paste/circos.test.delete/g' karyotype.renamer1.sh >circosconf.renamer.sh
sh circosconf.renamer.sh
#paste scaffolds into circos.test.delete


ROSTOCHIENSIS
cp -rf ../../26_GloboderaSynteny/otherGloboderaGenomes/rostochiensis/circos/ rostochiensis

#karyotype.both.txt
cp  ../karyotype.renamer1.sh .
sed -i 's|*/karyotype.both.txt| H.glycines.kary|g' karyotype.renamer1.sh
sh karyotype.renamer1.sh
cat H.glycines.kary subj.kary >karyotype.both.txt

#syntenic.ribbons.txt
cut -d " " -f 4- syntenic.ribbons.txt >rostoch.syntenicribbons.4paste
cut -d " " -f -3 syntenic.ribbons.txt >hgly.syntenicribbons.4paste
cp ../karyotype.renamer1.sh .
sed -i 's|*/karyotype.both.txt|hgly.syntenicribbons.4paste|g' karyotype.renamer1.sh
sh karyotype.renamer1.sh
paste hgly.syntenicribbons.4paste rostoch.syntenicribbons.4paste |tr "\t" " " >syntenic.ribbons.txt

#circos.test.delete
#COPY THE ellingtonae scaffolds
sed  's/hgly.syntenicribbons.4paste/circos.test.delete/g' karyotype.renamer1.sh >circosconf.renamer.sh
sh circosconf.renamer.sh
#paste scaffolds into circos.test.delete


##final change to all circos.test.delete
sed -i 's|/shared/software/GIF/programs/circos/0.69.2/|/work/GIF/software/programs/circos/0.69-4/|g' */circos.test.delete
```

#### effector mapping
```
mkdir 4_effectors/
cd 4_effectors/
cp ../../29_effectorMapping/effector.gmapped.filtered.gff3 .
cp ../1_genomeNgff/gff.renamer.sh .
sed -i 's/augustus.gff3/effector.gmapped.filtered.gff3/g' gff.renamer.sh
sh gff.renamer.sh &
```
#### transcript mapping
```
mkdir 5_transcripts
cd 5_transcripts
cp ../../30_transcriptMappingUpdated/nematode_transcripts_Trinity.fasta_sorted.gff .
cp ../1_genomeNgff/gff.renamer.sh .
sed -i 's/augustus.gff3/nematode_transcripts_Trinity.fasta_sorted.gff/g' gff.renamer.sh
sh gff.renamer.sh &
```
#### isoseq mapping

```
mkdir 6_isoseq
cd 6_isoseq/
cp ../../34_IsoSeq/isoseq_738_sorted.gff .
cp ../1_genomeNgff/gff.renamer.sh .
sed -i 's/augustus.gff3/isoseq_738_sorted.gff/g' gff.renamer.sh
sh gff.renamer.sh &
```
#### ccs mapping
```
mkdir 7_ccsReads
cd 7_ccsReads/
cp ../../36_deduplicate122016/GMAP/ccsReads.nocol9_sorted_sorted.gff .
cp ../1_genomeNgff/gff.renamer.sh .
sed -i 's/augustus.gff3/ccsReads.nocol9_sorted_sorted.gff/g' gff.renamer.sh
sh gff.renamer.sh &
```

#### NCBI nucleotide database H. glycines downloads, mapping
```
mkdir 8_Nucleotide
cd 8_Nucleotide/
cp ../../37_Align_NCBI_Nucleotide/Nucleotide_mapped_sorted.gff .
cp ../1_genomeNgff/gff.renamer.sh .
sed -i 's/augustus.gff3/Nucleotide_mapped_sorted.gff/g' gff.renamer.sh
sh gff.renamer.sh &
```

#### NCBI EST mapping
```
mkdir 9_EST
cd 9_EST/
cp ../../38_Align_NCBI_EST/ESTs_mapped_sorted_gff_sorted.gff .
cp ../1_genomeNgff/gff.renamer.sh .
sed -i 's/augustus.gff3/ESTs_mapped_sorted_gff_sorted.gff/g' gff.renamer.sh
sh gff.renamer.sh &
```

#### tandem duplication
```
mkdir 10_tandemDups
cd 10_tandemDups/
cp ../../40_tandemdups/split/sorted.final.renamed.out .
sed -i 's/>//g'  sorted.final.renamed.out
cp ../1_genomeNgff/gff.renamer.sh .
sed -i 's/augustus.gff3/sorted.final.renamed.out/g' gff.renamer.sh
sh gff.renamer.sh &

#rename the 1copy tandem genome
cp ../../40_tandemdups/split/1tandemcopyonly.masked.genome738 .
cp ../1_genomeNgff/scaffold.renamer.sh .
sed -i 's/genome738sl.polished.mitoFixed.fa/1tandemcopyonly.masked.genome738/g' scaffold.renamer.sh
sh scaffold.renamer.sh &

#rename the 1copy tandem gff
cp ../../40_tandemdups/split/1tandemcopyonly.masked.genome738.gff3 .
cp ../1_genomeNgff/gff.renamer.sh gff.renamer1.sh
sed -i 's/augustus.gff3/1tandemcopyonly.masked.genome738.gff3/g' gff.renamer1.sh
sh gff.renamer1.sh &
```

#### Syntenic tracks
```
mkdir 11_syntenicTracks
cd 11_syntenicTracks/
cp ../../46_SCN_syntenic_tracks/*gff .
cp ../1_genomeNgff/gff.renamer.sh .
sed 's/augustus.gff3/*gff/g' gff.renamer1.sh |less
cp ../1_genomeNgff/gff.renamer.sh .
sed -i 's/augustus.gff3/*gff/g' gff.renamer.sh
sh gff.renamer.sh &
```
#### Functional gff
```
mkdir 12_functional
cd 12_functional
cp /work/GIF/arnstrm/Baum/20151221_Baum_Hg_annotation/annotation_merging_20170419/augustus_swissprot_iprscan.gff3 .
cp ../1_genomeNgff/gff.renamer.sh .
sed -i 's/augustus.gff3/augustus_swissprot_iprscan.gff3/g' gff.renamer.sh
sh gff.renamer.sh &
```
#### helitrons gff
```
mkdir 13_helitron
cd 13_helitron/
ls -ltrh  ../../54_helitronScanner/
cp  ../../54_helitronScanner/helitrons.gff .
cp ../1_genomeNgff/gff.renamer.sh .
sed -i 's/augustus.gff3/helitrons.gff/g' gff.renamer.sh
sh gff.renamer.sh &
```
