# Checking to see if internal synteny still exists after scaffold merging based on synteny

### Synteny rerun
```
Arun has already ran these genomes through braker to call genes using the same genemark training.

ln -s /data021/GIF/arnstrm/Baum/GenePrediction_HgCamtech_20160809/G368_Sp1/autoAugPred_abinitio/predictions/G368_Hg_augustus.gff3
ln -s /data021/GIF/arnstrm/Baum/GenePrediction_HgCamtech_20160809/G368_Sp1/autoAugPred_abinitio/predictions/G368_Hg_augustus.aa
ln -s /data021/GIF/arnstrm/Baum/GenePrediction_HgCamtech_20160809/G2692_nSp1/autoAugPred_abinitio/predictions/G2692_v2_Hg_augustus.aa
ln -s /data021/GIF/arnstrm/Baum/GenePrediction_HgCamtech_20160809/G2692_nSp1/autoAugPred_abinitio/predictions/G2692_v2_Hg_augustus.gff3

sed -e 's/ID=//g' G368_Hg_augustus.gff3 -e 's/;/\t/g' -e 's/NODE_//g'|grep gene -|awk '{print $3 "\t" $9 "\t" substr($1, 0, 8) "\t" $4 "\t" $5 "\t" $7 "\t" "t"}' |sed 's/scaffold_//g' |sort -k 5 -V|awk '{print ">"$2 "\t" "1" "\t" $4 "\t" $5 "\t" $6 "\t" $3}'|sort -k 1 -V |sed -e 's/F+//g' -e 's/F-//g' -e 's/F//g' >augustus.368.gff3.reformat1
tr "\n" "\t" < G368_Hg_augustus.aa |sed 's/>/\n>/g'|sed 's/\t//g'|grep ".t1"|sed  's/.t1/\t/g'>3684join2
join -1 1 -2 1 augustus.368.gff3.reformat1 3684join2 >joined.386.3
sed 's/ /\n/6' joined.386.3 >operon.db

#generating for 2692 genome, genome.db
sed -e 's/ID=//g' G2692_v2_Hg_augustus.gff3 -e 's/;/\t/g' -e 's/NODE_//g'|grep gene -|awk '{print $3 "\t" $9 "\t" substr($1, 0, 8) "\t" $4 "\t" #$5 "\t" $7 "\t" "t"}' |sed 's/scaffold_//g' |sort -k 5 -V|awk '{print ">" $2 "\t" "1" "\t" $4 "\t" $5 "\t" $6 "\t" $3}'|sort -k 1 -V |sed -e 's/F+//g' -e 's/F-//g' -e 's/F//g' >augustus_G2692.gff3.reformat1
tr "\n" "\t" <G2692_v2_Hg_augustus.aa |sed 's/>/\n>/g'|sed 's/\t//g'|grep ".t1"|sed  's/.t1/\t/g'>26924join2
join -1 1 -2 1 augustus_G2692.gff3.reformat1 26924join2 >joined.2692.3
sed 's/ /\n/6' joined.2692.3 >genome.db
```

### iadhore setup

```
cp operon.db ../Nematode/
 cp genome.db ../Nematode/
 cp OPSCAN.o106726.condo iADHoRe_PLAZA.table
 cp iADHoRe_PLAZA.table ../Nematode/

#368 first, generating lst files and 368.ini
#/data021/GIF/remkv6/CamTechGenomeComparison/11finalsynth2692selfiadhore/Nematode/368
tr "\n" "\t" <../operon.db |sed 's/>/\n>/g'|awk '{print $1$5 "\t" $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls -1 *.lst|awk '$2=$1{print $1=substr($1,0,6), " 368/"$2}'>368.ini


##This was already generated
##moving to 2692/genome.db file reformat to .lst and .ini files
#/data021/GIF/remkv6/CamTechGenomeComparison/11finalsynth2692selfiadhore/Nematode/2692
tr "\n" "\t" < ../genome.db |sed 's/>/\n>/g'|awk '{print $1$5 "\t" $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls -1 *.lst|awk '$2=$1{print $1=substr($1,0,6), " 2692/"$2}'>2692.ini

#/data021/GIF/remkv6/CamTechGenomeComparison/06_iadhore_allscaf/Nematode
cat 2692/2692.ini 368/368.ini >Nematode.ini
#added to the Nematode.ini file
this goes at the top of the gene list for 2692 and a newline is added at the end of 2692's list before the genome= 368.
genome= 269

genome= 368

prob_cutoff=0.05
anchor_points=3
number_of_threads=16
visualizeAlignment=true
blast_table= iADHoRe_PLAZA.table
output_path= output
alignment_method=gg2
gap_size=15
cluster_gap=20
level_2_only=true
q_value=.9


#rearranging opscan output to form iadhore plaza table
#/data021/GIF/remkv6/CamTechGenomeComparison/06_iadhore_allscaf/opscan
##################################################################################
make sure to change the input file here next time you run

grep "CL"  OPSCAN.o104582.condo |awk '{print $3 "\t" $6}'  >iADHoRe_PLAZA.table
cp iADHoRe_PLAZA.table ../Nematode/
```

TSynteny still remains.

awk '$2=="269" {print $3"K"}' multiplicons.txt >STILLsyntenic
