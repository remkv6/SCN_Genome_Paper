# Need some functional annotation for the SCN genes
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/32_genePredictionComp/UnmaskedBraker/braker/genome738sl.polished.mitoFixed/InterproScan
Need to annotate the gene models with interproscan.
Arun suggested that I use blast2go rather than the command line version of interproscan

ln -s ../augustus.aa unmasked_738mitofixed.aa
fasta-splitter.pl --n-parts 4 unmasked_738mitofixed.aa
submitted the first two file-parts myself, arun said that he would get 3 and 4.  
A blast to the swissprot and uniprot are recommended, and additional GO analyses will be needed.

#
mkdir swissprot
mkdir uniprot
cd swissprot/
fasta-splitter.pl --n-parts 100 ../unmasked_738mitofixed.aa
for file in unmasked_738mitofixed.part-*; do   num=$(echo $file | cut -d "-" -f 2);   echo "blastp -query ${file} -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -max_target_seqs 1 -outfmt 6 > td_uniref.blastp_${num}.outfmt6" >> cmds.txt; done
cp ~/common_scripts/makeSLURMp.py .
modify makeSLURMp.py to load blast and parallel
for f in *.sub; do qsub $f; done

#uniprot
cd ../uniprot/
fasta-splitter.pl --n-parts 100 ../unmasked_738mitofixed.aa
for file in longest_orfs.cds.part*; do num=$(echo $file | cut -d "-" -f 2);   echo "blastx -query ${file} -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -max_target_seqs 1 -outfmt 6 > td_swissprot.blastx_${num}.outfmt6" >> cmds.txt; done
cp ../swissprot/makeSLURMp.py .
python makeSLURMp.py 16 cmds.txt
for f in *.sub; do qsub $f; done
```
```
Split Merged gene models into 4 and ran with blast2go using interproscan with all
Files have been transferred and unzipped here:
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/Interproscan fasta-splitter.pl –n-parts ../all.augustus.aa
```

### Blast analyses for swissprot
```
fasta-splitter.pl --n-parts 16 ../all.augustus.aa
for f in *part*; do echo "blastp -query "$f" -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > "$f".xml"; done >cmds.txt

#Created pbs script from above lines and ran in parallel.
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 96:00:00
#SBATCH -J cmds_0
#SBATCH -o cmds_0.o%j
#SBATCH -e cmds_0.e%j
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module use /work/GIF/software/modules
module load ncbi-blast
module load parallel
parallel --joblog cmds_progress_0.log --workdir $PWD <<FIL
blastp -query all.augustus.aa.part-01 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-01.xml
blastp -query all.augustus.aa.part-01.tab -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-01.tab.xml
blastp -query all.augustus.aa.part-02 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-02.xml
blastp -query all.augustus.aa.part-03 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-03.xml
blastp -query all.augustus.aa.part-04 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-04.xml
blastp -query all.augustus.aa.part-05 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-05.xml
blastp -query all.augustus.aa.part-06 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-06.xml
blastp -query all.augustus.aa.part-07 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-07.xml
blastp -query all.augustus.aa.part-08 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-08.xml
blastp -query all.augustus.aa.part-09 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-09.xml
blastp -query all.augustus.aa.part-10 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-10.xml
blastp -query all.augustus.aa.part-11 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-11.xml
blastp -query all.augustus.aa.part-12 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-12.xml
blastp -query all.augustus.aa.part-13 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-13.xml
blastp -query all.augustus.aa.part-14 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-14.xml
blastp -query all.augustus.aa.part-15 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-15.xml
blastp -query all.augustus.aa.part-16 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-16.xml
FIL
qstat -f ${PBS_JOBID} |head
```

### Blast to uniprot
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/Interproscan/uniprot


fasta-splitter.pl --n-parts 16 ../all.augustus.aa

for f in *part*; do echo "blastp -query "$f" -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > "$f".xml"; done >cmds.txt
#used this pbs script
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 96:00:00
#SBATCH -J cmds_0
#SBATCH -o cmds_0.o%j
#SBATCH -e cmds_0.e%j
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module use /work/GIF/software/modules
module load ncbi-blast
module load parallel
parallel --joblog cmds_progress_0.log --workdir $PWD <<FIL
blastp -query all.augustus.aa.part-01 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-01.xml
blastp -query all.augustus.aa.part-02 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-02.xml
blastp -query all.augustus.aa.part-03 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-03.xml
blastp -query all.augustus.aa.part-04 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-04.xml
blastp -query all.augustus.aa.part-05 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-05.xml
blastp -query all.augustus.aa.part-06 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-06.xml
blastp -query all.augustus.aa.part-07 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-07.xml
blastp -query all.augustus.aa.part-08 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-08.xml
blastp -query all.augustus.aa.part-09 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-09.xml
blastp -query all.augustus.aa.part-10 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-10.xml
blastp -query all.augustus.aa.part-11 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-11.xml
blastp -query all.augustus.aa.part-12 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-12.xml
blastp -query all.augustus.aa.part-13 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-13.xml
blastp -query all.augustus.aa.part-14 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-14.xml
blastp -query all.augustus.aa.part-15 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-15.xml
blastp -query all.augustus.aa.part-16 -db /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90 -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -outfmt 5 > all.augustus.aa.part-16.xml
FIL
qstat -f ${PBS_JOBID} |head
```

#### Convert BLAST xml to tab
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/Interproscan/swiss-prot
for f in *xml; do python ~/common_scripts/blastXML2Tab.py -o ${f%.*}.tab -c std $f; done
cat *.tab » augustus_pep_swissprot_blast.txt

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/Interproscan/uniprot
for f in *xml; do python ~/common_scripts/blastXML2Tab.py -o ${f%.*}.tab -c std $f; done
cat *tab >>augustus_pep_uniprot_blast.txt
```
#### Convert interproscan, uniref blast, and swissprot blast to gff3
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/52_functional


for f in *.xml; do awk 'NR>2' "$f" |head -n -1 >>all.interproscan.txt; done
cat <(awk 'NR<3' g15135.t1.xml) <(less all.interproscan.txt ) <(tail -n 1 g15135.t1.xml) >interproscan.combined.xml
#this did not work on the head node even after loading java
interproscan.sh --mode convert -f tsv,gff3 -i interproscan.combined.xml -b scnInterproscanResultsConverted

#how many genes were annotated with interproscan
less scnInterproscanResultsConverted.tsv |cut -f 1 |sort|uniq|wc
  22376   22376  215798

#because the maker-p had problems with the extra annotations, I had to join all of the annotations by hand.  This converts the multi-lined annotations for each gene into single lines for join.
less scnInterproscanResultsConverted.tsv |awk '{print $1}'|sort|uniq|while read line; do grep -w $line scnInterproscanResultsConverted.tsv|tr "\n" "\t" |awk '{print $line, $0}' >>scnInterproscanResultsConverted.oneLine.tsv;done &

#convert swissprot blast xml output to tab
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/Interproscan/swiss-prot
for f in *xml; do python ~/common_scripts/blastXML2Tab.py -o ${f%.*}.tab -c std $f; done
cat *tab >>augustus_pep_swissprot_blast.txt

#convert uniref blast xmls to tab
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/Interproscan/uniprot
for f in *xml; do python ~/common_scripts/blastXML2Tab.py -o ${f%.*}.tab -c std $f; done
cat *tab >>augustus_pep_uniprot_blast.txt

#convert tab to gff3 for swissprot
maker_functional_gff /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/Interproscan/swissprot/augustus_pep_swissprot_blast.txt ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 >swissprot.reformatted.function.gff3


#There was a problem with the maker script below.  It would not accept the uniref headers, so I had to convert them as below.
sed '/>/ s/Tax=/OS=/g' /work/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/uniref90/uniref90.fasta |sed '/>/s/RepID=/GN=/g' |sed '/>/s/$/ PE=/g' >modified.uniref90.fasta
#convert tab to gff3 for uniref
maker_functional_gff modMakerP/modified.uniref90.fasta ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/Interproscan/uniprot/augustus_pep_uniprot_blast.txt ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 >uniref.reformatted.function.gff3

#I lost the history to this, but I essentially used join on the first 8columns of the two gff files above by changing tabs to "#".
#That is how I generated the combined.uniref.swissprot.gff from unirefprot.reformatted.function.4join and swissprot.reformatted.function.4join.
#I then used genome tools -tidy to fix the errors that this generated in the gff.  

#this just extracts the annotation information from the gff above, and adds a column in column 1 for join -- this is a partial of the notes that are missing
awk '$3=="mRNA" {print $0,$9}' combined.uniref.swissprot.gff |sed -e 's/;/\t/1' -e 's/ID=//1' |awk '{print $9,$0}' |sed -e 's/\t/\tID=/8' -e 's/\t/;/9' >combined.uniref.swissprot.4join.gff

# this is essentially the gff in the wrong format.  
join  -a 1 -1 1 -2 1 <(awk '$3=="mRNA" {print $0,$9}' combined.uniref.swissprot.gff |sed -e 's/;/\t/1' -e 's/ID=//1' |awk '{print $9,$0}' |sed -e 's/\t/\tID=/8' -e 's/\t/;/9' |sort -k1,1V) <(less scnInterproscanResultsConverted.oneLine.tsv |sort -k1,1V) |uniq  >augustusMergedPutativeFuntions.gff3

#this changes the above file to tabular, removes the first column gene index, concatenates with the original augustus.gff3 and sorts.  --There will still be duplicate entries here.  
less -S augustusMergedPutativeFuntions.gff3 |sed 's/ /\t/1' | sed 's/ /\t/1' |sed 's/ /\t/1' |sed 's/ /\t/1' |sed 's/ /\t/1' |sed 's/ /\t/1' |sed 's/ /\t/1' |sed 's/ /\t/1' |sed 's/ /\t/1' |cut -f 2- |cat - ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 |sort -k1,1 -k4,4n > test.gff
#This cleaned up all of the duplicates I think
gt gff3 -tidy test.gff
mv test.gff augustusFunctionalAnnotation.gff3


################################
################################
I had trouble getting all three gffs merged, so arun did that for me.  Here is the functional annotation.
/work/GIF/arnstrm/Baum/20151221_Baum_Hg_annotation/annotation_merging_20170419/augustus_swissprot_iprscan.gff3
```
