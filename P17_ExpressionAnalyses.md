# a genome paper needs expression analyses, so here goes with Melissa's rnaseq

### Sebastians suggestions for draft paper
```
RNAseq:
Depends on what is going in the transcriptome paper. If it wont be its own section it should, as a minimum, be integrated it into all the sections. – are duplicated genes concordantly expressed, or have the diverged (can they be distinguished by only taking reads that map once at utr boundaries for example?)?. if they are diverged – how does this depend on direction (i.e. shared promoters)? – any evidence of HGT expression, transposons expression?

RE popseq:
Depends what the pops are and whether they correlate with pathotypes or not? As a minimum, are the effectors more diverse than other genes (snp density c.f. all others)? Re-tandem duplications, are they perfect copies or is one degenerating/evolving? Have TD genes accumulated more SNPs than non-TD (can we functionally link it to evolution?). same question but with transposon associated effectors (do we see any evidence that those associated with transposons have higher snp density than non-transposon associated effectors, or non-transposon associated genes?). Cross ref with orthologous gene cluster categories (are preferentially HG private genes more diverse in a population than those common to nematodes/common to CN)?
```

### Gather files and rename featurecounts file obtained from melissa's RNA-seq
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/20_Expression
ln -s ../17_EffectorProteinNetwork/Secretome.gene.list
ln -s ../17_EffectorProteinNetwork/Family1265.list
ln -s ../17_EffectorProteinNetwork/Family976.list
ln -s ../17_EffectorProteinNetwork/effector.pop.list
ln -s ../10_tandemDups/ontologenizer/tandem.gene.list

#have to update names in featurecounts file below to get some expression numbers
cp ../../45_RNA-seqJB/featureCounts/featureCounts_counts.txt .
cp ../1_genomeNgff/geneRenamer.sh .
cp ../1_genomeNgff/scaffold.renamer.sh .

sed -i 's/augustus.aa/featureCounts_counts.txt/g' geneRenamer.sh
sh geneRenamer.sh &

sed 's/genome738sl.polished.mitoFixed.fa/col2/g' scaffold.renamer.sh |sed 's/>//g' >scaffold.renamer1.sh
sh scaffold.renamer1.sh &
```
### Stringtie
```
Siva suggested stringtie because it could account for multiple mapped reads, as well as unique reads. Featurecount only does unique reads.


#/work/GIF/remkv6/Baum/CamTechGenomeComparison/45_RNA-seqJB
for f in *rnaseq_sorted.bam; do echo "stringtie "$f" -p 16 -G ../32_genePredictionComp/MergeIsoseqRnsaseq/braker/genome738sl.polished.mitoFixed1/augustus.gff3 -e -B -A "${f%.*}"_GeneAbundancepJ2_race3_Forrest > "${f%.*}"_GeneAbundance.gff"; done >stringtie.sh

cp ~/common_scripts/makeSLURMs.py .
vi makeSLURMs.py
python makeSLURMs.py 1 stringtie.sh
for f in stringtie*sub; do sbatch $f ; done
paste   <(sort -k1,1V pJ2_race3_Forrest_ACTTGA_L001_R1_001_fastq_gz_paired_fq_rnaseq_sorted_GeneAbundance.tab)  <(sort -k1,1V pJ2_s63_TGACCA_L001_R1_001_fastq_gz_paired_fq_rna
seq_sorted_GeneAbundance.tab|cut -f 7- )  <(sort -k1,1V pJ2_s63_TTAGGC_L001_R1_001_fastq_gz_paired_fq_rnaseq_sorted_GeneAbundance.tab|cut -f 7- )  <(sort -k1,1V ppJ2_PA3_ATCACG_L001_R1_001_fastq_gz_paired_fq
_rnaseq_sorted_GeneAbundance.tab|cut -f 7- )  <(sort -k1,1V ppJ2_PA3_TAGCTT_L001_R1_001_fastq_gz_paired_fq_rnaseq_sorted_GeneAbundance.tab|cut -f 7- ) >AllGeneAbundance.tab
#generates a tabular file with FPKM information averaged across the 5 rnaseq lines.
less -S AllGeneAbundance.tab |awk '{print $1,$2,$3,$4,$5,$6,($8+$11+$14+$17+$20)/5 }' >AllGeneAbundanceAverage.tab


#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/20_Expression
cp ../../45_RNA-seqJB/AllGeneAbundanceAverage.tab .
less geneRenamer.sh |sed 's|/|\t|g'|awk '{print $4,$5}' |sort -k1,1V |sed '1i\\' |paste - AllGeneAbundanceAverage.tab  |awk '{print $2,$4,$5,$6,$7,$8,$9}' >AllGeneAbundanceAverageMod.tab


#redid the effectors so that they represent the FIMO genes and the direct gmap aligments.
ln -s ../../../Baum/CamTechGenomeComparison/57_secretome/ontologenizer/population
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/10_tandemDups/ontologenizer
java -jar Ontologizer.jar -a simpleformat.ids -g go.obo -m Bonferroni -p population -s Transcripts2Genes.list
grep -v -w -f  <(sort population) <(cut -f 1 simpleformat.ids|sort|uniq) |sort|uniq >tandem.gene.list
wc -l ../tandem.gene.list
4774 ../tandem.gene.list
ln -s ../../../Baum/CamTechGenomeComparison/57_secretome/ontologenizer/population
```

### Re-evaluation of tandem orthologs
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/20_Expression/CallTandemOrthologs

#get the tandem gene fastas and self blast
awk '{print $1".t1"}' ../tandem.gene.list |cdbyank augustus.aa.cidx >tandem.gene.fasta
makeblastdb -in tandem.gene.fasta -dbtype prot -out tandem.gene.blastdb
blastp -db tandem.gene.blastdb -query tandem.gene.fasta -outfmt  '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -num_threads 15 -out tandem2tandemself.blast.out

#90% identity and 80% length filter to find orthologs.
paste <(awk '$13>(.8*$14) && $14>(.8*$13) &&$1!=$2 &&$3>90' tandem2tandemself.blast.out |awk '{print $1}' |sed 's/\.t1//g'  |while read line; do grep -w $line GenesExpressionSnpsLength.tab;done ) <(awk '$13>(.8*$14) && $14>(.8*$13) &&$1!=$2 &&$3>90' tandem2tandemself.blast.out |awk '{print $2}' |sed 's/\.t1//g' |while read line; do grep -w $line GenesExpressionSnpsLength.t
ab; done ) >TandemOrthologExpressionSnpsLength.tab


#Not in the best form yet
less TandemOrthologExpressionSnpsLength.tab |awk '$3==$13 && $7>5 && $17>5' |awk '$7<(.8*$17) || $17<(.8*$7) ' |less
 less TandemOrthologExpressionSnpsLength.tab |awk '$3==$13 && $7>5 && $17>5' |awk '$7<(.8*$17) || $17<(.8*$7) ' |sed 's/Hetgly.G/Hetgly.G\t/g' |awk '$13>$2' |sed 's/Hetgly.G\t/Hetgly.G/g' |awk '{print $7}' |summary.sh
Total:  2,817
Count:  58
Mean:   48
Median: 16
Min:    5
Max:    801
less TandemOrthologExpressionSnpsLength.tab |awk '$3==$13 && $7>5 && $17>5' |awk '$7<(.8*$17) || $17<(.8*$7) ' |sed 's/Hetgly.G/Hetgly.G\t/g' |awk '$13>$2' |sed 's/Hetgly.G\t/Hetgly.G/g' |awk '{print $17}' |summary.sh
Total:  2,891
Count:  58
Mean:   49
Median: 17
Min:    5
Max:    851

#LATER ADDITION, this associates the snpdensity and expression with the orthologs called above by blast -- needs viewed as network.
paste <(awk '{print $1}' tandemDuplicateOrtholog.network|while read line; do grep -w  $line GenesExpressionSnpsLength.tab;done) <(awk '{print $2}' tandemDuplicateOrtholog.network|while read line; do grep -w  $line GenesExpressionSnpsLength.tab;done) |awk '{print $1,$7,$10,$11,$17,$20}' >tandemDuplicateOrtholog.networkRedo
```

All else is contained within 05_ExpressionSNPDensityPlots.md
