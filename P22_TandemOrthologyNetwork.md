# Tandem duplication gene orthologs
```
The need here is to be able to make comparisons with the genes that are orthologous in tandemly duplicated regions. This paper is directly relevant to the comparisons that need to be made in the discussion.


#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/27_TandemRedo/01_CallTandemOrthologs

 #grab relevant files
ln -s ../TotalSectionsGene.list
ln -s ../../1_genomeNgff/augustus.aa.cidx
ln -s ../../1_genomeNgff/augustus.aa

#used cdbyank to get all of the proteins that are listit in the TotalSectionsGene.list above (all genes found in the tandem duplications)

blastp -db TandemGene.blastdb -query TandemGenes.fasta -outfmt  '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -num_threads 15 -out tandem2tandemself.blast.out
```
### Make the network
```

awk '$4>(.5*$13)  &&$1!=$2 &&$3>80' tandem2tandemself.blast.out |awk '{print $1,$2}' |sed 's/\./\t/2' |sed 's/\./\t/3' |awk '{print $1,$3}' |sort|uniq >Network3429

#attach snp density and expression in FPKM.
paste <(less Network3429 |awk '{print $1}' |while read line; do grep -w $line /work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/20_Expression/CallTandemOrthologs/GenesExpressionSnpsLength.tab;done |awk '{print $1,$7,$10}') <(less Network3429 |awk '{print $2}' |while read line; do grep -w $line /work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/20_Expression/CallTandemOrthologs/GenesExpressionSnpsLength.tab;done |awk '{print $1,$7,$10}' ) >Network3429TDOrthologs1

#attach the logfoldchange and the scaffold name
awk '$3=="gene"' ../../1_genomeNgff/fixed.augustus.gff3 |awk '{print $1,$9 }' |sed 's/ID=//g' |sed 's/;/\t/g' |awk '{print $1,$2}' |sort|uniq >ScaffoldIndex

paste <(less Network3429 |awk '{print $1}' |while read line; do grep -w $line /work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/20_Expression/CallTandemOrthologs/GenesExpressionSnpsLength.tab;done |awk '{print $1,$7,$10}') <(less Network3429 |awk '{print $1}' |while read line; do grep -w $line ../../26_ExpressionSets/diffexpr-resultsSortRenamed.tab;done |awk '{print $4}') <(less Network3429 |awk '{print $1}' |while read line; do grep -w $line ScaffoldIndex;done |awk '{print $1}' |uniq)  <(less Network3429 |awk '{print $2}' |while read line; do grep -w $line /work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/20_Expression/CallTandemOrthologs/GenesExpressionSnpsLength.tab;done |awk '{print $1,$7,$10}' )  <(less Network3429 |awk '{print $2}' |while read line; do grep -w $line ../../26_ExpressionSets/diffexpr-resultsSortRenamed.tab;done |awk '{print $4}') <(less Network3429 |awk '{print $2}' |while read line; do grep -w $line ScaffoldIndex;done |awk '{print $1}' |uniq ) >Network3429TDOrthologsAnnot
```
