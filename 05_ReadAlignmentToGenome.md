# Converting blasr ccs read and subread mapping to gff3
```

```
##  ccsreads
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/36_deduplicate122016/ccsMapping/blasrAlignment/ccs2genome.out

module load blasr
blasr reads_of_insert.fasta genome738sl.polished.mitoFixed.fa --nproc 16 --nCandidates 1 -m 4  --out ccs2genome.out --bestn 1 --unaligned  ccs2genomeUnaligned.fasta

less /work/GIF/remkv6/Baum/CamTechGenomeComparison/36_deduplicate122016/ccsMapping/blasrAlignment/ccs2genome.out |awk '{if($9==0){print $2,"Blasr","CCSRead",$10,$11,".","+",".","ID="$1";"} else {print $2,"Blasr","CCSRead",$10,$11,".","-",".","ID="$1";"}}' |tr " " "\t">ccs2genome.gff

/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/35_ReadAlignments/ccs
awk '{print $1}' ccs2genome.gff|sort -k1,1V > NamesSorted11V
#ran scaffoldRenamer.sh
paste NamesSorted11V <(sort -k1,1V ccs2genome.gff|awk '{print $2,$3,$4,$5,$6,$7,$8,$9}' ) |tr " " "\t" >ccs2genomeRenamed.gff
 sh ~/common_scripts/runTabix.sh ccs2genomeRenamed.gff

```
## preads
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/36_deduplicate122016/preadMapping/blasr
blasr preads4falcon.fasta genome738sl.polished.mitoFixed.fa --nproc 16 --nCandidates 1 -m 4  --out preads2genome.out --bestn 1 --unaligned  preads2genomeUnaligned.fasta

less preads2genome.out |awk '{if($9==0){print $2,"Blasr","PRead",$10,$11,".","+",".","ID="$1";"} else {print $2,"Blasr","[PRead",$10,$11,".","-",".","ID="$1";"}}' |tr " " "\t">preads2genome.gff

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/35_ReadAlignments/pread
awk '{print $1}' preads2genome.gff|sort -k1,1V > NamesSorted11V
paste NamesSorted11V <(sort -k1,1V preads2genome.gff|awk '{print $2,$3,$4,$5,$6,$7,$8,$9}' ) |tr " " "\t" >preads2genomeRenamed.gff
 sh ~/common_scripts/runTabix.sh preads2genomeRenamed.gff
```
## subreads
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/36_deduplicate122016/subreadMapping/blasrAlignment
less subreads2genome.out |awk '{if($9==0){print $2,"Blasr","SubRead",$10,$11,".","+",".","ID="$1";"} else {print $2,"Blasr","SubRead",$10,$11,".","-",".","ID="$1";"}}' |tr " " "\t">subreads2genome.gff

#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/35_ReadAlignments/subread
awk '{print $1}' subreads2genome.gff|sort -k1,1V > NamesSorted11V
#ran scaffoldRenamer.sh
paste NamesSorted11V <(sort -k1,1V subreads2genome.gff|awk '{print $2,$3,$4,$5,$6,$7,$8,$9}' ) |tr " " "\t" >subreads2genomeRenamed.gff
 sh ~/common_scripts/runTabix.sh subreads2genomeRenamed.gff
```
