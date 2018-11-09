### The introns are missing and some of the annotations are also.  Need to fix the 738 genome gff.

```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/44_FixGFFforJbrowse

/work/GIF/software/programs/maker/2.31.9/bin/maker_map_ids --prefix HetGly. --justify 5 augustus.gff3  >augustus.map

# had to create a new map, because I wanted the gene names to remain the same
sed 's|/|\t|g' ../1_genomeNgff/geneRenamer.sh |cut -f 2,3 |cat <(grep "\.t" augustus.map|cut -f 1) - |sort -k1,1V |awk -v col2=0 '{if(NF==2) {col2=$2;print $0} else {print $0,col2}}'  |tr " " "\t" |sed 's/\./\t/1' |awk '{if($2!="Hetgly") {print $1"."$2,$3"."$2} else {print $1,$2,$3 }}'  |awk '{if(NF>2) {print $1,$2"."$3} else {print $0}}' |tr " " "\t" >Augustus.Improved.map


```
