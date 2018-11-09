# Trying to get full length transposon and retrotransposon models

Does repeat family 976 have any other hallmarks of a transposon?

### LTR finder

```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/23_LTR_finder
/work/GIF/software/programs/ltrfinder/1.0.5/ltr_finder genome738sl.polished.mitoFixed.fa -w 1 2 >ltrFinderTable &

#checked to see if repeat family 976 gene (16B09 alignment) is part of an ltr retroelement.  it is not.

#making gff coordinates for ltr retroelements
paste <(grep -v "Sequence:" ltrFinderTable |grep "scaffold") <(less ltrFinderTable |grep "Location" | awk '{print $3,$5,$8}' |sed 's/Strand://g' |less) |awk '{print  $2,"ltrfinder","gene",$4,$5,$6,".",".",$2":"$4"-"$5}'|tr " " "\t" |sort|uniq> LtrRetroelement.gff

#How many potential LTR retroelements are there?
wc -l LtrRetroelement.gff
6220 LtrRetroelement.gff

#How many elements are there if we remove extend the overlaps that intersect?
bedtools merge -i <(sort -k1,1V -k4,5n LtrRetroelement.gff)  >LtrRetroelementMerge.gff
wc -l LtrRetroelementMerge.gff
569 LtrRetroelementMerge.gff

#how many genes are found within ltr retroelements?
 bedtools intersect -wo -a NonRedundantLtrRetroelementMerge.gff -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq |wc
   1401    1401   29421
 bedtools intersect -wo -a NonRedundantLtrRetroelementMerge.gff -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq |sed 's/ID=//g' >NonRedundantLtrRetroelementMerge.list

#How many effectors are in LTR retroelements?
bedtools intersect -wo -a NonRedundantLtrRetroelementMerge.gff -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq|sed 's/ID=//g' |grep -w -f - ../18_effectorRedo/4SebastianEffvsGenes.list |wc
      6      18     429
#what are they? really 5
bedtools intersect -wo -a LtrRetroelement.gff -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq|sed 's/ID=//g' |grep -w -f - ../18_effectorRedo/4SebastianEffvsGenes.list
====================================================================================
Hetgly.G000004589 GLAND14 Pioneer,endopeptidase
Hetgly.G000004590 GLAND14 Pioneer,endopeptidase
Hetgly.G000011436 GLAND5 Pioneer,11A06/2D01/24A12/16B09/30E03/22C12/4D06/29D09family
Hetgly.G000011440 GLAND5 Pioneer,11A06/2D01/24A12/16B09/30E03/22C12/4D06/29D09family
Hetgly.G000026369 4D06 Pioneer,29D09/11A06/2D01/16B09/30E03/22C12/24A12family
Hetgly.G000026369 GLAND6 Pioneer,4D06/29D09/11A06/2D01/16B09/30E03/22C12/24A12family
=====================================================================================

#How many genes in the effectorome (292 genes) are found in LTR retrotransposons
bedtools intersect -wo -b ../17_EffectorProteinNetwork/Intermixed.network.gff -a NonRedundantLtrRetroelementMerge.gff |cut -f 18 |sort|uniq|wc
     19      19     342

#How many tandemly duplicated genes are found in ltr retroelements?
bedtools intersect -wo -a NonRedundantLtrRetroelementMerge.gff -b <(grep -f ../10_tandemDups/tandem.gene.list ../1_genomeNgff/fixed.augustus.gff3 ) |awk '$12=="gene" {print $18 }' |sort|uniq|wc
    453     453   18393

#How many repeats are found within these LTR retroelements?
bedtools intersect -wo -a NonRedundantLtrRetroelementMerge.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200' |wc
   2615   57530  458701


#How many unique repeatmodeler repeats that are longer than 200bp and have more than one copy associated associated with retroelements.
bedtools intersect -wo -a NonRedundantLtrRetroelementMerge.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sort|uniq -c |sort -k1,1nr |awk '$1 >1' |wc
    255     510    8540

#creating grep index to make things faster
bedtools intersect -wo -a NonRedundantLtrRetroelementMerge.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sort|uniq -c |sort -k1,1nr |awk '$1 >1 {print $2}'|sed 's/"//g' |sed 's/Motif://g' |grep -w -f -  <(less ../24_IRF_DNATrans/consensi.fa.classified |grep ">" |sed 's/>//g' |awk '{print $1}' |sed 's/#/\t/g' ) >ClassifiedIndex


#those repeatmodeler elements are now classified for the bedtools overlap.
bedtools intersect -wo -a NonRedundantLtrRetroelementMerge.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sed 's/"//g' |sed 's/Motif://g' |while read line; do grep -w $line  ClassifiedIndex;done |paste <(bedtools intersect -wo -a NonRedundantLtrRetroelementMerge.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200' ) - >LtrRetroelementClassified.gff


#How many of the ltr retroelements have repeatmodeler support by SINE/LINE/LTR associations.
cat <(grep "LINE" NonRedundantLtrRetroelementMerge.gff ) <(grep "SINE" NonRedundantLtrRetroelementMerge.gff) <(grep "LTR" NonRedundantLtrRetroelementMerge.gff) |cut -f 1-9 |sort|uniq|wc
    569    5121   40823


#Is the above appropriate for determining if the elements are real? How many of each classified element is found between a pair of TIRs?
==============================================================================
  167 "Motif:rnd-4_family-436"        rnd-4_family-436        LTR/Gypsy
     96 "Motif:rnd-5_family-6085"       rnd-5_family-6085       LTR/Gypsy
     78 "Motif:rnd-3_family-113"        rnd-3_family-113        Unknown
     74 "Motif:rnd-4_family-946"        rnd-4_family-946        LTR/Gypsy
     56 "Motif:rnd-5_family-15104"      rnd-5_family-15104      Unknown
     47 "Motif:rnd-3_family-42" rnd-3_family-42 Unknown
     46 "Motif:rnd-5_family-1842"       rnd-5_family-1842       LTR/Gypsy
     42 "Motif:rnd-5_family-218"        rnd-5_family-218        Unknown
     40 "Motif:rnd-4_family-352"        rnd-4_family-352        Unknown
     37 "Motif:rnd-5_family-1919"       rnd-5_family-1919       LTR/Gypsy
     34 "Motif:rnd-5_family-11759"      rnd-5_family-11759      LTR/Pao
     33 "Motif:rnd-5_family-3897"       rnd-5_family-3897       LTR/Pao
     31 "Motif:rnd-4_family-1443"       rnd-4_family-1443       Unknown
     29 "Motif:rnd-5_family-2116"       rnd-5_family-2116       LTR/Gypsy-Cigr
     29 "Motif:rnd-5_family-7694"       rnd-5_family-7694       LTR/Pao
     27 "Motif:rnd-5_family-962"        rnd-5_family-962        Unknown
     26 "Motif:rnd-5_family-498"        rnd-5_family-498        Unknown
     26 "Motif:rnd-5_family-7050"       rnd-5_family-7050       LTR/Copia
     25 "Motif:rnd-5_family-2420"       rnd-5_family-2420       Unknown
     24 "Motif:rnd-4_family-24" rnd-4_family-24 Unknown
     22 "Motif:rnd-4_family-1273"       rnd-4_family-1273       LINE/CR1
     22 "Motif:rnd-4_family-1754"       rnd-4_family-1754       LTR/Pao
     21 "Motif:rnd-5_family-2771"       rnd-5_family-2771       LTR/Pao
     20 "Motif:rnd-3_family-228"        rnd-3_family-228        Unknown
     20 "Motif:rnd-5_family-993"        rnd-5_family-993        LTR/Gypsy
     19 "Motif:rnd-4_family-111"        rnd-4_family-111        DNA/hAT-hAT19
     19 "Motif:rnd-4_family-2211"       rnd-4_family-2211       LTR/Pao

==============================================================================
yes.  All ltr retroelements have been confirmed by retroelement designations


#How many duplicated buscos are found in ltr retroelements
bedtools intersect -wo -a <(cut -f 1-9 NonRedundantLtrRetroelementMerge.gff) -b <(cut -f 1-9 ../1_genomeNgff/DuplicatedBuscos.gff3) |cut -f 18 |sed 's/;/\t/g' |awk '{print $2}' |sort|uniq|wc
     22      22     550
```

### Inverted repeat finder, an attempt to distinguish DNA transposons
```
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/24_IRF_DNATrans
irf genome738sl.polished.mitoFixed.fa  2 3 5 80 10 40 500000 10000 -t7 20000 -ngs


#what is the most transposon dense scaffold and how many possible TIRs are there?
less testIRF |awk '{print $1,$4}' |grep -v "A" |grep -v "C" |grep -v "T" |grep -v "G" |sed 's/@/>/g' |tr "\n" "\t" |sed 's/>/\n>/g' |awk 'NF>1'|awk 'NF==185{print $0}' |sort -k1,1n |less


#This one has the most at 185 fields, or 184 pairs of ltrs.
>scaffold_19    10852 19825     18066 22417     16296 24151     24151 33706     36667 40349     34651 40357     33887 41289     35359 44930     35100 44930     36668 47021     73564 74002     88215 93483     95984 100692    93293 101652
    103466 103931   100602 109186   109186 112535   112535 117198   120879 121343   126966 131371   126408 136166   142422 142867   142034 148694   148265 152505   152505 156196   156419 165432   165905 169227   182472 184341   181546 186075   223779 230204   256453 258838   272652 273176   280490 284390   277405 286854   280913 289754   290717 291662   301715 312105   325990 326435   325602 332264   345965 354303   348845 356226   420520 427791   421207 429632   441826 446399   446988 457655   457655 463602   457889 478460   473595 493388   503121 508224   515428 520940   521949 522955   520784 523529   521949 528339   520806 528944   531186 536689   529449 538380   528344 538832   521511 541926   545887 551409   544758 553687   543498 553992   554511 560882   559920 560917   553993 562121   559370 562121   559900 566563   558889 567331   579683 585790   586799 587805   585634 588379   586799 593181   585656 593782   596037 601508   594290 603192   593186 603644   586324 606778   610688 616217   609560 618495   608299 618800   619319 625706   624737 625737   624189 626935   618801 626935   640497 641425   634797 642725   642231 650974   647219 650977   645197 654914   650977 656839   659234 659759   689655 716646   771220 773513

#The number of transposons that are still partially overlapping
less testIRF |awk '{print $1,$4}' |grep -v "A" |grep -v "C" |grep -v "T" |grep -v "G" |sed 's/@/>/g' |tr "\n" "\t" |sed 's/>/\n>/g' |awk 'NF>1' |awk '{print $1,$2,$3"\n"$1,$4,$5"\n"$1,$6,$7"\n"$1,$6,$7"\n"$1,$8,$9"\n"$1,$10,$11"\n"$1,$12,$13"\n"$1,$14,$15"\n"$1,$16,$17"\n"$1,$18,$19"\n"$1,$20,$21"\n"$1,$22,$23"\n"$1,$24,$25"\n"$1,$26,$27"\n"$1,$28,$29"\n"$1,$30,$31"\n"$1,$32,$33"\n"$1,$34,$35"\n"$1,$36,$37"\n"$1,$38,$39"\n"$1,$40,$41"\n"$1,$42,$43"\n"$1,$44,$45"\n"$1,$46,$47"\n"$1,$48,$49"\n"$1,$50,$51"\n"$1,$52,$53"\n"$1,$54,$55"\n"$1,$56,$57"\n"$1,$58,$59"\n"$1,$60,$61"\n"$1,$62,$63"\n"$1,$64,$65"\n"$1,$66,$67"\n"$1,$68,$69"\n"$1,$70,$71"\n"$1,$72,$73"\n"$1,$74,$75"\n"$1,$76,$77"\n"$1,$78,$79"\n"$1,$80,$81"\n"$1,$82,$83"\n"$1,$84,$85"\n"$1,$86,$87"\n"$1,$88,$89"\n"$1,$90,$91"\n"$1,$92,$93"\n"$1,$94,$95"\n"$1,$96,$97"\n"$1,$98,$99"\n"$1,$100,$101"\n"$1,$102,$103"\n"$1,$104,$105"\n"$1,$106,$107"\n"$1,$108,$109"\n"$1,$110,$111"\n"$1,$112,$113"\n"$1,$114,$115"\n"$1,$116,$117"\n"$1,$118,$119"\n"$1,$120,$121"\n"$1,$122,$123"\n"$1,$124,$125"\n"$1,$126,$127"\n"$1,$128,$129"\n"$1,$130,$131"\n"$1,$132,$133"\n"$1,$134,$135"\n"$1,$136,$137"\n"$1,$138,$139"\n"$1,$140,$141"\n"$1,$142,$143"\n"$1,$144,$145"\n"$1,$146,$147"\n"$1,$148,$149"\n"$1,$150,$151"\n"$1,$152,$153"\n"$1,$154,$155"\n"$1,$156,$157"\n"$1,$158,$159"\n"$1,$160,$161"\n"$1,$162,$163"\n"$1,$164,$165"\n"$1,$166,$167"\n"$1,$168,$169"\n"$1,$170,$171"\n"$1,$172,$173"\n"$1,$174,$175"\n"$1,$176,$177"\n"$1,$178,$179"\n"$1,$180,$181"\n"$1,$182,$183}' |awk 'NF>1' |sort|uniq|awk '{print $1,"IRF","gene",$2,$3,"+",".",".",$1":"$2"-"$3}' |sed 's/>//g' |tr " " "\t"|awk 'NF>1' |wc
   2345    7035   62630

#merging overlapping DNA transposons
bedtools merge -i <(sort -k1,1V -k4,5n IRF.gff) >IRFMerge.bed
awk '{if($2==0){print $1,"IRFFinder","gene","1",$3,"+",".",".",$1":"$2"-"$3}else {print $1,"IRFFinder","gene",$2,$3,"+",".",".",$1":"$2"-"$3}}' IRFMerge.bed|tr " " "\t"  >IRFMerge.gff

#how many distinct DNA transposons are there determined from IRF finder output.
wc IRFMerge.bed
1075  3225 27635 IRFMerge.bed


#How many overlaps are there with a repeatmodeler element longer than 200bp?
bedtools intersect -wo -a IRFMerge.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200' |wc
   4620  101640  810991



#how many repeatmodeler elements are associated with these TIRs more than as single time?
bedtools intersect -wo -a IRFMerge.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sort|uniq -c |sort -k1,1nr |awk '$1 >1' |wc     
    356     712   11936

#creating grep index to make things faster
bedtools intersect -wo -a IRFMerge.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sort|uniq -c |sort -k1,1nr |awk '$1 >1 {print $2}'|sed 's/"//g' |sed 's/Motif://g' |grep -w -f -  <(less consensi.fa.classified |grep ">" |sed 's/>//g' |awk '{print $1}' |sed 's/#/\t/g' ) >ClassifiedIndex


#those repeatmodeler elements are now classified for the bedtools overlap.
bedtools intersect -wo -a RoundedLargeForEdits.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sed 's/"//g' |sed 's/Motif://g' |while read line; do grep -w $line  ClassifiedIndex;done |paste <(bedtools intersect -wo -a RoundedLargeForEdits.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200' ) - >RoundedLargeForEditsClassified.gff

#How many pairs of TIRs have known DNA elements?
less IRFMergeClassified.gff |grep "DNA" |cut -f 1-9 |sort|uniq |wc
    429    3861   31021

#Is the above appropriate for determining if the elements are real? How many of each classified element is found between a pair of TIRs?
bedtools intersect -wo -a IRFMergeClassified.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sort|uniq -c |sort -k1,1nr |awk '$1 >1{print $2}' |sed 's/"//g' |sed 's/Motif://g' |while read line; do grep -w $line ClassifiedIndex ;done |paste <(bedtools intersect -wo -a IRFMergeClassified.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sort|uniq -c |sort -k1,1nr |awk '$1 >1' ) -|less

==================================================================================
    4754 "Motif:rnd-4_family-1273"       rnd-4_family-1273       LINE/CR1
   2117 "Motif:rnd-3_family-42" rnd-3_family-42 Unknown
   1223 "Motif:rnd-4_family-352"        rnd-4_family-352        Unknown
   1164 "Motif:rnd-4_family-65" rnd-3_family-228        Unknown
    933 "Motif:rnd-5_family-6399"       rnd-4_family-65 Unknown
    874 "Motif:rnd-4_family-976"        rnd-4_family-111        DNA/hAT-hAT19
    722 "Motif:rnd-4_family-4752"       rnd-5_family-962        Unknown
    677 "Motif:rnd-3_family-228"        rnd-4_family-24 Unknown
    534 "Motif:rnd-4_family-111"        rnd-3_family-7  DNA/hAT-hAT19
    504 "Motif:rnd-5_family-2773"       rnd-5_family-6050       Unknown
    336 "Motif:rnd-4_family-928"        rnd-3_family-146        Unknown
    318 "Motif:rnd-3_family-690"        rnd-3_family-847        Unknown
    259 "Motif:rnd-5_family-922"        rnd-3_family-964        Unknown
    222 "Motif:rnd-4_family-256"        rnd-4_family-265        DNA/Sola-2
    220 "Motif:rnd-5_family-11767"      rnd-4_family-976        Simple_repeat
    210 "Motif:rnd-4_family-861"        rnd-4_family-299        Unknown
    191 "Motif:rnd-5_family-7921"       rnd-5_family-2862       Unknown
    190 "Motif:rnd-5_family-2980"       rnd-2_family-24 Unknown
    181 "Motif:rnd-4_family-2103"       rnd-5_family-7376       Unknown
    144 "Motif:rnd-4_family-946"        rnd-5_family-606        DNA/hAT-Tip100
    138 "Motif:rnd-5_family-218"        rnd-4_family-1  DNA/hAT-hAT19
    137 "Motif:rnd-4_family-138"        rnd-5_family-3719       DNA/CMC-Transib
    125 "Motif:rnd-4_family-1443"       rnd-4_family-115        DNA/hAT-hATx
    117 "Motif:rnd-4_family-130"        rnd-4_family-333        Unknown
    116 "Motif:rnd-4_family-1978"       rnd-4_family-2996       Unknown
    111 "Motif:rnd-4_family-149"        rnd-5_family-13 Unknown
    110 "Motif:rnd-4_family-1330"       rnd-5_family-1534       Unknown
    108 "Motif:rnd-5_family-5257"       rnd-4_family-268        Unknown
    107 "Motif:rnd-4_family-570"        rnd-5_family-3177       Unknown
    102 "Motif:rnd-4_family-1054"       rnd-4_family-138        Unknown
    100 "Motif:rnd-5_family-2129"       rnd-4_family-301        Unknown
     94 "Motif:rnd-3_family-980"        rnd-4_family-1238       Unknown
     94 "Motif:rnd-4_family-436"        rnd-4_family-1265       LINE/L1-Tx1
     83 "Motif:rnd-5_family-1286"       rnd-4_family-256        DNA/MuLE-MuDR
     82 "Motif:rnd-3_family-847"        rnd-4_family-72 Unknown
     79 "Motif:rnd-5_family-3897"       rnd-5_family-2859       Unknown
     64 "Motif:rnd-4_family-599"        rnd-4_family-354        Unknown
     62 "Motif:rnd-4_family-71" rnd-5_family-2129       Unknown
     62 "Motif:rnd-5_family-3559"       rnd-4_family-1408       Unknown

=============================================================================
#no, going by known DNA elements does not seem to be appropriate.  There are lots of unknown that are enriched between TIRs.  Interestingly, rnd-4-family976 is 11th most enriched.

#Extracting all repeatmodeler elements that have more than 10 copies in the putative DNA elements.  Used those to grep the classified GFF.  
#How many DNA transposons now?
bedtools intersect -wo -a IRFMergeClassified.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sort|uniq -c |sort -k1,1nr |awk '$1 >1{print $2}' |sed 's/"//g' |sed 's/Motif://g' |while read line; do grep -w $line ClassifiedIndex ;done |paste <(bedtools intersect -wo -a IRFMergeClassified.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sort|uniq -c |sort -k1,1nr |awk '$1 >1' ) -|awk '$1>10' |awk '{print $3}' |grep -w -f - IRFMergeClassified.gff |cut -f 1-9|sort|uniq|wc
    897    8073   64903


bedtools intersect -wo -a IRFMergeClassified.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sort|uniq -c |sort -k1,1nr |awk '$1 >1{print $2}' |sed 's/"//g' |sed 's/Motif://g' |while read line; do grep -w $line ClassifiedIndex ;done |paste <(bedtools intersect -wo -a IRFMergeClassified.gff -b ../2_repeats/genome738sl.polished.mitoFixed.fa.out.gff |grep -v "(" |grep -v "rich" |awk '$22>200 {print $19}' |sort|uniq -c |sort -k1,1nr |awk '$1 >1' ) -|awk '$1>10' |awk '{print $3}' |grep -w -f - IRFMergeClassified.gff |cut -f 1-9|sort|uniq >SupportedIRFMergeClassified.gff

#Number of genes found in the supported DNA transposons
bedtools intersect -wo -a SupportedIRFMergeClassified.gff -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq|wc
   1915    1915   40215
 bedtools intersect -wo -a SupportedIRFMergeClassified.gff -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq|sed 's/ID=//g' >SupportedIRFMergeClassified.list

#How many duplicated buscos are found in DNA transposons
bedtools intersect -wo -a <(cut -f 1-9 IRFMergeClassified.gff) -b <(cut -f 1-9 ../1_genomeNgff/DuplicatedBuscos.gff3) |cut -f 18 |sed 's/;/\t/g' |awk '{print $2}' |sort|uniq|wc
     29      29     725



#extracting fasta for analysis   
bedtools intersect -wo -a SupportedIRFMergeClassified.gff -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene" {print $18}'  |sed 's/;/\t/g' |cut -f 1 |sed 's/ID=//g' |sort|uniq|cdbyank ../1_genomeNgff/gff2fasta.gene.fasta.cidx >GenesInDNATransposons.fasta
 cd-hit-est -i GenesInDNATransposons.fasta -d 30 -o GenesInDNATransposons.cd-hit -c .8
1915 finished   1364 clusters

bedtools intersect -F .8 -wo -a SupportedIRFMergeClassified.gff  -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene" {print $18}'  |sed 's/;/\t/g' |cut -f 1 |sed 's/ID=//g' |sort|uniq|awk '{print $1".t1"}' | cdbyank ../1_genomeNgff/augustus.aa.cidx > ProteinsInDNATransposons.fasta

cd-hit -i ProteinsInDNATransposons.fasta -d 30 -o ProteinsInDNATransposons.cd-hit -c .8 &
1409 finished   1055 clusters


#How many of these genes have also been called as tandemly duplicated.
bedtools intersect -wo -a SupportedIRFMergeClassified.gff -b <(grep -f ../10_tandemDups/tandem.gene.list ../1_genomeNgff/fixed.augustus.gff3 ) |awk '$12=="gene" {print $18 }' |sort|uniq|wc
    710     710   28882

#so about a seventh of the genes in tandem duplicates are
wc ../10_tandemDups/tandem.gene.list
 4774  4774 85932 ../10_tandemDups/tandem.gene.list

#How many genes from the effectorome are found in DNA transposons?  Lets make the gff first
#/work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/17_EffectorProteinNetwork
 sed 's/,/\t/g' ../17_EffectorProteinNetwork/Intermixed.network|cut -f 2 |sort|uniq |grep "Hetgly" |wc
    292     292    5256
sed 's/,/\t/g' Intermixed.network |cut -f 2 |sort|uniq|sed '/^$/d'|grep -w -f - <(sed 's/ID=/\t/g' ../1_genomeNgff/fixed.augustus.gff3 |sed 's/;/\t/g' ) |awk '$3=="gene"' |cut -f 1-10|sed 's/\t\t/\t/g' >Intermixed.network.gff
less Intermixed.network.gff |awk '$3=="gene"' |cut -f 9 |sort|uniq|wc
    292     292    5256
#How many gene from effectorome are found in DNA transposons?
bedtools intersect -wo -b ../17_EffectorProteinNetwork/Intermixed.network.gff -a SupportedIRFMergeClassified.gff |cut -f 18 |sort|uniq|wc
     52      52     936


#which effectors are found in DNA transposons?
 bedtools intersect -wo -a SupportedIRFMergeClassified.gff -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq|sed 's/ID=//g' |grep -w -f - ../18_effectorRedo/4SebastianEffvsGenes.list |less

=====================================================================================
Hetgly.G000003427 4D06 Pioneer,29D09/11A06/2D01/16B09/30E03/22C12/24A12family
Hetgly.G000003427 GLAND6 Pioneer,4D06/29D09/11A06/2D01/16B09/30E03/22C12/24A12family
Hetgly.G000003429 4D06 Pioneer,29D09/11A06/2D01/16B09/30E03/22C12/24A12family
Hetgly.G000003430 GLAND6 Pioneer,4D06/29D09/11A06/2D01/16B09/30E03/22C12/24A12family
Hetgly.G000006478 11A06 Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family
Hetgly.G000008340 22C12 Pioneer,16B09/30E03/11A06/2D01/24A12/4D06/29D09family
Hetgly.G000009412 11A06 Pioneer,2D01/24A12/16B09/30E03/22C12/4D06/29D09family
Hetgly.G000009412 24A12 Pioneer,2D01/11A06/16B09/22C12/30E03/4D06/29D09family
Hetgly.G000009412 2D01 Pioneer,11A06/24A12/16B09/30E03/22C12/4D06/29D09family
Hetgly.G000011436 GLAND5 Pioneer,11A06/2D01/24A12/16B09/30E03/22C12/4D06/29D09family
Hetgly.G000011440 GLAND5 Pioneer,11A06/2D01/24A12/16B09/30E03/22C12/4D06/29D09family
Hetgly.G000012497 16A01 Pioneer,30D08/21E12family
Hetgly.G000012497 21E12 Pioneer,30D08/16A01family
Hetgly.G000016164 10A07 Pioneer,20G04/13A06/6E07/32E03family
Hetgly.G000016164 20G04 Pioneer,10A07/27D09/13A06/6E07/32E03family
Hetgly.G000016164 27D09 Pioneer,10A07/20G04/10A07/6E07/13A06/32E03family
Hetgly.G000017662 17G01 Pioneer
Hetgly.G000018262 32E03 Pioneer,10A07/6E07/13A06/27D09/20G04family
Hetgly.G000020629 GLAND5 Pioneer,11A06/2D01/24A12/16B09/30E03/22C12/4D06/29D09family
Hetgly.G000022993 25A01 Pioneer,30G12/4G05family
Hetgly.G000022993 30G12 Pioneer,4G06/25A01family
Hetgly.G000022993 4G05 Pioneer,30G12/25A01family
Hetgly.G000024438 16A01 Pioneer,30D08/21E12family
Hetgly.G000024438 21E12 Pioneer,30D08/16A01family
Hetgly.G000024438 30D08 Pioneer,16A01/21E12family
Hetgly.G000024443 30D08 Pioneer,16A01/21E12family
Hetgly.G000024447 30D08 Pioneer,16A01/21E12family
===========================================================================================

#how many effector gene loci?
bedtools intersect -wo -a SupportedIRFMergeClassified.gff -b ../1_genomeNgff/fixed.augustus.gff3 |awk '$12=="gene"' |cut -f 18 |sed 's/;/\t/g' |awk '{print $1}' |sort|uniq|sed 's/ID=//g' |grep -w -f - ../18_effectorRedo/4SebastianEffvsGenes.list |awk '{print $1}' |sort|uniq|wc
     17      17     306
```
