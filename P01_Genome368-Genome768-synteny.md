# Repeatmasker, opscan, and iadhore needed redone to include all scaffolds

### Input files for opscan
```
#/data021/GIF/remkv6/CamTechGenomeComparison/06_iadhore_allscaf/opscan\\
#generating for 368 genome, operon.db
ln -s ../../augustus.368.gff3
ln -s ../../augustus368.aa
sed -e 's/ID=//g' augustus.368.gff3 -e 's/;/\t/g' -e 's/NODE_//g'|grep gene -|awk '{print $3 "\t" $9 "\t" substr($1, 0, 8) "\t" $4 "\t" $5 "\t" $7 "\t" "t"}' |sed 's/scaffold_//g' |sort -k 5 -V|awk '{print ">"$2 "\t" "1" "\t" $4 "\t" $5 "\t" $6 "\t" $3}'|sort -k 1 -V |sed -e 's/F+//g' -e 's/F-//g' -e 's/F//g' >augustus.368.gff3.reformat1
tr "\n" "\t" < augustus368.aa |sed 's/>/\n>/g'|sed 's/\t//g'|grep ".t1"|sed  's/.t1/\t/g'>3684join2
join -1 1 -2 1 augustus.368.gff3.reformat1 3684join2 >joined.386.3
sed 's/ /\n/6' joined.386.3 >operon.db

#generating for 2692 genome, genome.db
#used elegans masked 2692 genome that was used to generate the .aa and .gff3 files.
#this was previously generated and not recorded but the same method was used below.
#sed -e 's/ID=//g' G2692_Hg_augustus.gff3 -e 's/;/\t/g' -e 's/NODE_//g'|grep gene -|awk '{print $3 "\t" $9 "\t" substr($1, 0, 8) "\t" $4 "\t" #$5 "\t" $7 "\t" "t"}' |sed 's/scaffold_//g' |sort -k 5 -V|awk '{print ">" $2 "\t" "1" "\t" $4 "\t" $5 "\t" $6 "\t" $3}'|sort -k 1 -V |sed -e 's/F+//g' -e 's/F-//g' -e 's/F//g' >augustus_G2692.gff3.reformat1
#tr "\n" "\t" <G2692_Hg_augustus.aa |sed 's/>/\n>/g'|sed 's/\t//g'|grep ".t1"|sed  's/.t1/\t/g'>26924join2
#join -1 1 -2 1 augustus_G2692.gff3.reformat1 26924join2 >joined.2692.3
#sed 's/ /\n/6' joined.2692.3 >genome.db

#Same opscan pbs script as previous, default parameters
```

### File setup for iadhore

```
#/data021/GIF/remkv6/CamTechGenomeComparison/06_iadhore_allscaf/Nematode
cp ../opscan/genome.db .
cp ../opscan/operon.db .

#368 first, generating lst files and 368.ini
/data021/GIF/remkv6/CamTechGenomeComparison/06_iadhore_allscaf/Nematode/368
tr "\n" "\t" <../operon.db |sed 's/>/\n>/g'|awk '{print $1$5 " " $6}'|awk '{print >> $2 ".lst"; close($2)}'
sed -i 's/\t.*//g' *.lst
sed -i 's/>//g' *.lst
ls -1 *.lst|awk '$2=$1{print $1=substr($1,0,6), " 368/"$2}'>368.ini


###This was already generated
#moving to 2692/genome.db file reformat to .lst and .ini files
#/data021/GIF/remkv6/CamTechGenomeComparison/06_iadhore_allscaf/Nematode/2692
##tr "\n" "\t" < ../genome.db |sed 's/>/\n>/g'|awk '{print $1$5 " " $6}'|awk '{print >> $2 ".lst"; close($2)}'
##sed -i 's/\t.*//g' *.lst
##sed -i 's/>//g' *.lst
##ls -1 *.lst|awk '$2=$1{print $1=substr($1,0,6), " 2692/"$2}'>2692.ini

#/data021/GIF/remkv6/CamTechGenomeComparison/06_iadhore_allscaf/Nematode
cat 2692/2692.ini 368/368.ini >Nematode.ini
#added to the Nematode.ini file
this goes at the top of the gene list for 2692 and a newline is added at the end of 2692's list before the genome= 368.
#no spaces allowed

genome=269

genome=368

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


#rearranging opscan output to form iadhore plaza table
#/data021/GIF/remkv6/CamTechGenomeComparison/06_iadhore_allscaf/opscan
##################################################################################
make sure to change the input file here next time you run

grep "CL"  OPSCAN.o104582.condo |awk '{print $3 "\t" $6}'  >iADHoRe_PLAZA.table
cp iADHoRe_PLAZA.table ../Nematode/

#the iadhore pbs script
?
1
2
3
4
5
6
7
8
9
10
11
12
13
14

#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=1:00:00
#PBS -N iadhore
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -q debug
cd $PBS_O_WORKDIR
ulimit -s unlimited
module use /shared/software/GIF/modules
module load iAdHoRe/3.0.01
i-adhore Nematode.ini

##in case you need stats after job completion retain this as last line
ssh condo "qstat -f ${PBS_JOBID} |head"
```

### Post processing

```
#/data021/GIF/remkv6/CamTechGenomeComparison/06_iadhore_allscaf/Nematode/output/forgffsubtraction
?
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51

ln -s ../../../../augustus.368.gff3 #braker run prior to masking
ln -s ../../../../05_iadhore/Nematode/output/forgffsubtraction/G2692_Hg_augustus.gff3 #braker run after masking with elegans

#attaching scaffold lengths to muliplicons.txt
bioawk -c fastx '{ print $name, length($seq) }' <genome.2692.masked.fa|sed -e 's/F+//g' -e 's/F-//g' -e 's/F//g' -e 's/NODE_//g' -e 's/R//g'|awk '{print substr ($1,0,6) "K\t" $2}' >genome.2692.masked.chr.len
bioawk -c fastx '{ print $name, length($seq) }' <genome.368.fa|sed -e 's/F+//g' -e 's/F-//g' -e 's/F//g' -e 's/NODE_//g' |awk '{print substr ($1,0,6) "\t" $2}'>genome.368.chr.len

##they have different names now, so splitting by 269 and 368 is not necessary
cat genome.368.masked.chr.len genome.2692.masked.chr.len |sort -k 1,1 >cat.genome.masked.chr.len
awk '{print $3,$0}' ../multiplicons.txt|sort -k 1,1 -V >multiplicons.forchrlen.attach
join -a 1 -1 1 -2 1 multiplicons.forchrlen.attach cat.genome.masked.chr.len >len.1.done

awk '{print $6,$0}' len.1.done |sort -k 1,1 -V >multiplicons.forchrlen.attach.2
join -a 1 -1 1 -2 1 multiplicons.forchrlen.attach.2 cat.genome.masked.chr.len |cut -d" " -f 3->len.2.done
done!
#getting info from the gff, to obtain the scaffold positions of genes
grep gene augustus_G368.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print $1 "\t" $5 "\t" $6}'>368.gene.coord
#be sure to rename.
grep gene G2692_Hg_augustus.gff3 |sed -e 's/ID=//g' -e 's/;//g'|awk '{print $9 "\t" $0}'|awk '{print "K"substr($1,2,6) "\t" $5 "\t" $6}'>2692.gene.coord
cat 2692.gene.coord 368.gene.coord >both.gene.coord

#reformatting the genes.txt files to obtain the gene number (assigned by iadhore) and attaching them to the above
awk '$2==368 {print $1 "\t" $2 "\t" $3 "\t" $4}' ../genes.txt|sort -k 1,1 -V >genes.txt.368
awk '$2==269 {print $1 "\t" $2 "\t" $3 "\t" $4}' ../genes.txt|sort -k 1,1 -V  >genes.txt.269
cat genes.txt.269 genes.txt.368 >reformat.genes.txt

#no longer necessary to split
####join -a 2  -1 1 -2 1 genes.txt.368 368.gene.coord >position.368
####join -a 2  -1 1 -2 1 genes.txt.269 2692.gene.coord >position.269
join -a 1  -1 1 -2 1 reformat.genes.txt both.gene.coord >gene.positions

#no longer necessary to separate
####awk '{print $3$4,$0}' position.269|sort -k 1,1>position.269.fuse
####awk '{print $3$4,$0}' position.368|sort -k 1,1>position.368.fuse
awk '{print $9_$3, $0}' len.2.done |sort -k 1,1>4gene1.attach
awk '{print $4$3,$0}' gene.positions|sort -k 1,1>position.forjoin
join -a 1 -1 1 -2 1  4gene1.attach position.forjoin>gene.int1.attach
cut -d " " -f 2- gene.int1.attach| awk '{print $10$3,$0}'|tr " " "\t"|sort -k 1,1>4gene2.attach
join -a 1 -1 1 -2 1  4gene2.attach position.forjoin   >gene.int2.attach
cut -d " " -f 2- gene.int2.attach| awk '{print $11$5,$0}'|tr " " "\t"|sort -k 1,1>4gene3.attach
join -a 1 -1 1 -2 1  4gene3.attach position.forjoin   >gene.int3.attach
cut -d " " -f 2- gene.int3.attach| awk '{print $12$5,$0}'|tr " " "\t"|sort -k 1,1>4gene4.attach
join -a 1 -1 1 -2 1  4gene4.attach position.forjoin   >gene.int4.attach
tr " " "\t" <gene.int4.attach |cut -f 2-|head -n -1>master.table
done!
#Headerinformation (tab separated file and header)
id      genome_x        scaffold_x  parent  genome_y        scaffold_y  level   number_of_anchorpoints  profile_length  begin_x end_x   begin_y end_y   is_redundant     genome_x_scaffold_len     genome_y_scaffold_len     genome_x     scaffold_x     relativeGeneNumber     geneStart     geneEnd     genome_x     scaffold_x     relativeGeneNumber    geneStart     geneEnd     genome_y     scaffold_y     relativeGeneNumber    geneStart     geneEnd     genome_y     scaffold_y     relativeGeneNumber    geneStart     geneEnd

#location
/data021/GIF/remkv6/CamTechGenomeComparison/05_iadhore/Nematode/output/forgffsubtraction/master.table
###WARNING -- I had to manually add the scaffold length for multiplicon 441 for scaffold 2691K
```
### attaching gene and gene positions
```


##now attaching gene and gene positions #####269 gene 1 ####awk '$2==269' scaffoldlen.attached.mult|awk '{print $9“\t”$0}'>269.scaf.len.1 ####cut -f 1- 269.scaf.len.1 |awk '{print $4$1,$0}'|sort -k 1,1 >269.scaf.len.1.fusename ####join -a 1 -1 1 -2 1 269.scaf.len.1.fusename position.269.fuse |tr “ ” “\t”|cut -f 3→269.posattach

#####368 gene 1 ####awk '$2==368' scaffoldlen.attached.mult|awk '{print $9“\t”$0}'>368.scaf.len.1 ####cut -f 1- 368.scaf.len.1 |awk '{print $4$1,$0}'|sort -k 1,1 >368.scaf.len.1.fusename ####join -a 1 -1 1 -2 1 368.scaf.len.1.fusename position.368.fuse| tr “ ” “\t”|cut -f 3→368.posattach

#####join the two files ####cat 269.posattach 368.posattach >gene1added

#####269 gene 2 ####awk '$2==269' gene1added|awk '{print $10“\t”$0}' >269.scaf.len.2 ####cut -f 1- 269.scaf.len.2 |awk '{print $4$1,$0}'|sort -k 1,1 >269.scaf.len.2.fusename ####join -a 1 -1 1 -2 1 269.scaf.len.2.fusename position.269.fuse |tr “ ” “\t”|cut -f 3→269.posattach.2

#####368 gene 2 ####awk '$2==368' gene1added|awk '{print $10“\t”$0}' >368.scaf.len.2 ####cut -f 1- 368.scaf.len.2 |awk '{print $4$1,$0}'|sort -k 1,1 >368.scaf.len.2.fusename ####join -a 1 -1 1 -2 1 368.scaf.len.2.fusename position.368.fuse| tr “ ” “\t”|cut -f 3→368.posattach.2

####join the two files ####cat 269.posattach.2 368.posattach.2 >gene2added

#####269 gene 3 ####awk '$2==269' gene2added|awk '{print $11“\t”$0}' >269.scaf.len.3 ####cut -f 1- 269.scaf.len.3 |awk '{print $6$1,$0}'|sort -k 1,1 >269.scaf.len.3.fusename ####join -a 1 -1 1 -2 1 269.scaf.len.3.fusename position.269.fuse |tr “ ” “\t”|cut -f 3→269.posattach.3

#####368 gene 3 ####awk '$2==368' gene2added|awk '{print $11“\t”$0}' >368.scaf.len.3 ####cut -f 1- 368.scaf.len.3 |awk '{print $6$1,$0}'|sort -k 1,1 >368.scaf.len.3.fusename ####join -a 1 -1 1 -2 1 368.scaf.len.3.fusename position.368.fuse| tr “ ” “\t”|cut -f 3→368.posattach.3

#####join the two files ####cat 269.posattach.3 368.posattach.3 >gene3added

#####269 gene 4 ####awk '$2==269' gene3added|awk '{print $12“\t”$0}' >269.scaf.len.4 ####cut -f 1- 269.scaf.len.4 |awk '{print $6$1,$0}'|sort -k 1,1 >269.scaf.len.4.fusename ####join -a 1 -1 1 -2 1 269.scaf.len.4.fusename position.269.fuse |tr “ ” “\t”|cut -f 3→269.posattach.4

#####368 gene 4 ####awk '$2==368' gene3added|awk '{print $12“\t”$0}' >368.scaf.len.4 ####cut -f 1- 368.scaf.len.4 |awk '{print $6$1,$0}'|sort -k 1,1 >368.scaf.len.4.fusename ####join -a 1 -1 1 -2 1 368.scaf.len.4.fusename position.368.fuse| tr “ ” “\t”|cut -f 3→368.posattach.4

#####join the two files ####cat 269.posattach.4 368.posattach.4 >gene4added

```
### Postprocessing for scaffold extension from 5' and 3' ends
```
Trying to make larger scaffolds by extending 368 scaffolds where 2692 scaffolds will create a longer overall scaffold

?
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
173
174
175
176
177
178
179
180
181
182
183
184
185
186
187
188
189
190
191
192
193
194
195
196
197
198
199
200
201
202
203
204
205
206
207
208
209
210
211
212
213
214
215
216
217
218
219
220
221
222
223
224
225
226
227
228
229
230
231
232
233
234
235
236
237
238
239
240
241
242
243
244
245
246
247
248
249
250
251
252
253
254
255
256
257
258
259
260
261
262
263
264
265
266
267
268
269
270
271
272
273
274
275
276
277
278
279
280
281
282
283
284
285
286
287
288
289
290
291
292
293
294
295
296
297
298
299
300
301
302
303
304
305
306
307
308
309
310
311
312
313
314
315
316
317
318
319
320
321
322
323
324
325
326
327
328
329
330
331
332
333
334
335
336
337
338
339
340
341
342
343
344
345
346

#/data021/GIF/remkv6/CamTechGenomeComparison/10_scaffold_extend
less  ../../05_iadhore/Nematode/output/forgffsubtraction/scaffoldgenepositioncomp |awk '{print $1,$2,($4-$3),$4,$5,$6,$7,$8,$9}'>scaffold.directional.expansion
#header info for scaffold.directional.expansion space delimited
#2692_scaffold gene_start(exp_5') gene_expansion_3' 2692_scaf_len . 368_scaffold gene_start gene_end 368_scaf_len

#get those genomes
ln -s ../genome.2692.fa
ln -s ../genome.368.fa
ln -s ../06_iadhore_allscaf/Nematode/output/forgffsubtraction/master.table



#rename fasta seq to match my processed files
sed -e '/>/s/NODE_//g' -e '/>/s/F+//' -e '/>/s/F-//g' -e '/>/s/F//g' -e '/>/s/R//g' -e  '/>/s/^\(.\{7\}\).*$/\1/' genome.368.fa >genome.368.namemod.fa
sed -e '/>/s/NODE_//g' -e '/>/s/F+//' -e '/>/s/F-//g' -e '/>/s/F//g' -e '/>/s/R//g' -e '/>/s/$/K/g' genome.2692.fa >genome.2692.namemod.fa

awk '{print $3,$20,$27,$14, ".",$5,$32,$39,$15}' master.table|tr " " "\t"|head -n -1 >basics.table

#which ones add more that 5kb to the 5' end?# this may require modification later because the precursor file was changed.
awk '$2==269{print $0, ($7-$2)}' basics.table |awk '$10<-5000'|head -n -1>5pr_exp_table
#manually removed two from this list, and reformatted as below.  the bottom one is useful, the top is not
#000689K    32693   67047   67747   .   000457  80887   114858  116787
#000041K    91496   116278  431849  .   000113  227077  252205  277827 -- this one was rearranged like here and then put in basics.table1


#which ones expand more than 5kb on the 3' end of the 368 scaffolds
tr " " "\t"< basics.table|awk '$2==269{print $0,($9-$8),($4-$3)}'  |awk '$11>$10 && $10 >5000'>3pr_exp_table
tr " " "\t"< basics.table1|awk '$2==269{print $0,($9-$8),($4-$3)}'  |awk '$11>$10 && $10 >5000'>>3pr_exp_table

#Retain all clipped scaffold portions and rename.
awk '{print "samtools faidx genome.2692.namemod.fa " $1":0-"$2">>2692.5pr.merger.fa"}' 5pr_exp_table
#making merged scaffolds 5'
samtools faidx genome.2692.namemod.fa 000252K:0-157793>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000449K:0-75702>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000388K:0-99639>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000265K:0-154524>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000928K:0-35686>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000401K:0-67734>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000522K:0-66189>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000403K:0-38208>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000761K:0-34038>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000366K:0-90353>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000213K:0-159637>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000091K:0-154256>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000871K:0-32512>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000658K:0-48409>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 002077K:0-7754>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000968K:0-37372>>2692.5pr.merger.fa
samtools faidx genome.2692.namemod.fa 000710K:0-39278>>2692.5pr.merger.fa


awk '{print "samtools faidx genome.368.namemod.fa " $6":"$7"-"$9" >>368.5pr.merger2.fa"}' 5pr_exp_table

samtools faidx genome.368.namemod.fa 000050:20390-367144 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000069:10582-333063 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000340:20209-332218 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000026:22468-526853 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000647:20929-1457963 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000249:36860-170264 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000385:20796-125902 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000276:12196-178876 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000053:12734-353492 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000025:74551-530477 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000121:94055-262346 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000145:118450-239373 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000097:4313-288060 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000278:3709-135006 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 002077:1-376172 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000280:10843-539497 >>368.5pr.merger2.fa
samtools faidx genome.368.namemod.fa 000394:4441-119412 >>368.5pr.merger2.fa




#now with 3' ends
awk '{print "samtools faidx genome.368.namemod.fa " $6":0-"$8" >>368.3pr.merger.fa"}' 3pr_exp_table
samtools faidx genome.368.namemod.fa 000184:0-628671 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000211:0-844812 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000236:0-80929 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000943:0-48886 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000382:0-106834 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000406:0-106183 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000453:0-81547 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000639:0-64507 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000567:0-775702 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000520:0-161325 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000281:0-140370 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 001130:0-28589 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000402:0-114942 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000459:0-20538 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000313:0-148058 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000364:0-138376 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000814:0-194486 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 001878:0-229865 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000239:0-169135 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000338:0-107225 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000148:0-213362 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000054:0-1981830 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000147:0-175179 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000175:0-204405 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000699:0-26841 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000105:0-274172 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 002284:0-284112 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000052:0-331657 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000216:0-44477 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000652:0-249425 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000078:0-301454 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000071:0-448107 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000113:0-252205 >>368.3pr.merger.fa
samtools faidx genome.368.namemod.fa 000432:0-504107 >>368.3pr.merger.fa




#for 2692
awk '{print "samtools faidx genome.2692.namemod.fa " $1":"$3"-"$4" >>2692.3pr.merger2.fa"}' 3pr_exp_table
samtools faidx genome.2692.namemod.fa 000255K:51687-167724 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000009K:654992-668936 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000417K:13344-111288 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 001277K:15492-39921 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000382K:105087-122296 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 001351K:31739-37465 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000453K:74559-103214 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000639K:64493-77518 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000566K:12469-68095 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000946K:38523-53783 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000281K:140331-151964 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 001130K:26989-44762 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 001004K:39926-51259 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000459K:19612-122453 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000313K:147608-156712 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000364K:138367-151988 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000339K:34374-134369 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000093K:199820-279281 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000390K:48935-113166 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000338K:107206-134426 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000148K:203361-247129 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000023K:520796-546288 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000381K:60423-111996 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000219K:37824-125270 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000699K:21687-72092 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000645K:51095-72159 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000198K:98726-191576 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000296K:112565-155710 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000216K:42002-199035 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000757K:36802-66978 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000078K:298050-310391 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000377K:61221-123743 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000041K:116278-431849 >>2692.3pr.merger2.fa
samtools faidx genome.2692.namemod.fa 000935K:46657-54357 >>2692.3pr.merger2.fa



#collecting the orphaned fragments
awk '{print "samtools faidx genome.2692.namemod.fa " $1":0-"$3" >>2692.3pr.orphan.fa"}' 3pr_exp_table
samtools faidx genome.2692.namemod.fa 000255K:0-51687 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000009K:0-654992 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000417K:0-13344 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 001277K:0-15492 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000382K:0-105087 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 001351K:0-31739 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000453K:0-74559 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000639K:0-64493 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000566K:0-12469 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000946K:0-38523 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000281K:0-140331 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 001130K:0-26989 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 001004K:0-39926 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000459K:0-19612 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000313K:0-147608 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000364K:0-138367 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000339K:0-34374 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000093K:0-199820 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000390K:0-48935 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000338K:0-107206 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000148K:0-203361 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000023K:0-520796 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000381K:0-60423 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000219K:0-37824 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000699K:0-21687 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000645K:0-51095 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000198K:0-98726 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000296K:0-112565 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000216K:0-42002 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000757K:0-36802 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000078K:0-298050 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000377K:0-61221 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000041K:0-116278 >>2692.3pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000935K:0-46657 >>2692.3pr.orphan.fa


awk '{print "samtools faidx genome.368.namemod.fa " $6":"$7"-"$9" >>368.3pr.orphan.fa"}' 3pr_exp_table
samtools faidx genome.368.namemod.fa 000184:605177-656863 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000211:224850-858514 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000236:69156-177514 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000943:44871-53915 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000382:16341-122296 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000406:79122-111908 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000453:4880-103214 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000639:21194-77518 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000567:765226-830209 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000520:137312-176564 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000281:4176-151964 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 001130:104-44762 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000402:108563-120764 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000459:10051-122453 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000313:10042-156712 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000364:31312-151988 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000814:173315-291796 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 001878:48491-306196 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000239:142685-175931 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000338:17094-134426 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000148:193561-247129 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000054:1470342-1991137 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000147:120742-207122 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000175:170963-212039 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000699:2007-72092 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000105:245356-281157 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 002284:268327-334191 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000052:320850-361499 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000216:16990-199035 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000652:215219-262145 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000078:21774-310391 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000071:403873-465093 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000113:227077-277827 >>368.3pr.orphan.fa
samtools faidx genome.368.namemod.fa 000432:464163-511359 >>368.3pr.orphan.fa



awk '{print "samtools faidx genome.2692.namemod.fa " $1":"$3"-"$4" >>2692.5pr.orphan.fa"}' 5pr_exp_table
samtools faidx genome.2692.namemod.fa 000252K:169286-169907 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000449K:101083-104727 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000388K:118564-120484 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000265K:166988-178366 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000928K:50486-54840 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000401K:74370-115377 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000522K:91438-92777 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000403K:91155-115715 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000761K:75595-76737 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000366K:115396-126472 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000213K:174146-197045 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000091K:200050-244787 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000871K:44771-58808 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000658K:69972-75641 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 002077K:20637-21235 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000968K:61420-61422 >>2692.5pr.orphan.fa
samtools faidx genome.2692.namemod.fa 000710K:59992-71208 >>2692.5pr.orphan.fa



awk '{print "samtools faidx genome.368.namemod.fa " $6":0-"$8" >>368.5pr.orphan.fa"}' 5pr_exp_table

samtools faidx genome.368.namemod.fa 000050:0-36422 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000069:0-36238 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000340:0-39132 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000026:0-34928 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000647:0-35743 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000249:0-43497 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000385:0-50628 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000276:0-60964 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000053:0-55174 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000025:0-99595 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000121:0-108564 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000145:0-163066 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000097:0-20870 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000278:0-25270 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 002077:0-13832 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000280:0-33712 >>368.5pr.orphan.fa
samtools faidx genome.368.namemod.fa 000394:0-29719 >>368.5pr.orphan.fa



--------------------------------------------------------------------------------------double overlaps done manually
001310K 14904   21001   38806   .       000508  73929   80041   89362 9321 17805
001580K 13404   32383   32695   .       000508  5885    24863   89362 -7519
000233K 118569  134007  155419  .       000169  1       15970   518262 -118568
000590K 2203    40580   83237   .       000169  446633  484239  518262 34023 42657
000378K 47801   106824  123454  .       000140  990     58626   222655 -46811
000378K 5654    26470   123454  .       000140  125373  146162  222655 76493 96984
000069K 72028   243782  333063  .       000036  274073  426553  444794 18241 89281
000234K 115105  162759  163747  .       000036  91973   140733  444794 -23132
000525K 12050   28762   92097   .       000393  101486  118198  125516 7318 63335
001006K 31531   43663   51209   .       000393  5510    17640   125516 -26021

samtools faidx genome.2692.namemod.fa 001580K:0-13404 >>complex1
samtools faidx genome.368.namemod.fa 000508:5885-80041 >>complex1
samtools faidx genome.2692.namemod.fa 001310K:21001-38806 >>complex1

samtools faidx genome.2692.namemod.fa 000233K:0-118569 >>complex2
samtools faidx genome.368.namemod.fa 000169:0-484239 >>complex2
samtools faidx genome.2692.namemod.fa 000590K:40580-83237 >>complex2

samtools faidx genome.2692.namemod.fa 000378K:0-47801 >>complex3
samtools faidx genome.368.namemod.fa 000140:990-146162 >>complex3
samtools faidx genome.2692.namemod.fa 000378K:26470-123454 >>complex3

samtools faidx genome.2692.namemod.fa 000234K >complex4
samtools faidx genome.368.namemod.fa 000036 >>complex4
samtools faidx genome.2692.namemod.fa 000069K >>complex4

samtools faidx genome.2692.namemod.fa 001006K:0-31531 >complex5
samtools faidx genome.368.namemod.fa 000393:5510-118198 >>complex5
samtools faidx genome.2692.namemod.fa 000525K:28762-92097 >>complex5

samtools faidx genome.2692.namemod.fa 001310K:0-14904 >complex.orphan
samtools faidx genome.2692.namemod.fa 000233K:134007-155419 >>complex.orphan
samtools faidx genome.2692.namemod.fa 000378K:26470-47801 >>complex.orphan
samtools faidx genome.368.namemod.fa 000169:484239-518262  >>complex.orphan
samtools faidx genome.368.namemod.fa 000140:146162-222655 >>complex.orphan
samtools faidx genome.2692.namemod.fa 000069K:0-72028 >>complex.orphan
samtools faidx genome.368.namemod.fa 000525K:0-12050 >>complex.orphan
samtools faidx genome.2692.namemod.fa 001006K:43663-51209 >>complex.orphan
samtools faidx genome.368.namemod.fa 000036:0-91973 >>complex.orphan
samtools faidx genome.368.namemod.fa 000036:426553-444794 >>complex.orphan
------------------------------------------------------------------------
#Manually curated the orphans list to remove duplicates and placed them in orphansamtools
cat 2692.3pr.orphan.fa 368.3pr.orphan.fa 2692.5pr.orphan.fa 368.5pr.orphan.fa complex.orphan>orphans.origname.fa
awk '/^>/{print ">" ++i"orphan"; next}{print}' < orphans.origname.fa>orphans.renamed

#concatenating merged scaffolds, and renaming
paste <(bioawk -c fastx '{print ">"$name,$seq}' 2692.5pr.merger.fa) <(bioawk -c fastx '{print $seq}' 368.5pr.merger2.fa) |sed 's/\t/\n/g'>5pr.merger.only
paste <(bioawk -c fastx '{print ">"$name,$seq}' 368.3pr.merger.fa) <(bioawk -c fastx '{print $seq}' 2692.3pr.merger2.fa) |sed 's/\t/\n/g'>3pr.merger.only

cat 3pr.merger.only 5pr.merger.only complex1 complex2 complex3 complex4 complex5>merged.fa
awk '/^>/{print ">" ++i"syntmer"; next}{print}' < merged.fa>merged.renamed


cat 5pr_exp_table 3pr_exp_table combo_table|awk '{print $1}'>2692.remove
cat 5pr_exp_table 3pr_exp_table combo_table|awk '{print $6}'>368.remove

cdbfasta genome.2692.masked.namemod.fa
cdbfasta genome.368.masked.namemod.fa

#removing scaffolds that were modified
grep -F -x -v -f 2692.remove <(bioawk -c fastx '{print $name}' genome.2692.namemod.fa) | cdbyank genome.2692.namemod.fa.cidx > syntattach.2692.namemod.fa
grep -F -x -v -f 368.remove <(bioawk -c fastx '{print $name}' genome.368.namemod.fa) | cdbyank genome.368.namemod.fa.cidx > syntattach.368.namemod.fa
cat syntattach.368.namemod.fa merged.renamed >finalsynt.368.genome.fa
cat syntattach.2692.namemod.fa orphans.renamed >finalsynt.2692.genome.fa

#making 368 scaffold name list including those merged to remove from 2692 genome
grep ">" ../genome.368.fa |sed -e's/NODE_//g' -e 's/F+//g' -e 's/F-//g' -e 's/F//g' -e 's/_length_.*//g' -e 's/,/\n>/g' -e 's/_/\n>/g' -e 's/$/K/g'>368.remove.from.2692

#removing those scaffolds
grep -F -x -v -f 368.remove.from.2692 <(bioawk -c fastx '{print $name}' syntattach.2692.namemod.fa) | cdbyank finalsynt.2692.genome.fa.cidx > finalsynt.2692.namemod.nodup.fa
###fixing error in subtractions caused by ">"
sort -V 368.remove.from.2692|sed -e 's/K//g' -e 's/>//g' |awk '{print $0"K"}'|grep -F -x -v -f - <(bioawk -c fastx '{print $name}' syntattach.2692.namemod.fa) | cdbyank finalsynt.2692.genome.fa.cidx > finalsynt.2692.namemod.nodup.fa-test
```

### Assembly information before and after synteny

```
 Final results show a slight improvement in the overall scaffold length, but scaffold number has stayed the same in 368 and gone up in 2692 (although it looks lower now because it reflects the duplicate removal.

   1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
105

---------------- Information for assembly 'finalsynt.2692.namemod.nodup.fa-test' ----------------


                                         Number of scaffolds       2075
                                     Total size of scaffolds   79152658
                                            Longest scaffold     303693
                                           Shortest scaffold        213
                                 Number of scaffolds > 1K nt       2064  99.5%
                                Number of scaffolds > 10K nt       1826  88.0%
                               Number of scaffolds > 100K nt         96   4.6%
                                 Number of scaffolds > 1M nt          0   0.0%
                                Number of scaffolds > 10M nt          0   0.0%
                                          Mean scaffold size      38146
                                        Median scaffold size      31199
                                         N50 scaffold length      48551
                                          L50 scaffold count        490
                                                 scaffold %A      31.21
                                                 scaffold %C      18.80
                                                 scaffold %G      18.77
                                                 scaffold %T      31.23
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs       0.0%
              Percentage of assembly in unscaffolded contigs     100.0%
                      Average number of contigs per scaffold        1.0
Average length of break (>25 Ns) between contigs in scaffold          0

                                           Number of contigs       2075
                              Number of contigs in scaffolds          0
                          Number of contigs not in scaffolds       2075
                                       Total size of contigs   79152658
                                              Longest contig     303693
                                             Shortest contig        213
                                   Number of contigs > 1K nt       2064  99.5%
                                  Number of contigs > 10K nt       1826  88.0%
                                 Number of contigs > 100K nt         96   4.6%
                                   Number of contigs > 1M nt          0   0.0%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size      38146
                                          Median contig size      31199
                                           N50 contig length      48551
                                            L50 contig count        490
                                                   contig %A      31.21
                                                   contig %C      18.80
                                                   contig %G      18.77
                                                   contig %T      31.23
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0



---------------- Information for assembly 'finalsynt.368.genome.fa' ----------------


                                         Number of scaffolds        368
                                     Total size of scaffolds   99552098
                                            Longest scaffold    2007323
                                           Shortest scaffold      40046
                                 Number of scaffolds > 1K nt        368 100.0%
                                Number of scaffolds > 10K nt        368 100.0%
                               Number of scaffolds > 100K nt        276  75.0%
                                 Number of scaffolds > 1M nt          7   1.9%
                                Number of scaffolds > 10M nt          0   0.0%
                                          Mean scaffold size     270522
                                        Median scaffold size     181957
                                         N50 scaffold length     403247
                                          L50 scaffold count         74
                                                 scaffold %A      31.06
                                                 scaffold %C      18.93
                                                 scaffold %G      18.95
                                                 scaffold %T      31.06
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs       0.0%
              Percentage of assembly in unscaffolded contigs     100.0%
                      Average number of contigs per scaffold        1.0
Average length of break (>25 Ns) between contigs in scaffold          0

                                           Number of contigs        368
                              Number of contigs in scaffolds          0
                          Number of contigs not in scaffolds        368
                                       Total size of contigs   99552098
                                              Longest contig    2007323
                                             Shortest contig      40046
                                   Number of contigs > 1K nt        368 100.0%
                                  Number of contigs > 10K nt        368 100.0%
                                 Number of contigs > 100K nt        276  75.0%
                                   Number of contigs > 1M nt          7   1.9%
                                  Number of contigs > 10M nt          0   0.0%
                                            Mean contig size     270522
                                          Median contig size     181957
                                           N50 contig length     403247
                                            L50 contig count         74
                                                   contig %A      31.06
                                                   contig %C      18.93
                                                   contig %G      18.95
                                                   contig %T      31.06
                                                   contig %N       0.00
                                           contig %non-ACGTN       0.00
                               Number of contig non-ACGTN nt          0

Now braker needs to be run again so that the gene positions will be correct for the newly extended scaffolds.
I have asked Arun for his help to run it again without running genemark again. (no new training) 
```
