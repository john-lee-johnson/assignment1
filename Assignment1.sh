#!/bin/sh
#-------------------setting dir-----------------------------------------
ref=mm9_RefSeq.bed
k27ac=K27Ac_D7_Th1-W200-G200-E200-island.bedgraph
k4me3=K4m3_Th1_72h-W200-G200-E200-island.bedgraph
#---------------------finding overlapping genes---------------------------
bedtools intersect -a ${ref} -b ${k27ac} -u -bed > 1.bed 
cut -f4 -d$'\t' 1.bed | uniq -u > 1.txt 
wc -l 1.txt
bedtools intersect -a ${ref} -b ${k4me3} -u -bed > 2.bed
cut -f4 -d$'\t' 2.bed | uniq -u > 2.txt 
wc -l 2.txt


#---------------------generating TSS file-------------------------------
#TSS are found by taking the first nucleotide of the gene (start for + strand, and last for - strand)
awk 'OFS="\t"{if ($6=="+") {start=$2; end=start+1} else if ($6=="-") {start=$3-1;end=$3} print $1,start,end,$4,$5,$6}' ${ref} > tss.bed

#---------------------generating 5kbp windows around TSS--------------------
bedtools slop -i tss.bed -g mm9.genome -b 5000 > tss5000.bed

#---------------------finding TSS overlapping histone marks------------------------
bedtools intersect -a tss5000.bed -b ${k27ac} -bed > 3.bed
bedtools intersect -a tss5000.bed -b ${k4me3} -bed > 4.bed

#--------------------listing genes with TSS that have overlapping histone marks-----
cut -f4 -d$'\t' 3.bed | uniq -u > 3gene.txt 
wc -l 3gene.txt
cut -f4 -d$'\t' 4.bed | uniq -u > 4gene.txt 
wc -l 4gene.txt

#---------------------finding overlapping marks------------------------------
bedtools intersect -a 1.bed -b ${k4me3} -bed > 5.bed
bedtools intersect -a 2.bed -b ${k27ac} -bed > 6.bed

#--------------------listing genes with that have overlapping histone marks-----
cut -f4 -d$'\t' 5.bed | uniq -u > 5gene.txt 
wc -l 5gene.txt
cut -f4 -d$'\t' 6.bed | uniq -u > 6gene.txt 
wc -l 6gene.txt

#---------------------sorting and merging-------------------------------
#All overlapping peaks are merged 
sort -k1,1 -k2,2n 3.bed > 3sort.bed
sort -k1,1 -k2,2n 4.bed > 4sort.bed
sort -k1,1 -k2,2n 5.bed > 5sort.bed
sort -k1,1 -k2,2n 6.bed > 6sort.bed
bedtools merge -i 3sort.bed > 3merge.bed
bedtools merge -i 4sort.bed > 4merge.bed
bedtools merge -i 5sort.bed > 5merge.bed
bedtools merge -i 6sort.bed > 6merge.bed



