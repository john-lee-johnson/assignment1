#!/bin/sh
#-------------------setting dir-----------------------------------------
ref=mm9_RefSeq.bed
k27ac=K27Ac_D7_Th1-W200-G200-E200-island.bedgraph
k4me3=K4m3_Th1_72h-W200-G200-E200-island.bedgraph
#---------------------finding overlapping genes---------------------------
bedtools intersect -a ${ref} -b ${k27ac} -u -bed | awk '!seen[$4] {print}
     {++seen[$4]}' > 1.bed 
wc -l 1.bed

sort 1.bed > 1sort.bed

bedtools intersect -a ${ref} -b ${k4me3} -u -bed | awk '!seen[$4] {print}
     {++seen[$4]}' > 2.bed
wc -l 2.bed
sort 2.bed > 2sort.bed

#---------------------generating TSS file-------------------------------
#TSS are found by taking the first nucleotide of the gene (start for + strand, and last for - strand)
awk 'OFS="\t"{if ($6=="+") {start=$2; end=start+1} else if ($6=="-") {start=$3-1;end=$3} print $1,start,end,$4,$5,$6}' ${ref} > tss.bed

#---------------------generating 5kbp windows around TSS--------------------
bedtools slop -i tss.bed -g mm9.genome -b 5000 > tss5000.bed

#---------------------finding TSS overlapping histone marks------------------------
bedtools intersect -a tss5000.bed -b ${k27ac} -bed > 3.bed
bedtools intersect -a tss5000.bed -b ${k4me3} -bed > 4.bed

#--------------------listing genes with TSS that have overlapping histone marks-----
bedtools intersect -a tss5000.bed -b ${k27ac} -wa -bed | awk '!seen[$4] {print}
     {++seen[$4]}' > 3gene.bed
wc -l 3gene.bed
bedtools intersect -a tss5000.bed -b ${k4me3} -wa -bed | awk '!seen[$4] {print}
     {++seen[$4]}'> 4gene.bed

#---------------------finding overlapping marks------------------------------
bedtools intersect -a 1.bed -b ${k4me3} -bed | awk '!seen[$4] {print}
     {++seen[$4]}'> 5gene.bed
bedtools intersect -a 2.bed -b ${k27ac} -bed | awk '!seen[$4] {print}
     {++seen[$4]}'> 6gene.bed

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



