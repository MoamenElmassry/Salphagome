
## Run Magic-BLAST to identify phages in Salmonella WGS

``
magicblast -db $BLAST_DB -sra $i -no_unaligned -num_threads $THREADS -out $OUT_DIR/$i.sam
``

## Calculate coverage of the phages

```
for i in *.sam; do 
samtools view -Sb $i > $i.bam; done
for i in *.sam.bam; do 
samtools sort $i -o $i.sorted.bam; done
for i in `ls *.bam`; do 
samtools index $i $i.bai; done
for i in *.bam; do 
bamCoverage --binSize 100 -b $i --outFileFormat bedgraph -o $i.bedgraph; done
for f in `ls *.bedgraph`; do
echo ""
echo $f
echo "---"
awk -F\t 'NR>1{if ($4 < 3) {arrlow[$1]++} else {arrhi[$1]++}} END {for (a in arrlow) { if(length(arrlow[a]) == 0){arrlow[a] = 0} if(length(arrhi[a]) == 0){arrhi[a] = 0} print a, 100*(arrhi[a])/(arrlow[a] + arrhi[a]),"%" }}' $f; done
```

## Calculate percentage of WGS that carries each phage
```
library(plyr)
library(reshape2)
library(reshape)
data=read.csv("coverage.csv",header=T,check.names = F)

data$cov[data$cov < 0.5] <- 0
data$cov[data$cov >= 0.5] <- 1
samples_count <- count(data, c('seq','strain'))

cov_count=aggregate(cov~seq+strain, data, sum)

merged=merge(samples_count, cov_count, all.x = TRUE)
merged$percentage=merged$cov*100/merged$freq
casted=cast(merged, seq~strain,value="percentage")
write.csv(casted,"new calculation.csv")
```
