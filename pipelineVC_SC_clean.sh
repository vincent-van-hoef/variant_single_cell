#!/bin/bash
# This is a script to count alleles in a certain region over multiple bam files. You need to supply a reference file and a BED file of the region of interest. The input are BAM files, the output is a file with per line one allele (so multiple lines for one location are possible) and the number of reads covering that allele per BSAM file. Output of this can be visualized with plot.R (and Panel.R)
function reshape {
while read line; do
	alleles=$(echo $line | awk '{n=split($3, a, ",")}END{print n}')
	for j in $(seq 1 $alleles); do
		echo $line | awk '{for(i=3; i<=NF; i++){split($i, a, ","); $i=a['$j']}} 1'
	done
done < $1
}

REF=~/projects/Marine_Lab/data/ref/hg19/hg19.fa
BED=~/projects/Marine_Lab/data/ref/trusight_tumor_reportedregions.bed

# Count reads per allele and extract relevant information, ignore INDELS with -I
samtools mpileup -ARBQ0 -I -d 1000000000 -f $REF -l $BED -t INFO/AD,AD -vu -b $1 | bcftools query -H -f '%CHROM\t%POS\t%REF,%ALT\t[%AD\t]\n' > "${1}_tmp"
# wide to long data transformation
reshape "${1}_tmp" > "${1}_tmpLong" 
awk 'NR==1 {$1=""}1' "${1}_tmpLong" | tr " " "\t" | awk '{if (NR==1) print $0; else print "\t"$0}'| cut -f1 --complement > "${1}.out"
rm "${1}_tmp"
rm "${1}_tmpLong"


#plot.R
#panel.R
