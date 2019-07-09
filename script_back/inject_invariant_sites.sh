#!/bin/bash

# injects missing invariant sites into a msprime vcf
# ks may 2019

# takes a msprime vcf as input
vcf="data/msprime/test_vcf_var_only.vcf"
outfile="vcf_out.txt"

# parse the length of the sequence from the header
seq_length=$(grep length data/msprime/test_vcf_var_only.vcf | sed 's/.*length=//g' | sed 's/>//g')

# count the number of samples in the file
n_samples=$(grep "#CHROM" data/msprime/test_vcf_var_only.vcf | awk '{print NF; exit}')

# there are 9 non-genotype columns
let n_samples=n_samples-9

sites=$(seq 1 $seq_length)

echo "$n_samples samples"
echo "$seq_length sites"

# create a VCF with all invariant sites

row=".\tA\tT\t.\tPASS\t.\tGT"

rm vcf_blank.vcf

for i in $(seq 1 $seq_length)
do

	gt=$(printf '0|0\t%.0s' $(eval echo "{1..$n_samples}"))

	newline="1\t$i\t$row\t$gt"
	echo $newline >> vcf_blank_spaces.vcf

done

# fix the delimiters in the blank file

gsed 's/\s/\t/g' vcf_blank_spaces.vcf | gsed 's/^[ \t]*//;s/[ \t]*$//' > vcf_blank.vcf
rm vcf_blank_spaces.vcf


# create header/headless msprime vcf
rm vcf_out.vcf
grep "#" $vcf > vcf_out.vcf

for i in $(seq 1 $seq_length)
do

	targetsite="1\t$i\t"
	varexists=$(grep $targetsite $vcf | wc -c)
	
	#echo "$i $varexists"
	
	if [ "$varexists" -gt 1 ]; then
		#echo "variant is present"
		grep $targetsite $vcf >> vcf_out.vcf
	else
		#echo "variant is missing"
		grep $targetsite vcf_blank.vcf >> vcf_out.vcf
	fi
done
