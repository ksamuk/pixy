#!/bin/bash

# injects missing invariant sites into a msprime vcf
# ks may 2019

# takes a msprime vcf as input
vcf=$1
outfile=$2

# parse the length of the sequence from the header
seq_length=$(grep length $vcf | sed 's/.*length=//g' | sed 's/>//g')

# count the number of samples in the file
n_samples=$(grep "#CHROM" $vcf  | awk '{print NF; exit}')

# there are 9 non-genotype columns
let n_samples=n_samples-9

echo "injecting invariant sites into $vcf..."
echo "$n_samples samples"
echo "$seq_length sites"

# all the possible sites
sites=$(seq 1 $seq_length)
echo "$sites" | sort -V > all_sites.tmp

# sites with variants

var_sites=$(grep -v "#" $vcf | awk '{print $2}')
echo "$var_sites" | sort -V > var_sites.tmp

# find sites that need invar
diff --new-line-format="" --unchanged-line-format="" all_sites.tmp var_sites.tmp | sort -V > invar_sites.tmp

invar_sites=$(cat invar_sites.tmp)


# create a VCF with all invariant sites

# the start of a blank row
row=".\tA\tT\t.\tPASS\t.\tGT"

while read site
do

	gt=$(printf '0|0\t%.0s' $(eval echo "{1..$n_samples}"))

	newline="1\t$site\t$row\t$gt"
	echo $newline >> vcf_blank_spaces.vcf

done < invar_sites.tmp

rm invar_sites.tmp
rm var_sites.tmp
rm all_sites.tmp

gsed 's/\s/\t/g' vcf_blank_spaces.vcf | gsed 's/^[ \t]*//;s/[ \t]*$//' > vcf_blank.vcf

grep "#" $vcf > vcf_header.vcf

grep -v "#" $vcf > vcf_variants.vcf

cat vcf_blank.vcf vcf_variants.vcf > test_tmp.vcf

cat test_tmp.vcf | sort -k1V,1 -k2n,2 > test.vcf

cat vcf_header.vcf test.vcf | grep . | gzip -c > $outfile 

rm vcf_header.vcf vcf_blank_spaces.vcf vcf_variants.vcf test_tmp.vcf test.vcf

echo "wrote to $outfile"

