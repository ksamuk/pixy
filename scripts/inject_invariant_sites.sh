#!/bin/bash

# injects missing invariant sites into a msprime vcf
# ks may 2019

vcfdir=$1
outdir=$2

mkdir -p inject_tmp
mkdir -p $outdir

# list the vcf files
ls $vcfdir/*.vcf > inject_tmp/vcf_files.tmp

# exemplar vcf file for building a fake invariants sites vcf
vcfex=$(cat inject_tmp/vcf_files.tmp | head -n 1) 

# parse the length of the sequence from the header
seq_length=$(grep length $vcfex | sed 's/.*length=//g' | sed 's/>//g')

# count the number of samples in the file
n_samples=$(grep "#CHROM" $vcfex  | awk '{print NF; exit}')

# there are 9 non-genotype columns
let n_samples=n_samples-9

# all the possible sites
sites=$(seq 1 $seq_length)
echo "$sites" | sort -V > inject_tmp/all_sites.tmp

while read vcf
do

# takes a msprime vcf as input

outfile=$(echo $vcf | sed 's/.vcf/_invar.vcf.gz/g')

echo "injecting invariant sites into $vcf..."
echo "$n_samples samples"
echo "$seq_length sites"

# sites with variants

var_sites=$(grep -v "#" $vcf | awk '{print $2}')
echo "$var_sites" | sort -V > inject_tmp/var_sites.tmp

# find sites that need invar
diff --new-line-format="" --unchanged-line-format="" inject_tmp/all_sites.tmp inject_tmp/var_sites.tmp | sort -V > inject_tmp/invar_sites.tmp

invar_sites=$(cat inject_tmp/invar_sites.tmp)

# create a VCF with all invariant sites

# the start of a blank row
row=".\tA\tT\t.\tPASS\t.\tGT"

while read site
do

	gt=$(printf '0|0\t%.0s' $(eval echo "{1..$n_samples}"))

	newline="1\t$site\t$row\t$gt"
	echo $newline >> inject_tmp/vcf_blank_spaces.vcf

done < inject_tmp/invar_sites.tmp

rm inject_tmp/invar_sites.tmp
rm inject_tmp/var_sites.tmp

gsed 's/\s/\t/g' inject_tmp/vcf_blank_spaces.vcf | gsed 's/^[ \t]*//;s/[ \t]*$//' > inject_tmp/vcf_blank.vcf

grep "#" $vcf > inject_tmp/vcf_header.vcf

grep -v "#" $vcf > inject_tmp/vcf_variants.vcf

cat inject_tmp/vcf_blank.vcf inject_tmp/vcf_variants.vcf > inject_tmp/test_tmp.vcf

cat inject_tmp/test_tmp.vcf | sort -k1V,1 -k2n,2 > inject_tmp/test.vcf

cat inject_tmp/vcf_header.vcf inject_tmp/test.vcf | grep . | gzip -c > $outfile 

rm inject_tmp/vcf_header.vcf inject_tmp/vcf_blank_spaces.vcf inject_tmp/vcf_variants.vcf inject_tmp/test_tmp.vcf inject_tmp/test.vcf

echo "wrote to $outfile"

done < inject_tmp/vcf_files.tmp

#clean up
cp $vcfdir/*_invar.vcf.gz $outdir

rm $vcfdir/*_invar.vcf.gz 

rm -r inject_tmp


