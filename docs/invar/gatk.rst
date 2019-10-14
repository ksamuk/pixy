************
Quickstart Guide for Generating PIXY Input
************

pixy facilitates the correct computation of π and dxy by allowing users to input data in a standard file format: Variant Call Format (VCF). AllSites VCFs contain invariant (AKA monomorphic) sites in addition to variant sites. Several commonly used, well-documented programs can generate AllSites VCFs. Here, we provide a quickstart guide for generating AllSites VCFs using two of the most widely-used tools for variant discovery and manipulation: BCFtools and GATK. Either of these tools can be run given only a set of aligned data (BAM files) and the reference sequence used to align them.

.. Utilizing genomic intervals for improved runtime::
    If generation of an AllSites VCF is time-consuming, we recommend parallelizing      your pipeline by breaking analyses down into smaller genomic regions. In our test datasets, we ran individual chromosomes in            parallel. Depending on the size of chromosomes in a dataset, it may be beneficial to break down chromosomes into smaller intervals      for variant calling. Genomic intervals can be specified using the -L parameter in GATK or the -r parameter in bcftools mpileup.

Generating AllSites VCFs using BCFtools (mpileup/call)
===================

BCFtools mpileup can be used to produce genotype likelihoods, and this operation can be followed by bcftools call to call SNPs/INDELS. BCFtools offers a number of flexible options documented here: https://samtools.github.io/bcftools/bcftools-man.html#call

In this example, we call mpileup and pipe the output to call variants and generate and AllSites VCF::

    bcftools mpileup -f <reference.fa> -b <bamlist.txt> -r <X> | bcftools call -m -Oz -f GQ -o <output>

Notes on the options selected here:

* b points mpileup to a list of BAM files contained in a text file.
* r specifies the genomic region. In this example, we specify the X chromosome, but mpileup provides a variety of options for specifying the region (CHR|CHR:POS|CHR:FROM-TO|CHR:FROM-[,…]). Alternatively, -R <file> can be used to read regions from a provided file.
* Oz specifies compressed VCF output
* f GQ indicates that the FORMAT field GQ should be output for each sample

Generating AllSites VCFs using GATK
===================

GATK recommends first calling variants per-sample using HaplotypeCaller in GVCF mode (Step 1 below). Next, GenomicsDBImport consolidates information from GVCF files across samples to improve the efficiency joint genotyping (Step 2 below). In the 3rd step, GenotypeGVCFs produces a set of jointly-called SNPs and INDELS ready for filtering and analysis. We recommend consulting the full GATK documentation found here: https://software.broadinstitute.org/gatk/

Step1 - HaplotypeCaller (this step should be run for each BAM file):: console

    gatk-4.0.7.0/gatk --java-options "-Xmx4G" HaplotypeCaller \
    -R <reference.fa> -I <input.bam> -O <output.g.vcf> -ERC GVCF -L <X>

Step2 - GenomicsDBImport:: console

    gatk-4.0.7.0/gatk --java-options "-Xmx4g" GenomicsDBImport \
    -V $file1 -V $file2 --genomicsdb-workspace-path <allsamples_genomicsdb> \
    -L <X>

Step3 - GenotypeGVCFs:: console

    gatk-4.0.7.0/gatk --java-options "-Xmx4g" GenotypeGVCFs \
    -R <reference.fa> -V gendb://<allsamples_genomicsdb> \
    -all-sites -L <X> -O <output-allsites.vcf.gz>

Notes on the options selected above:

* ERC GVCF sets the mode for emitting reference confidence scores. GVCF format specifies condensed non-variant blocks.
* L specifies the genomic region. In this example, we specify the X chromosome
* V specifies the variant data for input. In the case of GenomicsDBImport, this is GVCF file(s). In the case of GenotypeGVCFs, the variant data for joint genotyping is provided as a GenomicsDB workspace created with GenomicsDBImport in the previous step.
* all-sites indicates that the final VCF should include loci found to be non-variant after genotyping. The most important parameter for our purposes.

Important: In GATK v4.x, we have found that you must specify an interval (-L) when running GenotypeGVCFs. Without a designated interval, it appears to encounter missing reference confidence blocks, causing it to fail. This is true even when the problematic blocks are outside of the GenomicsDB interval being passed to it. We recommend always specifying an interval (-L) to avoid such issues.

.. note::
    In the example above, we use GATK v4, but AllSites VCFs can also be easily generated in GATK v3.x by running GenotypeGVCFs with the “-allSites” parameter. (Note the slightly different syntax from “-all-sites” in GATK v4).
