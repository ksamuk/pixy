*******************************
Generating invariant sites VCFs
*******************************

pixy facilitates the correct computation of π and dxy by allowing users to input data in a standard file format: Variant Call Format (VCF). AllSites VCFs contain invariant (AKA monomorphic) sites in addition to variant sites. Several commonly used, well-documented programs can generate AllSites VCFs. 

Here, we provide a quickstart guide for generating AllSites VCFs using two of the most widely-used tools for variant discovery and manipulation: BCFtools and GATK. Either of these tools can be run given only a set of aligned data (BAM files) and the reference sequence used to align them. The examples below were verified against BCFtools 1.23.x and GATK 4.6.x; older releases of either tool will work, but we recommend at minimum BCFtools 1.20+ and GATK 4.6+ (GATK 4.6.0.0 fixed a serious CRAM-writing bug affecting 4.3–4.5).

**BCFtools and GATK are also well-equipped to filter VCFs, and we recommend taking advantage of this to filter your data prior to analysis with pixy.**

.. warning::
    **Filter invariant and variant sites separately.** Population-genetic filters (e.g. on minor allele count, allele frequency, or quality metrics that only exist at variant sites) will silently drop the invariant sites that pixy needs. The safe pattern is to split your AllSites VCF into variant-only and invariant-only files, apply the relevant filters to each, then concatenate them back together (``bcftools concat -a``) before passing to pixy.

.. note::
    **Utilizing genomic intervals for improved runtime:** If generation of an AllSites VCF is time-consuming, we recommend parallelizing your pipeline by breaking analyses down into smaller genomic regions. In our test datasets, we ran individual chromosomes in parallel. Depending on the size of chromosomes in a dataset, it may be beneficial to break down chromosomes into smaller intervals for variant calling. Genomic intervals can be specified using the -L parameter in GATK or the -r parameter in bcftools mpileup.

Generating AllSites VCFs using BCFtools (mpileup/call)
======================================================

BCFtools mpileup can be used to produce genotype likelihoods, and this operation can be followed by bcftools call to call SNPs/INDELS. BCFtools offers a number of flexible options documented here: https://samtools.github.io/bcftools/bcftools-man.html#call

In this example, we call mpileup and pipe the output to call variants and generate an AllSites VCF::

    bcftools mpileup -Ou -f <reference.fa> -b <bamlist.txt> -r <X> --max-depth 1000 \
      | bcftools call -m -Oz -f GQ -o <output.vcf.gz>

Notes on the options selected here:

* ``-b`` points mpileup to a list of BAM files contained in a text file.
* ``-r`` specifies the genomic region. In this example, we specify the X chromosome, but mpileup provides a variety of options for specifying the region (``CHR|CHR:POS|CHR:FROM-TO|CHR:FROM-[,…]``). Alternatively, ``-R <file>`` can be used to read regions from a provided file.
* ``-Ou`` passes uncompressed BCF between ``mpileup`` and ``call``, avoiding a wasteful VCF-text round-trip and noticeably speeding up the pipeline.
* ``--max-depth`` raises mpileup's per-file cap (default 250 reads per position), which is too low for most modern coverages — set this to roughly 2–3× your expected mean depth.
* ``-m`` selects the multiallelic-and-rare-variant caller, which is the recommended default. The older consensus caller (``-c``) is retained for backward compatibility only.
* ``-Oz`` specifies bgzipped VCF output.
* ``-f GQ`` indicates that the FORMAT field GQ should be output for each sample.

.. note::
    **gVCF block output.** If storage is a concern, ``bcftools call`` also supports a gVCF-style output via ``-g/--gvcf INT``, which collapses runs of homozygous-REF calls into blocks above a per-sample depth threshold. Blocks can be expanded back to one row per site with ``bcftools convert --gvcf2vcf`` before running pixy.

Generating AllSites VCFs using GATK
===================================

GATK recommends first calling variants per-sample using HaplotypeCaller in GVCF mode (Step 1 below). Next, GenomicsDBImport consolidates information from GVCF files across samples to improve the efficiency of joint genotyping (Step 2 below). In the 3rd step, GenotypeGVCFs produces a set of jointly-called SNPs and INDELs ready for filtering and analysis. We recommend consulting the full GATK documentation at https://gatk.broadinstitute.org/.

Step 1 - HaplotypeCaller (this step should be run for each BAM file)::

    gatk --java-options "-Xmx4G" HaplotypeCaller \
        -R <reference.fa> -I <input.bam> -O <output.g.vcf.gz> \
        -ERC GVCF -L <X>

Step 2 - GenomicsDBImport::

    gatk --java-options "-Xmx4g" GenomicsDBImport \
        -V $file1 -V $file2 \
        --genomicsdb-workspace-path <allsamples_genomicsdb> \
        -L <X>

Step 3 - GenotypeGVCFs::

    gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R <reference.fa> -V gendb://<allsamples_genomicsdb> \
        --all-sites -L <X> -O <output-allsites.vcf.gz>

Notes on the options selected above:

* ``-ERC GVCF`` sets the mode for emitting reference confidence scores. GVCF format specifies condensed non-variant blocks. ``-ERC BP_RESOLUTION`` is also accepted by downstream tools and emits one record per reference position, at the cost of much larger intermediate files.
* ``-L`` specifies the genomic region. In this example, we specify the X chromosome.
* ``-V`` specifies the variant data for input. In the case of GenomicsDBImport, this is GVCF file(s). In the case of GenotypeGVCFs, the variant data for joint genotyping is provided as a GenomicsDB workspace created with GenomicsDBImport in the previous step.
* ``--all-sites`` indicates that the final VCF should include loci found to be non-variant after genotyping. **The most important parameter for our purposes.**

Important: In GATK v4.x, we have found that you must specify an interval (``-L``) when running GenotypeGVCFs. Without a designated interval, it appears to encounter missing reference confidence blocks, causing it to fail. This is true even when the problematic blocks are outside of the GenomicsDB interval being passed to it. We recommend always specifying an interval (``-L``) to avoid such issues.

.. note::
    **Skipping the all-sites decompression.** If disk space is tight, the joint-called GVCF emitted by ``GenotypeGVCFs`` (without ``--all-sites``) can be passed directly to pixy with the ``--gvcf`` flag. ``pixy`` will expand the invariant reference blocks (``ALT=<NON_REF>``, ``INFO/END``) into per-site rows at read time, producing identical results to a fully-decompressed all-sites VCF without the intermediate file. The same flag works on bcftools gVCF-style output (``bcftools call -g``); the ``bcftools convert --gvcf2vcf`` round-trip is then optional.

.. warning::
    **``--all-sites`` is not compatible with VQSR.** The 0/0 (invariant) records emitted by ``GenotypeGVCFs --all-sites`` typically carry only ``DP`` in the INFO field and lack the annotations that ``VariantRecalibrator`` relies on (``QD``, ``FS``, ``SOR``, ``MQ``, ``MQRankSum``, ``ReadPosRankSum``, ``InbreedingCoeff``). For pixy input you should hard-filter instead, and — per the warning at the top of this page — filter variant and invariant sites separately before concatenating them back together.

.. note::
    **Performance with many alt alleles.** If your callset contains many sites with a high number of alternate alleles, ``GenotypeGVCFs --all-sites`` can slow down dramatically. Splitting the genome into smaller intervals and running them in parallel is usually the cleanest fix.
