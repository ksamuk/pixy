***********************************
Step by Step Installation and Usage
***********************************

.. note::
    pixy is currently only available for Linux and macOS systems.
    
 
1. Generate a VCF with Invariant Sites and perform filtering
============================================================



If you did not already generate an 'allsites' VCF (VCF with invariant sites), `see the guide here <https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html>`_.

.. note::
    When working with whole genome data, we suggest you generate *separate invariant sites VCFs for each chromosome*. This is to prevent
    memory limitation issues down the line. This is less of an issue for reduced representation sequencing, naturally.

We recommend using the standard tools dedicated to performing such operations on VCFs: VCFtools and BCFtools (both available on bioconda). The examples below use VCFtools because its filter flags are the most legible, but BCFtools (``bcftools view``, ``bcftools filter``) can perform all of the same operations and is the more actively-maintained tool — use whichever you prefer.

Site-level filtration
---------------------
A full treatment of VCF filtering for population genetics is beyond our scope here, but for a good start see: https://speciationgenomics.github.io/filtering_vcfs/.

At minimum, we recommend filtering your VCF (a) using the GATK best practices site filtration procedure (if applicable) and (b) performing further filtering for missingness, site quality score, and mean depth (minimum and maximum).

Here is an example using VCFtools. The specific values (especially for min/max-meanDP) will vary based on your dataset:

.. code:: console

    vcftools --gzvcf my_vcf.vcf.gz \
    --remove-indels \
    --max-missing 0.8 \
    --min-meanDP 20 \
    --max-meanDP 500 \
    --recode --stdout | bgzip -c > my_filtered_vcf.vcf.gz

.. note::
    ``GenotypeGVCFs`` only assigns QUAL scores to invariant sites if there are no missing genotypes. If you wish to filter on QUAL (``--minQ``), invariant and variant sites must be filtered **separately**, with the QUAL filter applied only to variant sites (see below for details).
 
 
Optional: Population genetic filters
------------------------------------
Depending on your goal, you might also consider filtering out sites with strong HWE violations (try --hwe 0.001 with VCFtools), unusually high observed heterozygosity, or allelic depth imbalances. See this paper https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12613 for more details on these considerations. 

These last two considerations are particularly important if your study organism has high levels of gene duplication (e.g. re-diploidized after whole genome duplication as in many plant and fish species). 

If your VCF contains both variant and invariant sites (as it should at this point), applying population genetic based filters will result in the loss of your invariant sites. To avoid this, filter the invariant and variant sites separately and concatenate the two resulting files. Below is an example of one way to achieve this using VCFtools and BCFtools:
 
 .. code:: console

    #!/bin/bash
    # requires bcftools/bgzip/tabix and vcftools

    # create a filtered VCF containing only invariant sites
    vcftools --gzvcf test.vcf.gz \
    --max-maf 0 \
    [add other filters for invariant sites here] \ 
    --recode --stdout | bgzip -c > test_invariant.vcf.gz

    # create a filtered VCF containing only variant sites
    vcftools --gzvcf test.vcf.gz \
    --mac 1 \
    [add other filters for variant sites here] \ 
    --recode --stdout | bgzip -c > test_variant.vcf.gz

    # index both vcfs using tabix
    tabix test_invariant.vcf.gz
    tabix test_variant.vcf.gz

    # combine the two VCFs using bcftools concat
    bcftools concat \
    --allow-overlaps \
    test_variant.vcf.gz test_invariant.vcf.gz \
    -O z -o test_filtered.vcf.gz
 


2. Install a conda distribution
===============================
If you haven't already, install Miniforge (https://github.com/conda-forge/miniforge) or Miniconda (https://docs.anaconda.com/miniconda/). We recommend **Miniforge** — it ships with the conda-forge channel pre-configured, is fully free for any use, and has a smaller footprint than the full Anaconda distribution. The full Anaconda distribution will also work but has more restrictive licensing for commercial/large-organization use.

3. Create a New Environment
===========================
Create and activate a new conda environment for working with pixy.
``pixy`` supports Python 3.9, 3.10 and 3.11 (3.12 is not yet supported):

.. code:: console

    conda create -n "pixy" python=3.11
    conda activate pixy

4. Install pixy
===============
Install pixy via the conda-forge channel, plus ``samtools`` from bioconda for the bgzip/tabix tools you'll need to prepare your VCF.

.. code:: console

    conda install --yes -c conda-forge pixy
    conda install --yes -c bioconda samtools

(``samtools`` pulls in ``htslib`` as a dependency, so there is no need to install it separately.)

To see a list of arguments and test the pixy installation, type:

.. code:: console

    pixy --help

.. note::
    If you have trouble installing in a Python 3.11 environment, try rolling
    back to Python 3.9.


5. Create a populations file
============================
Create a populations file. This is a headerless, tab-separated file where the first column contains sample names (exactly as represented in the VCF), and the second column contains population names (these can be anything, but should be consistent!).

For example:

.. parsed-literal::
    ERS223827	BFS
    ERS223759	BFS
    ERS223750	BFS
    ERS223967	AFS
    ERS223970	AFS
    ERS223924	AFS
    ERS224300	AFS
    ERS224168	KES
    ERS224314	KES

6. Compress and Index your VCF
==============================

pixy requires its input VCF to be **bgzip-compressed** (not plain gzip) and to have a **tabix index** alongside it. If you have not already done so:

.. code:: console

    bgzip [your.file.vcf]
    tabix -p vcf [your.file.vcf.gz]

If your file was compressed with plain ``gzip`` (for example by an earlier step in your pipeline), decompress it and re-compress with ``bgzip`` — tabix will refuse to index a plain-gzip file.

7. Run pixy
===========

Run pixy! The example below uses the small simulated *D. melanogaster* VCF that ships with the docs (see :doc:`../example_data` for details on how it was generated and what to expect from it):

.. code:: console

    pixy --stats pi fst dxy \
    --vcf docs/example_data/dmel_2pop.vcf.gz \
    --populations docs/example_data/dmel_populations.txt \
    --window_size 2000 \
    --n_cores 4 \
    --chromosomes 'chr2L'

As of pixy 2.0, you can also request Watterson's θ and Tajima's *D*:

.. code:: console

    pixy --stats pi fst dxy watterson_theta tajima_d \
    --vcf docs/example_data/dmel_2pop.vcf.gz \
    --populations docs/example_data/dmel_populations.txt \
    --window_size 2000 \
    --n_cores 4 \
    --chromosomes 'chr2L'

.. note::
    By default, ``pixy`` ignores multiallelic sites and INDELs even if
    they are left in the VCF after pre-filtering. Pass
    ``--include_multiallelic_snps`` to include sites with more than two
    alleles in the calculation (new in 2.0).

8. Interpret your results
=========================

pixy writes one tab-separated file per requested statistic, named ``<output_prefix>_<stat>.txt`` (e.g. ``pixy_pi.txt``, ``pixy_fst.txt``). For a column-by-column description of each output file see :doc:`../output`, and for ready-to-run R recipes for plotting π, dxy, and Fst across the genome see :doc:`../plotting`.


9. Stay up to date
==================

You can keep pixy up to date by re-running:

.. code:: console

    conda install --yes -c conda-forge pixy
 
You can check that you have the latest version via:
 
 .. code:: console
    
    pixy --version

And comparing the version number to the one listed here: https://anaconda.org/conda-forge/pixy.
