************
Step by Step Installation and Usage
************

.. note::
    pixy is currently only available for Linux and macOS systems.
    
 
1. Generate a VCF with Invariant Sites and perform filtering
======



If you did not already generate an 'allsites' VCF (VCF with invariant sites), `see the guide here <https://pixy.readthedocs.io/en/1.0.0.beta1/generating_invar/generating_invar.html>`_.

.. note::
    When working with whole genome data, we suggest you generate *separate invariant sites VCFs for each chromosome*. This is to prevent
    memory limitation issues down the line. This is less of an issue for reduced representation sequencing, naturally.

We recommend using the standard tools dedicated to performing such operations on VCFs: VCFtools and BCFtools (both available on bioconda).

Site-level filtration
------------------------
A full treatment of VCF filtering for population genetics is beyond our scope here, but for a good start see: https://speciationgenomics.github.io/filtering_vcfs/.

At minumum, we reccomend filtering your VCF (a) using the GATK best practices site filtration procedure (if applicable) and (b) performing further filtering for missingness, site quality score, and mean depth (minimum and maxmium). 

Here is an example using VCFtools. The specific values (especially for min/max-meanDP) will vary based on your dataset: 

.. code:: console

    vcftools --gzvcf my_vcf.vcf.gz \
    --remove-indels \
    --max-missing 0.8 \
    --minQ 30 \
    --min-meanDP 20 \
    --max-meanDP 500 \
    --recode --stdout | gzip -c > my_filtered_vcf.vcf.gz
 
 
Optional: Population genetic filters
------------------------
Depending on your goal, you might also consider filtering out sites with strong HWE violations (try --hwe 0.001 with VCFtools), unusually high observed heterozygosity, or allelic depth imbalances. See this paper https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12613 for more details on these considerations. 

These last two considerations are particularly important if your study organism has high levels of gene duplication (e.g. re-diploidized after whole genome duplication as in many plant and fish species). 

If your VCF contains both variant and invariant sites (as it should at this point), applying population genetic based filters will result in the loss of your invariant sites. To avoid this, filter the invariant and variant sites separately and concatenate the two resulting files. Below is an example of one way to achieve this using VCFtool and BCFtools:
 
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
 


2. Install Anaconda
======
If you haven't already, install Anaconda https://docs.anaconda.com/anaconda/install/ 

3. Create a New Environment
======
Create and activate a new conda environment for working with pixy:

.. code:: console

    conda create --name pixy
    conda activate pixy

4. Install pixy
======
Install pixy via the conda-forge channel. Also install the required htslib package from bioconda.

.. code:: console

    conda install --yes -c conda-forge pixy
    conda install --yes -c bioconda htslib

To see a list of arguments and test the pixy installation, type:

.. code:: console

    pixy --help


5. Create a populations file
======
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

    
6. Run pixy
======

Run pixy! An example is shown below.

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations data/vcf/ag1000/Ag1000_sampleIDs_popfile.txt \
    --window_size 10000 \
    --n_cores 4 \
    --chromosomes 'X' 

.. note::
    pixy ignores non-biallelic sites and INDELs, even if they are left in the VCF after pre-filtering. 

7. Profit
======

Parse the output files and enjoy your unbiased estimates of pi and dxy!


8. Stay up to date
======

You can keep pixy up to date by re-running:

.. code:: console

    conda install --yes -c conda-forge pixy
 
You can check that you have the latest version via:
 
 .. code:: console
    
    pixy --version

And comparing the version number to the one listed here: https://anaconda.org/conda-forge/pixy.
