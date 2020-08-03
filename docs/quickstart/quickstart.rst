************
Step by Step Installation and Usage
************

.. note::
    pixy is currently only available for Linux and macOS systems.
    
 
1. Generate a VCF with Invariant Sites and perform any pre-filtering
======
If you did not already generate an 'allsites' VCF (VCF with invariant sites), see the guide here: https://pixy.readthedocs.io/en/latest/invar/allsitesQuickstart.html 

.. note::
    When working with whole genome data, we suggest you generate *separate invariant sites VCFs for each chromosome*. This is to prevent
    memory limitation issues down the line. This is less of an issue for reduced representation sequencing, naturally.

Note that while pixy provides limited genotype-level filtering, VCF filtering can be complex. Thus, we recommend applying one of the several tools dedicated to performing such operations on VCFs, such as BCFtools: http://samtools.github.io/bcftools/bcftools.html

Site-level filtration
------------------------
The goal of site-level filtration is to remove sites that show evidence of sequencing errors, duplication, or other issues with mapping. You will likely want to a least apply site-level filters based on the number of missing individuals, site quality score, and depth (minimum and maxmium). This guide https://speciationgenomics.github.io/filtering_vcfs/ is a good starting place. Note that if you want to filter on minor allele frequency, you will need to be mindful that your invariant sites (which all have MAF = 0) are not removed as part of this process (e.g. by performing two different filtering passes and joining the results).

Here is an example using VCFtools. The specific values (especially for min/max-meanDP) will vary based on your dataset: 

.. code:: console

    vcftools --gzvcf my_vcf.vcf.gz \
    --remove-indels \
    --max-missing 0.8 \
    --minQ 30 \
    --min-meanDP 10 \
    --max-meanDP 100 \
    --recode --stdout | gzip -c > my_filtered_vcf.vcf.gz

You might also want to filter out sites with strong HWE violations (try --hwe 0.001 with VCFtools), unusually high observed heterozygosity, or allelic depth imbalances. See this paper https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12613 for more details on these considerations. These last two considerations are particularly important if your study organism has high levels of paralogy (e.g. re-diploidized after whole genome duplication as in many plant and fish species). Again, be mindful that your invariant sites will also be affected by these filters.


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
Install pixy via the conda-forge channel. 

.. code:: console

    conda install --yes -c conda-forge pixy

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
    --zarr_path data/zarr/ag1000 \
    --chromosomes 'X' \
    --window_size 10000 \
    --populations data/vcf/ag1000/Ag1000_sampleIDs_popfile.txt \
    --variant_filter_expression 'DP>=10,GQ>=20,RGQ>=20' \
    --invariant_filter_expression 'DP>=10,RGQ>=20' \
    --outfile_prefix output/pixy_out

If your VCF is pre-filtered, you can also bypass genotype filtration:

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --zarr_path data/zarr/ag1000 \
    --chromosomes 'X' \
    --window_size 10000 \
    --populations data/vcf/ag1000/Ag1000_sampleIDs_popfile.txt \
    --bypass_filtration yes \
    --outfile_prefix output/pixy_out

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
