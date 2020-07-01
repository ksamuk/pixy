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

Note that while pixy provides some limited filtering expressions, VCF filtering can be complex. Thus, we recommend applying one of the several tools dedicated to performing such operations on VCFs, such as BCFtools: http://samtools.github.io/bcftools/bcftools.html

.. note::
    You will likely want to a least apply site filters based on the number of missing individuals, minor allele frequency, site quality score, and depth. This guide https://speciationgenomics.github.io/filtering_vcfs/ is a good starting place. You might also want to remove sites with HWE violations, unusually high osberved heterozygosity, or allelic depth imbalances. See this paper https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12613 for more details on these considerations.



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

.. note::
    pixy ignores non-biallelic sites. If you want to compute pi with polyallelic sites and/or INDELs, please let us know! 

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
