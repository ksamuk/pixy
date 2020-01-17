************
Step by Step Installation and Usage
************

.. note::
    pixy is currently only available for Linux and macOS systems.
    
 
1. Generate a VCF with Invariant Sites
======
If you did not already generate an 'allsites' VCF (VCF with invariant sites), see the guide here: https://pixy.readthedocs.io/en/latest/invar/allsitesQuickstart.html

.. note::
    When working with whole genome data, we suggest you generate *separate invariant sites VCFs for each chromosome*. This is to prevent
    memory limitation issues down the line. This is less of an issue for reduced representation sequencing, naturally.

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
Create a populations file. This is a headerless, tab-separated file where the first column are sample names (exactly as represented in the VCF), and the second column are population names (these can be anything).

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

    pixy --interval_start 1 \
    --interval_end 100000 \
    --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --zarr_path data/vcf/ag1000/chrX_36Ag_allsites \
    --chromosome X \
    --window_size 10000 \
    --populations data/vcf/ag1000/Ag1000_sampleIDs_popfile.txt \
    --variant_filter_expression DP>=10,GQ>=20,RGQ>=20 \
    --invariant_filter_expression DP>=10,RGQ>=20 \
    --outfile_prefix output/pixy_out

7. Profit
======

Parse the output files and enjoy your unbiased estimates of pi and dxy!
