************
Usage Examples
************

Basic usage 
----------------

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --window_size 10000 \
    --n_cores 2

Focusing on a specific genomic interval
----------------

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --chromosomes 'X' \
    --interval_start 15000 \
    --interval_end 400000 \
    --window_size 10000

Using a .BED file to define windows (can be any size, and non-uniform).
----------------

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --bed_file genomic_windows.bed

Using a sites file to limit calculations to specific sites 
----------------

With a list of 4-fold degenerate sites, this could be used to obtain 4-fold degenerate pi.

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --sites_file valid_sites.txt 

Extracting site-level estimates of pi, fst and dxy for a region
----------------

Note: site level estimates will be much slower to calculate than windowed estimates.

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --window_size 1 \
    --chromosomes 'X' \
    --interval_start 15000 \
    --interval_end 400000 