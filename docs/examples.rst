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

Using a sites file to exclude unwanted sites (e.g. for 4-fold degenerate pi)
----------------

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --sites_file valid_sites.txt 