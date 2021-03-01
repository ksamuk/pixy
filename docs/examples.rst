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

Using a .BED file to define windows
----------------

When defined with a BED file, windows can be any size, and non-uniform.

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --bed_file genomic_windows.bed

Using a sites file to limit calculations to specific sites 
----------------

Sites files define individual sites that will be included in the calculations. This can be used to obtain 4-fold degenerate pi, or statistics for particular classes of genomic elements (genes, introns, etc.).

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --sites_file valid_sites.txt 

Site-level estimates
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