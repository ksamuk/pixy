**************
Usage Examples
**************

Basic usage
-----------

Compute π, F\ :sub:`ST`, and d\ :sub:`xy` in 10 kb windows across the
whole VCF:

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --window_size 10000 \
    --n_cores 2

Including the new 2.0 statistics
--------------------------------

Watterson's θ and Tajima's *D* are also computed in windows. They use
the same populations file as π:

.. code:: console

    pixy --stats pi watterson_theta tajima_d \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --window_size 10000 \
    --n_cores 4

Focusing on a specific genomic interval
---------------------------------------

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --chromosomes 'X' \
    --interval_start 15000 \
    --interval_end 400000 \
    --window_size 10000

Using a .BED file to define windows
-----------------------------------

When defined with a BED file, windows can be any size and non-uniform.

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --bed_file genomic_windows.bed

Using a sites file to limit calculations to specific sites
----------------------------------------------------------

Sites files define individual sites that will be included in the
calculations. Useful for obtaining 4-fold degenerate π, or statistics
for particular classes of genomic elements (genes, introns, etc.).

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --sites_file valid_sites.txt

Including multiallelic SNPs
---------------------------

By default ``pixy`` considers only biallelic sites. Pass
``--include_multiallelic_snps`` to include sites with more than two
alleles in the calculation (added in 2.0):

.. code:: console

    pixy --stats pi dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --window_size 10000 \
    --include_multiallelic_snps

Switching the F\ :sub:`ST` estimator
------------------------------------

``pixy`` defaults to Weir & Cockerham's F\ :sub:`ST`. Pass
``--fst_type hudson`` to use the Hudson (1992) / Bhatia *et al.* (2013)
estimator instead:

.. code:: console

    pixy --stats fst \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --window_size 10000 \
    --fst_type hudson

Site-level estimates
--------------------

.. note::

    Site-level estimates are much slower to calculate than windowed
    estimates. Use ``--window_size 1`` with a narrow interval for
    targeted lookups rather than genome-wide site-level scans.

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --window_size 1 \
    --chromosomes 'X' \
    --interval_start 15000 \
    --interval_end 400000
