************
Usage Examples
************

Basic usage with genotype filtering
----------------

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --zarr_path data/zarr \
    --window_size 10000 \
    --populations Ag1000_sampleIDs_popfile.txt \
    --variant_filter_expression 'DP>=10,GQ>=20,RGQ>=20' \
    --invariant_filter_expression 'DP>=10,RGQ>=20' \
    --outfile_prefix output/pixy_out

Using a pre-filtered VCF (or simulated data with no FORMAT fields)
----------------

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --zarr_path data/zarr \
    --window_size 10000 \
    --populations Ag1000_sampleIDs_popfile.txt \
    --bypass_filtration yes \
    --outfile_prefix output/pixy_out
    
Focusing on a specific genomic interval
----------------

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --zarr_path data/zarr \
    --chromosome X \
    --interval_start 15000 \
    --interval_end 400000 \
    --window_size 10000 \
    --populations Ag1000_sampleIDs_popfile.txt \
    --bypass_filtration yes \
    --outfile_prefix output/pixy_out

