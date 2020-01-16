************
Arguments
************

Below is a list of arguments that pixy accepts.

* ``--version`` Print the pixy version number.

* ``--stats [fst,dxy,pi]`` Which statistics to calculate from the VCF (pi, dxy, and/or fst, separated by spaces)', required=True).

* ``--vcf`` Path to the input VCF.

* ``--zarr_path`` Folder in which to build the Zarr array.

* ``--regenerate_zarr [yes, no]`` Force regeneration of the Zarr array.

* ``--populations`` Path to the populations file.

* ``--window_size`` Window size in base pairs over which to calculate pi/dxy.

* ``--chromosome`` Target chromosome (as annotated in the CHROM field).

* ``--interval_start`` The start of the interval over which to calculate pi/dxy.

* ``--interval_end`` The end of the interval over which to calculate pi/dxy.

* ``--variant_filter_expression`` A comma separated list of filters (e.g. DP>=10,GQ>=20) to apply to SNPs.

* ``--invariant_filter_expression`` A comma separated list of filters (e.g. DP>=10,RGQ>=20) to apply to invariant sites.

* ``--outfile_prefix`` Path and prefix for the output file, e.g. path/to/outfile.

* ``--bypass_filtration`` [yes,no] Bypass all variant filtration (for data lacking FORMAT annotations, use with extreme .caution)

An example:

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
