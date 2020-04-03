************
Arguments
************

Below is a list of arguments that pixy accepts.

--help       Print the help message. 
--version       Print the pixy version number.
--stats         Which statistics to calculate from the VCF 
               (pi, dxy, and/or fst, separated by spaces)
--vcf           Path to the input VCF.
--zarr_path            Folder in which to build the Zarr array(s).
--reuse_zarr           Use existing Zarr array(s) (saves time if re-running). [yes,no] 
--populations            Path to the populations file. See quick start for format.
--window_size           Window size in base pairs over which to calculate pi/dxy.
--chromosomes            A single-quoted, comma separated list of chromosome(s) (e.g. 'X,1,2'). Defaults to all chromosomes in the VCF.
--interval_start            The start of the interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome.
--interval_end            The end of the interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome.
--variant_filter_expression            A comma separated list of filters contained in single quotes.
                                       (e.g. 'DP>=10,GQ>=20') to apply to SNPs.
--invariant_filter_expression          A comma separated list of filters contained in single quotes.
                                       (e.g. 'DP>=10,RGQ>=20') to apply to invariant sites.
--outfile_prefix            Path and prefix for the output file. Output files will be named like: 
                            path/to/outfile_pi_[popname].txt
--bypass_filtration            Bypass all variant filtration (for data lacking FORMAT annotations, 
                                use with extreme caution!)

An example:

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --zarr_path data/zarr \
    --window_size 10000 \
    --populations Ag1000_sampleIDs_popfile.txt \
    --variant_filter_expression 'DP>=10,GQ>=20,RGQ>=20' \
    --invariant_filter_expression 'DP>=10,RGQ>=20' \
    --outfile_prefix output/pixy_out
