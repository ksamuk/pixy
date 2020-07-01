************
Arguments
************

Below is a list of required and optional arguments that pixy accepts. 

--help       Print the help message. 
--version       Print the pixy version number.
--stats         **Required.** Which statistics to calculate from the VCF 
               (pi, dxy, and/or fst, separated by spaces)
--vcf           **Required.** Path to the input VCF.
--zarr_path            **Required.** Folder in which to build the Zarr array(s).
--reuse_zarr           Use existing Zarr array(s) (saves time if re-running). [yes,no] 
--populations            **Required.** Path to the populations file. See quick start for format.
--chromosomes            A single-quoted, comma separated list of chromosome(s) (e.g. 'X,1,2'). Defaults to all chromosomes in the VCF.
--window_size           Window size in base pairs over which to calculate pi/dxy. Defaults to the whole chromosome.
--interval_start            The start of a specific interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome.
--interval_end            The end of a specific interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome.

--variant_filter_expression           **Required.** A comma separated list of filters contained in single quotes.
                                       (e.g. 'DP>=10,GQ>=20') to apply to SNPs.
--invariant_filter_expression          **Required.** A comma separated list of filters contained in single quotes.
                                       (e.g. 'DP>=10,RGQ>=20') to apply to invariant sites.
--bypass_filtration            Bypass all variant filtration (for data lacking FORMAT annotations, 
                                use with extreme caution!) [yes,no]
--bypass_invariant_check            Bypass the check for invariant sites. Use with caution!
--fst_maf_filter       Minor allele frequency filter for FST calculations, with value 0.0-1.0. Sites with MAF less than this value will be excluded.
--outfile_prefix            **Required.** Path and prefix for the output file. Output files will be named like: 
                            path/to/outfile_pi_[popname].txt
--n_cores           Number of CPU cores to use for pi and dxy calculations.

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
