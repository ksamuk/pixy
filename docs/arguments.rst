************
Arguments
************

Option Lists
------------

For listing command-line options:

-a            command-line option "a"
-b file       options can have arguments
              and long descriptions
--long        options can be long also
--input=file  long options can also have
              arguments

--very-long-option
              The description can also start on the next line.

              The description may contain multiple body elements,
              regardless of where it starts.

-x, -y, -z    Multiple options are an "option group".
-v, --verbose  Commonly-seen: short & long options.
-1 file, --one=file, --two file
              Multiple options with arguments.
/V            DOS/VMS-style options too

There must be at least two spaces between the option and the description.
    
Below is a list of arguments that pixy accepts.

<dt>--version</dt> Print the pixy version number.

* ``--stats [fst,dxy,pi]`` Which statistics to calculate from the VCF (pi, dxy, and/or fst, separated by spaces)', required=True).

* ``--vcf [path/to/vcf]`` Path to the input VCF.

* ``--zarr_path [path/to/zarr/folder]`` Folder in which to build the Zarr array.

* ``--regenerate_zarr [yes, no]`` Force regeneration of the Zarr array.

* ``--populations [population_file.txt]`` Path to the populations file. See quick start for format.

* ``--window_size [integer]`` Window size in base pairs over which to calculate pi/dxy.

* ``--chromosome [string]`` Target chromosome (precisely as annotated in the CHROM field).

* ``--interval_start [integer]`` The start of the interval over which to calculate pi/dxy.

* ``--interval_end [integer]`` The end of the interval over which to calculate pi/dxy.

* ``--variant_filter_expression [string]`` A comma separated list of filters (e.g. DP>=10,GQ>=20) to apply to SNPs.

* ``--invariant_filter_expression [string]`` A comma separated list of filters (e.g. DP>=10,RGQ>=20) to apply to invariant sites.

* ``--outfile_prefix [path/to/zarr/folder]`` Path and prefix for the output file. Output files will be named like: path/to/outfile_pi_[popname].txt

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
