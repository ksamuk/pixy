************
Arguments
************

Below is a list of required and optional arguments that pixy accepts. 

--help       Print the help message. 
--version       Print the pixy version number.
--stats         **Required.** Which statistics to calculate from the VCF 
               (pi, dxy, and/or fst, separated by spaces)
--vcf           **Required.** Path to the input VCF.
--populations            **Required.** Path to the populations file. See quick start for format.
--chromosomes            A single-quoted, comma separated list of chromosome(s) (e.g. 'X,1,2'). Defaults to all chromosomes in the VCF.
--window_size           Window size in base pairs over which to calculate pi/dxy. Defaults to whole chromosomes/contigs.
--interval_start            The start of a specific interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome.
--interval_end            The end of a specific interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome.
--bed_file           Path to a headerless .BED file containing regions (chrom, chromStart, chromEnd) over which to calculate summary statistics
--sites_file           Path to a tab separated file containing a headerless list of sites (CHROM, POS) at which all statistics will be calculated (exclusively)
--bypass_invariant_check            Bypass the check for invariant sites. Use with caution!
--output_folder           Folder where output will be written, e.g. path/to/output_folder, defaults to current working directory.
--output_prefix           Optional prefix for output file(s), e.g. \'output\' will result in writing to [output folder]/output_pi.txt, defaults to \'pixy\'.
--n_cores           Number of cores to use for parallel processing. Defaults to 2.
--chunk_size           Approximate number of sites to read from VCF at any given time.  Defaults to 100000. Smaller numbers can reduce memory use.

An example:

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --window_size 10000 \
    --populations Ag1000_sampleIDs_popfile.txt \
    --output_folder output
