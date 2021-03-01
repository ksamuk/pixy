************
Arguments
************

Below is a list of required and optional arguments that pixy accepts. 

Required:

--stats        Which statistics to calculate from the VCF 
               (pi, dxy, and/or fst, separated by spaces)
--vcf           Path to the input VCF (bgzipped and tabix indexed).
--populations   Path to the populations file. See quick start for format.

In addition, one of either:

--window_size           Window size in base pairs over which to calculate pi/dxy. Defaults to whole chromosomes/contigs.
--bed_file           Path to a headerless .BED file containing regions (chrom, chromStart, chromEnd) over which to calculate summary statistics

Optional arguments:

--n_cores           Number of CPUs to utilize for parallel processing (default=1).
--output_folder           Folder where output will be written, e.g. path/to/output_folder, defaults to current working directory.
--output_prefix           Optional prefix for output file(s), e.g. \'output\' will result in writing to [output folder]/output_pi.txt, defaults to \'pixy\'.
--chromosomes            A single-quoted, comma separated list of chromosome(s) (e.g. 'X,1,2'). Defaults to all chromosomes in the VCF.
--interval_start            The start of a specific interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome. Defaults to 1.
--interval_end            The end of a specific interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome. Defaults to the last position for a chromosome.
--sites_file           Path to a tab separated file containing a headerless list of sites (CHROM, POS) to (exclusively) include in calculations 
--chunk_size           Approximate number of sites to read from VCF at any given time.  Defaults to 100000. Smaller numbers can reduce memory use.
--bypass_invariant_check            Bypass the check for invariant sites. Use with caution!
--version       Print the pixy version number.
--citation      Print the citation for pixy.
--help       Print the help message. 
--silent     Suppress all console output (flag, no value required).

An example:

.. code:: console

    pixy --stats pi fst dxy \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --window_size 10000 \
    --n_cores 4 \
    --output_folder output \
    --output_prefix pixy_output
