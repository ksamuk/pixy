*********
Arguments
*********

Below is the full list of arguments accepted by ``pixy``. You can also see the
canonical help text at any time with ``pixy --help``.

Required
========

**--stats [pi|dxy|fst|watterson_theta|tajima_d]**
    Which statistics to calculate from the VCF. Provide one or more, separated
    by spaces. For example, ``--stats pi fst`` will compute π and F\ :sub:`ST`.

**--vcf [path/to/vcf.vcf.gz]**
    Path to the input VCF. The VCF must be bgzipped and indexed with ``tabix``
    (``.tbi``) or with ``bcftools index`` (``.csi``).

**--populations [path/to/populations_file.txt]**
    Path to a headerless, tab-separated populations file with columns
    ``[SampleID Population]``. See :doc:`companions` for the exact format.

In addition, one of
===================

**--window_size [integer]**
    Window size, in base pairs, over which to calculate statistics. Window
    coordinates are determined automatically across the chromosomes selected.

**--bed_file [path/to/regions.bed]**
    Path to a headerless BED file containing custom regions
    (``chrom``, ``chromStart``, ``chromEnd``) over which to compute the
    statistics. Useful when windows must match a specific genomic feature set
    (genes, introns, etc.) and may be heterogeneous in size. Coordinates follow
    the BED standard (0-based, half-open): ``chromStart`` is inclusive and
    ``chromEnd`` is exclusive. The first base of a chromosome is ``0``, so a row
    ``chr1\t0\t1000`` selects the first 1000 bases (1-based positions 1..1000).
    Pixy reports the corresponding window in ``window_pos_1`` / ``window_pos_2``
    using 1-based inclusive coordinates.

Optional
========

**--n_cores [integer]**
    Number of CPUs to use for parallel processing (default: ``1``).

**--output_folder [path/to/output/folder]**
    Folder where output will be written. Defaults to the current working
    directory.

**--output_prefix [prefix]**
    Prefix for output file(s). For example ``--output_prefix run1`` produces
    ``[output_folder]/run1_pi.txt`` and so on. Defaults to ``pixy``.

**--chromosomes ['list,of,chromosomes']**
    A single-quoted, comma-separated list of chromosomes (e.g. ``'X,1,2'``).
    Defaults to all chromosomes in the VCF.

**--interval_start [integer]**
    Start position of an interval to restrict the analysis to. Only valid
    when computing over a single chromosome. Defaults to position 1.

**--interval_end [integer]**
    End position of an interval to restrict the analysis to. Only valid
    when computing over a single chromosome. Defaults to the last position
    on the chromosome.

**--sites_file [path/to/sites_file.txt]**
    Path to a headerless, tab-separated file with columns
    ``[CHROM POS]`` defining the sites over which statistics should be
    calculated. Can be combined with ``--window_size`` and ``--bed_file``.

**--chunk_size [integer]**
    Approximate number of sites to read from the VCF at a time
    (default: ``100000``). Smaller values reduce memory use; larger values
    reduce I/O overhead.

**--fst_type [wc|hudson]**
    F\ :sub:`ST` estimator to use: ``wc`` (Weir & Cockerham 1984) or
    ``hudson`` (Hudson 1992 / Bhatia et al. 2013). Defaults to ``wc``.

**--fst_components**
    Include the F\ :sub:`ST` estimator components in the F\ :sub:`ST`
    output table. For ``--fst_type wc``, this adds ``wc_fst_a``,
    ``wc_fst_b``, and ``wc_fst_c``. For ``--fst_type hudson``, this
    adds ``hudson_fst_num`` and ``hudson_fst_den``.

**--tajima_components**
    Include the Tajima's *D* aggregation components in the Tajima's *D*
    output table. This adds ``tajima_d_s_counts``, a comma-separated list
    of ``observed_alleles:segregating_sites`` pairs that can be summed
    across windows to recompute the Tajima's *D* denominator.

**--include_multiallelic_snps**
    Include sites with more than two alleles in the calculation. Disabled by
    default because biallelic-only mode is slightly faster. Added in
    ``pixy 2.0``.

**--bypass_invariant_check**
    Skip the check that ensures invariant sites are present in the VCF.
    Disabled by default. *Use with extreme caution.* Without invariant sites,
    estimates of π and d\ :sub:`xy` are systematically biased and will be
    wrong unless your data are simulated.

**--gvcf**
    Treat the input as a joint-called GVCF in which runs of consecutive
    invariant positions are collapsed into block records
    (``ALT=<NON_REF>`` and ``INFO/END`` covering ``[POS, END]``).
    ``pixy`` expands each block back into per-site rows at read time so
    that the per-window callable-site denominator is identical to the one
    you would get from a fully-decompressed all-sites VCF. Required when
    feeding a GVCF directly — running on a GVCF without ``--gvcf`` would
    silently drop the invariant blocks. Incompatible with ``--wisp_bed``.

**--gvcf_max_block_size [integer]**
    Maximum expected GVCF block length in bp (default: ``100000``). Used
    only when ``--gvcf`` is set: ``pixy`` widens each tabix query on the
    left by this amount so blocks that start before the current window
    are not missed. The default is generous for typical GATK-emitted
    block lengths; raise it if your GVCF contains longer blocks.

**--silent**
    Suppress all console output.

**--version**
    Print the ``pixy`` version number and exit.

**--citation**
    Print the ``pixy`` citation and exit.

**--help**
    Print the full help message and exit.

Example
=======

A typical multi-statistic run including the new Watterson's θ and Tajima's
*D* estimators:

.. code:: console

    pixy --stats pi fst dxy watterson_theta tajima_d \
    --vcf data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz \
    --populations Ag1000_sampleIDs_popfile.txt \
    --window_size 10000 \
    --n_cores 4 \
    --output_folder output \
    --output_prefix pixy_output
