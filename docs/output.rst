*************************
Understanding pixy output
*************************

Output file contents
====================

``pixy`` writes one output file per summary statistic that you request.
Each file is a tab-separated table, named
``[output_prefix]_[stat].txt`` (e.g. ``pixy_pi.txt``,
``pixy_watterson_theta.txt``). The columns in each file are documented
below.

Within-population nucleotide diversity (pi)
-------------------------------------------

File: ``[prefix]_pi.txt``

``pop``
    The ID of the population from the populations file.

``chromosome``
    The chromosome or contig.

``window_pos_1``
    First position of the genomic window.

``window_pos_2``
    Last position of the genomic window.

``avg_pi``
    Average per-site nucleotide diversity for the window. ``pixy``
    computes the *weighted average nucleotide diversity per site* over
    all sites in the window, where weights are determined by the number
    of genotyped samples at each site.

``no_sites``
    Total number of sites in the window with at least one valid
    genotype. This is included for the user and is not directly used in
    any calculation.

``count_diffs``
    Raw number of pairwise differences between all genotypes in the
    window. This is the numerator of ``avg_pi``.

``count_comparisons``
    Raw number of non-missing pairwise comparisons between all genotypes
    in the window. This is the denominator of ``avg_pi``.

``count_missing``
    Raw number of missing pairwise comparisons in the window.

Between-population nucleotide divergence (dxy)
----------------------------------------------

File: ``[prefix]_dxy.txt``

``pop1``, ``pop2``
    The IDs of the two populations being compared.

``chromosome``, ``window_pos_1``, ``window_pos_2``
    The chromosome and window coordinates (as for ``pi``).

``avg_dxy``
    Average per-site nucleotide divergence for the window.

``no_sites``
    Number of sites in the window with at least one valid genotype in
    *both* populations.

``count_diffs``, ``count_comparisons``, ``count_missing``
    Raw numerator, denominator and missing-comparison counts, defined
    analogously to the ``pi`` file but across the two populations.

F\ :sub:`ST` (fst)
------------------

File: ``[prefix]_fst.txt``

``pop1``, ``pop2``
    The IDs of the two populations being compared.

``chromosome``, ``window_pos_1``, ``window_pos_2``
    Window coordinates.

``avg_wc_fst`` *or* ``avg_hudson_fst``
    The window-averaged F\ :sub:`ST`, per SNP (not per site). Which
    column name you see depends on ``--fst_type``: ``wc`` (Weir &
    Cockerham 1984, the default) produces ``avg_wc_fst``; ``hudson``
    (Hudson 1992 / Bhatia *et al.* 2013) produces ``avg_hudson_fst``.

``no_snps``
    Total number of variable sites (SNPs) in the window.

``wc_fst_a``, ``wc_fst_b``, ``wc_fst_c``
    Present when ``--fst_components`` and ``--fst_type wc`` are used.
    These are the summed Weir & Cockerham variance components for all
    SNPs in the window.

``hudson_fst_num``, ``hudson_fst_den``
    Present when ``--fst_components`` and ``--fst_type hudson`` are used.
    These are the summed numerator and denominator terms for all SNPs in
    the window.

Watterson's θ (watterson_theta)
-------------------------------

File: ``[prefix]_watterson_theta.txt``

Watterson's θ is an estimator of the population mutation rate computed
from the number of segregating sites. The unbiased estimator implemented
in ``pixy`` corrects for missing data and arbitrary ploidy. See Bailey,
Stevison & Samuk (2025) for the derivation.

``pop``
    The ID of the population.

``chromosome``, ``window_pos_1``, ``window_pos_2``
    Window coordinates.

``avg_watterson_theta``
    Per-site Watterson's θ for the window.

``no_sites``
    Total number of sites in the window with at least one valid
    genotype in the focal population.

``raw_watterson_theta``
    Sum of per-site θ contributions over the window, used as the
    numerator of ``avg_watterson_theta``.

``no_var_sites``
    Number of segregating (variant) sites in the window that
    contributed to the estimate.

``weighted_no_sites``
    Sum of the per-site weights across the window. Used as the
    denominator of ``avg_watterson_theta``.

Tajima's *D* (tajima_d)
-----------------------

File: ``[prefix]_tajima_d.txt``

Tajima's *D* contrasts two estimators of θ — π (based on pairwise
differences) and Watterson's θ (based on segregating sites) — to detect
departures from neutrality. ``pixy`` reports the unbiased estimator
described in Bailey, Stevison & Samuk (2025), which handles missing
data correctly.

``pop``
    The ID of the population.

``chromosome``, ``window_pos_1``, ``window_pos_2``
    Window coordinates.

``tajima_d``
    Tajima's *D* for the window.

``no_sites``
    Total number of sites in the window with at least one valid
    genotype.

``raw_pi``, ``raw_watterson_theta``
    The unbiased per-site π and Watterson's θ values used to compute
    *D* (the numerator is ``raw_pi - raw_watterson_theta``).

``tajima_d_stdev``
    Standard deviation of the *D* statistic over the window (the
    denominator).

Working with pixy output data
=============================

For plotting workflows (long-format conversion, multi-statistic
panels, genome-wide plots) see :doc:`plotting`.

Post-hoc aggregating
--------------------

If you want to combine information across windows after the fact
(e.g. by averaging), **do not** simply average the per-window summary
statistics. Instead, sum the raw counts and recompute the ratio. For
``pi`` and ``dxy``:

.. parsed-literal::

    (window 1 count_diffs + window 2 count_diffs) /
    (window 1 count_comparisons + window 2 count_comparisons)

The same principle applies to Watterson's θ — sum the
``raw_watterson_theta`` and ``weighted_no_sites`` columns across
windows and divide. For Tajima's *D*, recompute from the raw π and θ
contributions; do not average ``tajima_d`` values directly across
windows.

For F\ :sub:`ST`, run with ``--fst_components`` and recompute from the
summed estimator components. For Weir & Cockerham F\ :sub:`ST`:

.. parsed-literal::

    sum(wc_fst_a) / (sum(wc_fst_a) + sum(wc_fst_b) + sum(wc_fst_c))

For Hudson F\ :sub:`ST`:

.. parsed-literal::

    sum(hudson_fst_num) / sum(hudson_fst_den)
