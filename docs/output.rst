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
    Auxiliary effective-site count that downweights sites by the fraction
    of observed alleles. This column is provided as a missingness
    diagnostic and is not used as the denominator of
    ``avg_watterson_theta``.

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
    The raw, unscaled π and Watterson's θ component sums used to compute
    *D* (the numerator is ``raw_pi - raw_watterson_theta``).

``tajima_d_stdev``
    Standard deviation of the *D* statistic over the window (the
    denominator).

``tajima_d_s_counts``
    Optional column emitted only with ``--tajima_components``. This is a
    comma-separated list of ``observed_alleles:segregating_sites`` pairs
    used to recompute ``tajima_d_stdev`` exactly when aggregating windows
    after running ``pixy``.

Working with pixy output data
=============================

For plotting workflows (long-format conversion, multi-statistic
panels, genome-wide plots) see :doc:`plotting`.

Post-hoc aggregating
--------------------

If you want to combine information across windows after the fact
(e.g. by averaging), **do not** simply average the per-window summary
statistics. Instead, sum the raw counts or components and recompute the
statistic.

Pi (pi)
^^^^^^^

For π, sum ``count_diffs`` and ``count_comparisons`` across windows, then
divide:

.. parsed-literal::

    sum(count_diffs) / sum(count_comparisons)

Dxy (dxy)
^^^^^^^^^

For d\ :sub:`xy`, use the same ratio as π:

.. parsed-literal::

    sum(count_diffs) / sum(count_comparisons)

FST (fst)
^^^^^^^^^

For F\ :sub:`ST`, run with ``--fst_components`` and recompute from the
summed estimator components. For Weir & Cockerham F\ :sub:`ST`:

.. parsed-literal::

    sum(wc_fst_a) / (sum(wc_fst_a) + sum(wc_fst_b) + sum(wc_fst_c))

For Hudson F\ :sub:`ST`:

.. parsed-literal::

    sum(hudson_fst_num) / sum(hudson_fst_den)

Watterson's θ (watterson_theta)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Watterson's θ can also be aggregated across windows. Sum
``raw_watterson_theta`` across windows and divide by the sum of
``no_sites``:

.. parsed-literal::

    sum(raw_watterson_theta) / sum(no_sites)

For Watterson's θ, ``no_var_sites`` and ``weighted_no_sites`` can also be
summed across windows as descriptive columns, but neither is the
denominator of ``avg_watterson_theta``.

Tajima's D (tajima_d)
^^^^^^^^^^^^^^^^^^^^^

Tajima's *D* cannot be exactly aggregated post hoc from the default output
columns because ``tajima_d_stdev`` is not additive. Do not average
``tajima_d`` values across windows. To enable exact post-hoc aggregation,
run with ``--tajima_components`` and sum ``tajima_d_s_counts`` by observed
allele count across the windows being combined.

Each ``tajima_d_s_counts`` entry is one or more
``observed_alleles:segregating_sites`` pairs. Here ``observed_alleles`` is
the number of called alleles at those segregating sites, not the number of
diploid individuals. For example, if two windows report
``tajima_d_s_counts`` as ``8:3,10:12`` and ``8:2,12:5``, combine them as
``8:5,10:12,12:5`` before recomputing the denominator.

This R helper mirrors the denominator calculation used by ``pixy``. Pass
``aggregate_tajima_d()`` a data frame containing the windows to combine:

.. code:: r

    parse_tajima_d_s_counts <- function(value) {
      if (length(value) == 0 || is.na(value)) {
        return(setNames(numeric(), character()))
      }

      value <- as.character(value)
      if (value == "" || value == "NA") {
        return(setNames(numeric(), character()))
      }

      counts <- setNames(numeric(), character())
      for (item in strsplit(value, ",", fixed = TRUE)[[1]]) {
        pair <- strsplit(item, ":", fixed = TRUE)[[1]]
        n <- pair[1]
        old <- counts[n]
        if (is.na(old)) {
          old <- 0
        }
        counts[n] <- old + as.numeric(pair[2])
      }
      counts
    }

    combine_tajima_d_s_counts <- function(values) {
      total <- setNames(numeric(), character())
      for (value in values) {
        counts <- parse_tajima_d_s_counts(value)
        for (n in names(counts)) {
          old <- total[n]
          if (is.na(old)) {
            old <- 0
          }
          total[n] <- old + counts[n]
        }
      }
      total
    }

    calc_tajima_d_stdev <- function(s_counts) {
      stdev <- 0
      for (n_name in names(s_counts)) {
        n <- as.integer(n_name)
        s <- as.numeric(s_counts[n_name])
        if (is.na(n) || n < 2 || s <= 0) {
          next
        }

        i <- seq_len(n - 1)
        a1 <- sum(1 / i)
        a2 <- sum(1 / (i^2))
        b1 <- (n + 1) / (3 * (n - 1))
        b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
        c1 <- b1 - (1 / a1)
        c2 <- b2 - ((n + 2) / (a1 * n)) + (a2 / (a1^2))
        e1 <- c1 / a1
        e2 <- c2 / (a1^2 + a2)

        stdev <- stdev + sqrt((e1 * s) + (e2 * s * (s - 1)))
      }
      stdev
    }

    aggregate_tajima_d <- function(rows) {
      raw_pi <- sum(rows$raw_pi, na.rm = TRUE)
      raw_watterson_theta <- sum(rows$raw_watterson_theta, na.rm = TRUE)
      s_counts <- combine_tajima_d_s_counts(rows$tajima_d_s_counts)
      tajima_d_stdev <- calc_tajima_d_stdev(s_counts)
      tajima_d <- if (tajima_d_stdev <= 0) {
        NA_real_
      } else {
        (raw_pi - raw_watterson_theta) / tajima_d_stdev
      }

      data.frame(
        tajima_d = tajima_d,
        no_sites = sum(rows$no_sites, na.rm = TRUE),
        raw_pi = raw_pi,
        raw_watterson_theta = raw_watterson_theta,
        tajima_d_stdev = tajima_d_stdev
      )
    }

The final calculation is:

.. parsed-literal::

    (sum(raw_pi) - sum(raw_watterson_theta)) / recomputed_tajima_d_stdev

If ``recomputed_tajima_d_stdev`` is zero, the aggregated Tajima's *D* is
undefined and should be reported as ``NA``. Sum ``no_sites`` across windows
for the aggregated site count. The important detail is that
``tajima_d_stdev`` itself is not summable; ``tajima_d_s_counts`` is the
summable denominator component.

When ``pixy`` internally chunks a large requested window, it uses the same
segregating-site counts by observed allele count to recompute
``tajima_d_stdev`` for the requested output window.
