*************************
Interpreting Tajima's *D*
*************************

Tajima's *D* is easy to compute and easy to over-read. This page covers what
``pixy``'s *D* does and does not tell you, what ``pixy`` counts when computing
it, and how to use it on real data.

The short version
=================

- *D* is sound for **comparing** windows, populations, or groups computed the
  same way on the same data. That is what it is for.
- *D* is **not a p-value**. Its absolute value is not a test of neutrality:
  under neutrality its expectation is not 0 and its standard deviation is not 1.
- ``pixy`` counts **mutations**, not segregating sites. On biallelic data this is
  identical to the classic estimator; at multiallelic sites it departs by design
  from tools that count sites.
- Multiallelic sites in real data are enriched for genotyping artifacts, and
  θ\ :sub:`W` and *D* weight them more heavily than π does. Filter accordingly.

Why *D* is not a p-value
========================

Tajima's *D* is often taught as "0 under neutrality, negative after a sweep,
positive under balancing selection," with ``|D| > 2`` as a rule of thumb for
significance. That framing relies on assumptions that essentially never hold for
whole-genome resequencing data, and it is not specific to ``pixy`` — it applies
to every implementation of the statistic.

**The expectation is not 0.** Tajima's *D* is a ratio whose denominator is itself
a random variable, so E[*D*] ≠ 0 even under a strictly neutral, infinite-sites
coalescent. Tajima (1989) used a beta approximation for significance precisely
because of this. Separately, real sequence data violates the infinite-sites
assumption: when a site is hit by mutation more than once, both π and θ\ :sub:`W`
undercount, and *D* settles at a small negative **finite-sites floor** rather
than at 0. That floor deepens with missingness and depends on θ.

**The standard deviation is not 1.** The variance term (Tajima 1989) is derived
for a *single non-recombining locus*. Real windows contain many recombining
genealogies, which average out and shrink the spread of *D* substantially.

.. note::

   Versions 2.1.0 through 2.2.1 additionally inflated the denominator whenever
   missing data made the per-site observed-allele count ragged, shrinking
   \|*D*\| toward 0 — badly enough to matter (a true *D* of −1.84 was reported
   as −1.38 at 20% missing genotypes). This was a regression of the fix for
   `issue #160 <https://github.com/ksamuk/pixy/issues/160>`_ and is corrected;
   see :doc:`changelog`. Data with no missing genotypes was never affected. If
   you have Tajima's D from an affected version computed on data with missing
   genotypes, recompute it.

In ``pixy``'s polyploid validation (θ ≈ 0.038, 10 kb windows, JC69 mutation) the
measured values were:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Regime
     - mean *D*
     - sd *D*
   * - No recombination, no recurrent mutation
     - ≈ −0.11
     - ≈ 0.87
   * - No recombination, recurrent mutation
     - ≈ −0.08 to −0.14
     - ≈ 0.83
   * - Recombination at r = μ (realistic)
     - finite-sites floor
     - ≈ 0.19

These are **not universal constants** — they depend on θ, window size,
recombination rate, sample size, and missingness. They are reported here to make
one point: the quantity you would need to divide by to turn *D* into a z-score is
not 1, is not knowable from the formula, and varies several-fold across regimes.

**What this means in practice.** Use *D* as a **relative** measure. Comparing *D*
across windows within a dataset, or between populations processed identically, is
exactly what the statistic supports and what ``pixy`` is built for. Reading a
single window's *D* = −0.4 as evidence of a sweep is not supported. If you need
calibrated significance, you need a null simulated to match your data's θ,
ploidy, missingness, and per-window recombination rate — an analysis ``pixy``
does not perform for you.

What ``pixy`` counts
====================

Watterson's estimator is derived from the number of **mutations on the
genealogy**: E[η] = a\ :sub:`n`\ θ follows from Poisson mutation along the
coalescent tree. The familiar "number of segregating sites" is a shorthand for
that quantity, exact only under strict infinite sites, where every mutation
creates a new site and the two counts coincide.

A site carrying three or four alleles is, by definition, a place where that
shorthand breaks. ``pixy`` therefore counts a site with *k* observed alleles as
*k* − 1 mutations — Tajima's ``s*``, the minimum (parsimony) number of mutations
per site (Tajima 1996) — in both θ\ :sub:`W` and the *D* variance term.

Consequences worth knowing:

- **On biallelic data nothing changes.** Every biallelic site has *k* − 1 = 1, so
  ``pixy`` reduces exactly to the classic estimator and agrees with
  ``scikit-allel``.
- **On multiallelic data ``pixy`` deliberately differs** from ``scikit-allel``,
  ``vcftools``, and other site-counting implementations. This is a design choice,
  not a discrepancy to reconcile. Report it in your methods.
- **A site fixed for the alternate allele is not segregating.** Every sample being
  homozygous ``1/1`` means one observed allele and zero mutations, even though the
  site is a variant relative to the reference. ``pixy`` counts it as 0. (Versions
  through 2.2.1 incorrectly counted such sites as segregating; see
  :doc:`changelog`.)
- **θ\ :sub:`W` carries a small downward bias under recurrent mutation.** ``s*``
  is a parsimony *minimum*: homoplasy and back-mutation are invisible to it. In
  the validation regime above this left θ\ :sub:`W` ≈ 4% below 4N\ :sub:`e`\ μ,
  flat across ploidy. No model-free estimator can remove this; doing so requires
  committing to a specific mutation model.

Multiallelic sites on real data
===============================

``pixy``'s multiallelic validation simulates a finite-sites Jukes–Cantor process,
in which **every** third or fourth allele is a genuine recurrent mutation. Real
multiallelic calls are not so clean. Collapsed paralogs, mismapping, and
indel-adjacent miscalls all manufacture spurious extra alleles, and they do so
preferentially in the repetitive regions where mapping is hardest.

This matters more for θ\ :sub:`W` and *D* than for π. Because a *k*-allele site
contributes *k* − 1 mutations, a spurious **third** allele contributes twice what
a spurious second allele would, whereas any single site's contribution to π is
bounded. The validation establishes that the estimators are correct *given* the
genotypes; it does not establish that a given set of multiallelic calls has
earned that trust.

How much will this affect me?
-----------------------------

The size of the multiallelic contribution scales with θ and with ploidy. From
``pixy``'s θ sweep, here is how much diversity biallelic-only π *misses* — a
useful proxy for how much multiallelic sites matter in your data:

.. list-table::
   :header-rows: 1
   :widths: 34 33 33

   * - θ
     - diploid
     - octoploid
   * - 0.005
     - −0.8%
     - −0.6%
   * - 0.01
     - −1.5%
     - −2.0%
   * - 0.025
     - −3.8%
     - −5.4%
   * - 0.05
     - −7.1%
     - −10.3%
   * - 0.10
     - −12.5%
     - −17.6%

At human-like diversity (θ ≈ 0.001) multiallelic sites are a rounding error;
enable the flag and move on. At high diversity, in polyploids, or with large
samples, they matter enough to be worth the care described below.

A practical checklist
---------------------

1. **Run both modes and diff them.** Run ``pixy`` with and without
   ``--include_multiallelic_snps``. The difference *is* the multiallelic
   contribution. If it is far larger than the table above predicts for your θ,
   your multiallelic calls are telling you about your pipeline, not your
   population.
2. **Check the allele spectrum.** Genuine recurrent mutation produces mostly
   3-allele sites and rarely 4-allele ones. A pile-up of 4-allele sites, or
   multiallelic sites clustered in repeats or beside indels, is a mapping
   signature rather than a biological one.
3. **Filter multiallelic sites at least as hard as biallelic ones.** Depth caps,
   mapping quality, and repeat masking all apply. The *k* − 1 weighting means a
   spurious triallelic site costs you double.
4. **Compare the two modes across windows**, not their absolute values.

Reporting in a methods section
==============================

Because ``pixy``'s multiallelic θ\ :sub:`W` and *D* differ by design from
site-counting implementations, state what you used:

    Watterson's θ and Tajima's *D* were computed with ``pixy`` vX, which counts
    mutations (Tajima's ``s*``, the minimum mutation count Σ(*k* − 1)) rather than
    segregating sites. On biallelic data this is identical to the standard
    estimators; at multiallelic sites it departs by design from implementations
    that count sites (e.g. ``scikit-allel``, ``vcftools``). θ estimates carry an
    irreducible ~4% downward bias under recurrent mutation.

References
==========

- Tajima, F. (1989). Statistical method for testing the neutral mutation
  hypothesis by DNA polymorphism. *Genetics* 123: 585–595.
- Tajima, F. (1996). The amount of DNA polymorphism maintained in a finite
  population when the neutral mutation rate varies among sites. *Genetics*
  143: 1457–1465.
- Watterson, G.A. (1975). On the number of segregating sites in genetical models
  without recombination. *Theoretical Population Biology* 7: 256–276.
- Bailey, N., Stevison, L. & Samuk, K. (2025). Correcting for bias in estimates
  of θ\ :sub:`W` and Tajima's *D* from missing data in next-generation
  sequencing. *Molecular Ecology Resources* e14104.
- Roychoudhury, A. & Wakeley, J. (2010). Sufficiency of the number of segregating
  sites in the limit under finite-sites mutation. *Theoretical Population
  Biology* 78: 118–122.
