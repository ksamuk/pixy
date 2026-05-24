************
Example data
************

A small all-sites VCF and matching populations file are included with the
``pixy`` source tree so users can verify their installation and reproduce
the figures shown in the documentation.

Download
========

The files live under ``docs/example_data/`` in the ``pixy`` GitHub
repository. They are deliberately tiny so they can be downloaded
quickly and committed to the repo without bloating clones.

.. list-table::
   :header-rows: 1
   :widths: 32 12 56

   * - File
     - Size
     - Description
   * - `dmel_2pop.vcf.gz <https://github.com/ksamuk/pixy/raw/master/docs/example_data/dmel_2pop.vcf.gz>`_
     - ~200 KB
     - bgzipped all-sites VCF, 50,000 sites on a single chromosome
       (``chr2L``), 20 diploid samples split across two populations.
   * - `dmel_2pop.vcf.gz.tbi <https://github.com/ksamuk/pixy/raw/master/docs/example_data/dmel_2pop.vcf.gz.tbi>`_
     - <1 KB
     - tabix index for the VCF.
   * - `dmel_populations.txt <https://github.com/ksamuk/pixy/raw/master/docs/example_data/dmel_populations.txt>`_
     - <1 KB
     - Populations file: 10 samples in ``pop_A``, 10 in ``pop_B``.

Or from the command line:

.. code:: console

    BASE=https://github.com/ksamuk/pixy/raw/master/docs/example_data
    curl -LO $BASE/dmel_2pop.vcf.gz
    curl -LO $BASE/dmel_2pop.vcf.gz.tbi
    curl -LO $BASE/dmel_populations.txt

How the data were generated
===========================

The VCF was simulated with `vcfsim
<https://github.com/samuk-lab/vcfsim>`_, a thin wrapper around
`msprime <https://tskit.dev/msprime/>`_ for producing all-sites VCFs.
The demographic model is a single split: an ancestral population *C*
splits into two derived populations *A* and *B* at ``--div_time``
generations before present.

Parameters were taken from the
`stdpopsim <https://popsim-consortium.github.io/stdpopsim-docs/>`_
*Drosophila melanogaster* species model:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Parameter
     - Value
     - Source
   * - Mutation rate (μ)
     - 5.49 × 10\ :sup:`-9` per site per generation
     - Keightley *et al.* 2014 (stdpopsim DroMel)
   * - Effective population size (N\ :sub:`e`)
     - 1.65 × 10\ :sup:`6`
     - Li & Stephan 2006 (stdpopsim DroMel)
   * - Divergence time
     - 500,000 generations
     - Chosen to give a moderate F\ :sub:`ST` (≈ 0.1–0.2)
   * - Sequence length
     - 50,000 bp
     - Single chromosome ``chr2L``
   * - Samples per population
     - 10 diploid individuals
     -
   * - Ploidy
     - 2
     -
   * - Missing data
     - 0%
     - Clean baseline; see notes below.

The exact command used was:

.. code:: console

    vcfsim \
      --seed 2024 \
      --replicates 1 \
      --chromosome chr2L \
      --sequence_length 50000 \
      --ploidy 2 \
      --Ne 1650000 \
      --mu 5.49e-9 \
      --population_mode 2 \
      --div_time 500000 \
      --sample_size 20 \
      --percent_missing_sites 0 \
      --percent_missing_genotypes 0 \
      --output_file dmel_2pop

The output ``dmel_2pop2024.vcf`` was then renamed, bgzipped, and indexed:

.. code:: console

    mv dmel_2pop2024.vcf dmel_2pop.vcf
    bgzip dmel_2pop.vcf
    tabix -p vcf dmel_2pop.vcf.gz

The populations file was constructed by hand based on ``vcfsim``'s
sample ordering. With ``--population_mode 2 --sample_size 20``, the
samples emitted are ``tsk_1`` through ``tsk_20``: the first 10
(``tsk_1``–``tsk_10``) are drawn from population *A* and the next 10
(``tsk_11``–``tsk_20``) from population *B*. See the ``vcfsim``
source for the exact assignment.

Running pixy on the example
===========================

.. code:: console

    pixy --stats pi fst dxy watterson_theta tajima_d \
      --vcf dmel_2pop.vcf.gz \
      --populations dmel_populations.txt \
      --window_size 5000 \
      --n_cores 2 \
      --output_prefix dmel_example

This produces ``dmel_example_pi.txt``, ``dmel_example_fst.txt``,
``dmel_example_dxy.txt``, ``dmel_example_watterson_theta.txt`` and
``dmel_example_tajima_d.txt`` in your working directory. See
:doc:`output` for the meaning of each column.

Expected results
----------------

With the parameters above (μ = 5.49×10\ :sup:`-9`,
N\ :sub:`e` = 1.65×10\ :sup:`6`), the infinite-sites expectation
for per-site heterozygosity in either daughter population is
θ = 4·N\ :sub:`e`·μ ≈ 0.0362. ``vcfsim`` (via ``msprime``) uses the
JC69 finite-sites mutation model by default, so the same site can
be hit by repeated or back-mutations and the expected *observable*
proportion of pairwise differences is slightly lower:

.. parsed-literal::

    π ≈ θ / (1 + 4θ/3) ≈ 0.035

and d\ :sub:`xy` between the two populations is slightly higher
because of the divergence. F\ :sub:`ST` should fall in the
0.1–0.2 range. With only 50 kb of sequence and 10 individuals per
population there is substantial sampling noise, so per-window
estimates will vary widely; aggregating across windows (see
:doc:`output`) recovers values close to the expectations.

.. note::
    These are theoretical values from neutral coalescent theory.
    Empirical π in *D. melanogaster* is lower (≈ 0.005–0.01) because
    selection and linkage reduce diversity. The simulation used here
    is neutral, so it reflects the full neutral expectation rather
    than what is observed in real flies.

Variants of the same data
-------------------------

For testing how ``pixy`` handles imperfect data, you can regenerate
the VCF with missing data by changing the missingness arguments to
``vcfsim``. A common pair is uniform 5% site-level and 2% genotype-level
dropout:

.. code:: console

    vcfsim ... \
      --percent_missing_sites 0.05 \
      --percent_missing_genotypes 0.02 \
      --output_file dmel_2pop_missing

This is exactly the scenario ``pixy`` was designed for: comparing the
``avg_pi`` produced by ``pixy`` against the theoretical value above
demonstrates that the unbiased estimator handles missing data
correctly, while the equivalent ``VCFtools`` ``--site-pi`` call gives a
deflated estimate.
