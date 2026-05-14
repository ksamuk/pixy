.. pixy documentation master file, created by
   sphinx-quickstart on Tue Aug 13 10:19:37 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html

    <div align="center"><h1> pixy 2.0.0.beta14 </h1></div>

.. image:: images/pixy_logo.png
   :width: 200
   :align: center

What is pixy?
=============
``pixy`` is a command-line tool for painlessly computing **unbiased** estimators of population genetic summary statistics that measure genetic variation within (π, θ\ :sub:`W`, Tajima's *D*) and between (d\ :sub:`xy`, F\ :sub:`ST`) populations from a VCF.

Many tools for computing these summary statistics from VCFs produce biased estimates in the presence of missing data. This is because they often make the simplifying assumption that if a genotype is missing, it is homozygous reference (``0/0``) by state. See the `pixy paper <https://doi.org/10.1111/1755-0998.13326>`_ and the `Watterson's θ / Tajima's D follow-up paper <https://doi.org/10.1111/1755-0998.14104>`_ for the details.

As of version 2.0, ``pixy`` also supports organisms of arbitrary (and variable) ploidy, multiallelic sites, and both ``.tbi`` and ``.csi`` VCF indexes.

.. toctree::
   :caption: Documentation
   :maxdepth: -1

   about
   installation
   arguments
   companions
   contributing
   changelog


.. toctree::
    :maxdepth: -1
    :caption: Guides

    generating_invar/generating_invar
    guide/pixy_guide
    examples
    example_data
    output
    plotting

How should I cite pixy?
=======================
If you use pixy in your research, please cite the manuscript below, as well the Zenodo DOI of specific version of pixy used for your project.

**Manuscript:**
Korunes, K.L. and Samuk, K. (2021), pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data. Molecular Ecology Resources. Accepted Author Manuscript. https://doi.org/10.1111/1755-0998.13326

And, if using the unbiased estimator of Tajima's D or Watterson's theta:

Bailey, N., Stevison, L., & Samuk, K. (2025). Correcting for bias in estimates of θw and Tajima’s D from missing data in next-generation sequencing. Molecular Ecology Resources, e14104. https://doi.org/10.1111/1755-0998.14104

**Zenodo DOI for various versions of pixy:**
Go to https://zenodo.org/record/4432294 and find the DOI that matches the version used (the current version is shown first). 



