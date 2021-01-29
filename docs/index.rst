.. pixy documentation master file, created by
   sphinx-quickstart on Tue Aug 13 10:19:37 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pixy's documentation!
================================
.. image:: pixy_logo.png
   :width: 200
   :align: center

.. centered::
   *pixy: Unbiased estimates of pi, dxy (and Fst) from VCFs containing invariant sites.*

What is pixy?
===============
pixy is a command line tool for calculating the population genetic summary statistics **pi** (average per site heterozygosity) and **dxy** (average number of nucleotide differences between populations per site) from a VCF file. Many major software methods for computing pi and dxy from VCFs produce biased estimates in the presence of missing data. This is because these methods make the (often implicit) simplifying assumption that if a site (or genotype) is missing, it counts as a "0" (an invariant site), BUT also contributes to the denominator (the total number of sites and/or genotypes). This can result in a substantial deflation in estimates of pi and dxy, a bias that scales with the amount of missing data. pixy is specifically designed to work with VCFs containing invariant sites and provide unbiased estimates in the presence of missing data. See pixy's paper (https://doi.org/10.1111/1755-0998.13326) for more details.

How should I cite pixy?
=======================
If you use pixy in your research, please cite the manuscript below, as well the Zenodo DOI of specific version of pixy used for your project..

**Manuscript:**
Korunes, K.L. and Samuk, K. (2021), pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data. Molecular Ecology Resources. Accepted Author Manuscript. https://doi.org/10.1111/1755-0998.13326

**Zenodo DOI for various versions of pixy:**
Go to https://zenodo.org/record/4432294 and find the DOI that matches the version used (the current version is shown first). 


.. toctree::
   :caption: Utility documentation
   :maxdepth: 2

   about
   installation
   quickstart
   arguments
   examples

Tutorials
==================

.. toctree::
    :maxdepth: 2
    :caption: Generating pixy Input

    invar/allsitesQuickstart
    
.. toctree::
    :maxdepth: 2
    :caption: Step by Step Installation & Usage

    quickstart/quickstart
