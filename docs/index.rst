.. pixy documentation master file, created by
   sphinx-quickstart on Tue Aug 13 10:19:37 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html

    <div align="center"><h1> pixy 2.0.0.beta4 </h1></div>

.. image:: images/pixy_logo.png
   :width: 200
   :align: center
   
What is pixy?
------------
pixy is a command-line tool for painlessly computing unbiased estimators of population genetic summary statistics that measure genetic variation within (π, θW, Tajima’s D) and between (dxy, FST) populations from a VCF.

Many tools for computing summary statistics using VCFs produce biased estimates in the presence of missing data. This is because these methods often make the simplifying assumption that if a genotype is missing, it is homozygous reference (0/0) by state. See pixy's paper (https://doi.org/10.1111/1755-0998.13326) for more details.

.. toctree::
   :caption: Documentation
   :maxdepth: -1

   about
   installation
   arguments
   companions
   changelog
   

.. toctree::
    :maxdepth: -1
    :caption: Guides

    generating_invar/generating_invar
    guide/pixy_guide
    examples
    output
    plotting

How should I cite pixy?
------------
If you use pixy in your research, please cite the manuscript below, as well the Zenodo DOI of specific version of pixy used for your project.

**Manuscript:**
Korunes, K.L. and Samuk, K. (2021), pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data. Molecular Ecology Resources. Accepted Author Manuscript. https://doi.org/10.1111/1755-0998.13326

And, if using the unbiased estimator of Tajima's D or Watterson's theta:

Bailey, N., Stevison, L., & Samuk, K. (2025). Correcting for bias in estimates of θw and Tajima’s D from missing data in next-generation sequencing. Molecular Ecology Resources, e14104. https://doi.org/10.1111/1755-0998.14104

**Zenodo DOI for various versions of pixy:**
Go to https://zenodo.org/record/4432294 and find the DOI that matches the version used (the current version is shown first). 



