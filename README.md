``pixy``<img src="https://raw.githubusercontent.com/ksamuk/pixy/master/docs/images/pixy_logo.png" align="right" width="20%">
====================

[![DOI](https://zenodo.org/badge/181987337.svg)](https://zenodo.org/badge/latestdoi/181987337) ![version](https://img.shields.io/badge/version-1.0.1.beta1-blue) [![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs) 

`pixy` is a command-line tool for painless and unbiased estimates of nucleotide diversity within (π) and between (d<sub>xy</sub>) populations from a VCF. In particular, pixy facilitates the use of VCFs containing invariant (monomorphic) sites, which are essential for the correct computation of π and d<sub>xy</sub> in the face of missing data (i.e. always).

The [manuscript describing pixy](https://doi.org/10.1111/1755-0998.13326) is now published in Molecular Ecology Resources.

## Authors
Kieran Samuk and Katharine Korunes <p>
Duke University
 
## Citation
If you use `pixy` in your research, please cite the manuscript below, as well the [Zenodo DOI](https://doi.org/10.5281/zenodo.4432294) of specific version of pixy used for your project..

**Manuscript**<br>
Korunes, K.L. and Samuk, K. (2021), pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data. Molecular Ecology Resources. Accepted Author Manuscript. https://doi.org/10.1111/1755-0998.13326

**Zenodo DOI for various versions of pixy**<br>
Go to https://zenodo.org/record/4596999 and find the DOI that matches the version used (the current version is shown first).

## Documentation

https://pixy.readthedocs.io/

## Installation

`pixy` is currently available for installation on Linux/OSX systems via conda-forge. To install pixy using conda, you will first need to add conda-forge as a channel (if you haven't already):
```
conda config --add channels conda-forge
```

Then install pixy and htslib:
```
conda install -c conda-forge pixy
conda install -c bioconda htslib
```

You can test your pixy installation by running:
```
pixy --help
```

For information in installing conda, see here:

anaconda (more features and initial modules): https://docs.anaconda.com/anaconda/install/

miniconda (lighter weight): https://docs.conda.io/en/latest/miniconda.html
