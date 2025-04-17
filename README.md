``pixy``<img src="https://raw.githubusercontent.com/ksamuk/pixy/master/docs/images/pixy_logo.png" align="right" width="20%">
====================

[![DOI](https://zenodo.org/badge/181987337.svg)](https://zenodo.org/badge/latestdoi/181987337)
[![Anaconda Version](https://anaconda.org/conda-forge/pixy/badges/version.svg)](https://anaconda.org/conda-forge/pixy)
[![Python Versions](https://img.shields.io/badge/python-3.9_|_3.10_|_3.11-blue)](https://github.com/ksamuk/pixy)

[![CI](https://github.com/ksamuk/pixy/actions/workflows/python_package.yml/badge.svg?branch=master)](https://github.com/ksamuk/pixy/actions/workflows/python_package.yml?query=branch%3Amaster)
[![mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![Poetry](https://img.shields.io/endpoint?url=https://python-poetry.org/badge/v0.json)](https://python-poetry.org/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)

[![Install with Conda](https://img.shields.io/badge/Install%20with-conda-brightgreen.svg)](https://anaconda.org/conda-forge/pixy)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`pixy` is a command-line tool for painlessly computing unbiased estimators of population genetic summary statistics that measure genetic variation within (π, θw, Tajima’s D) and between (d<sub>xy</sub>, F<sub>ST</sub>) populations from a VCF. In particular, pixy facilitates the use of VCFs containing invariant (monomorphic) sites, which are **essential** for the correct computation of π and d<sub>xy</sub> in the face of missing data (i.e. always).

The [manuscript describing pixy](https://doi.org/10.1111/1755-0998.13326) is now published in Molecular Ecology Resources. As of version 2.0.0, a new manuscript describing the unbiased estimators of Watterson's theta and Tajima's D used by pixy is also [published in the same journal](https://doi.org/10.1111/1755-0998.14104).

## Authors
Kieran Samuk (UC Riverside) and Katharine Korunes (Duke University) <p>

## Citation
If you use `pixy` in your research, please cite the manuscripts below, and the [Zenodo DOI](https://zenodo.org/badge/latestdoi/181987337) of the specific version of pixy used for your project.

**Manuscripts**<br>
Korunes, K.L. and Samuk, K. (2021), pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13326

And, if using the unbiased estimator of Tajima's D or Watterson's theta:<br>

Bailey, N., Stevison, L., & Samuk, K. (2025). Correcting for bias in estimates of θw and Tajima’s D from missing data in next-generation sequencing. Molecular Ecology Resources, e14104. https://doi.org/10.1111/1755-0998.14104

**Zenodo DOI for various versions of pixy**<br>
Go to https://zenodo.org/badge/latestdoi/181987337 and find the DOI that matches the version used (the current version is shown first).

## Supported Organisms and Data Formats

As of version 2.0.0, pixy supports organisms with arbitrary ploidy as well as multiallelic sites. VCFs need to be compressed with bgzip and indexed with tabix. Both .tbi and .csi indexes are supported.

## Documentation

https://pixy.readthedocs.io/

## Installation

`pixy` is currently available for installation on Linux/OSX systems via conda, and [hosted on conda-forge](https://anaconda.org/conda-forge/pixy). To install pixy using conda, you will first need to add conda-forge as a channel (if you haven't already):
```
conda config --add channels conda-forge
```
Then, create and activate a new conda environment for pixy:
```
conda create -n "pixy" python=3.11
conda activate pixy
```

Then install pixy, htslib, and samtools 1.21:
```
conda install -c conda-forge pixy
conda install -c bioconda htslib
conda install -c bioconda samtools=1.21
```

You can test your pixy installation by running:
```
pixy --help
```
If you have trouble installing pixy in an environment using python 3.11, try rolling back to python 3.9.

For information in installing conda, see here:

anaconda (more features and initial modules): https://docs.anaconda.com/anaconda/install/

miniconda (lighter weight): https://docs.conda.io/en/latest/miniconda.html

## A note on accuracy
We have made every effort to ensure that pixy provides accurate and unbiased results. As described in the paper, we use population genetic simulations, where the true value of parameters is exactly known, to assess the performance of pixy. However, because of the huge biological and methodological parameter space around preparing VCFs, it is not possible to guarantee that pixy will specifically work for your organism of interest. As such, it is ultimately up to the investigator to check that pixy is performing as expected for their use case, e.g. by simulating their data-generation process, including missingness. 

## Contribute to pixy
We are very open to pull requests for new features or bugfixes. If a pull request implements a new substantial feature or fixes a substantial bug, we would be happy to considering including contributors as authors on future manuscripts decscribing new versions of pixy. See [CONTRIUBTING.md](https://github.com/ksamuk/pixy/blob/master/CONTRIBUTING.md) on how to establish a development environment for working locally.

## Development Roadmap (Planned Features as of April 2025)
- Computation of summary statistics from genotype likelihoods
- Simplified alternative to "All-Sites VCF" workflow
- Python 3.12 support
- Reduced/simplified dependencies

### Completed in pixy 2.0.0
- Update to handle GATK missing data formats (reverted due to changes in GATK)
- Support for multiallelic sites
- Support for .csi indexes
- Support for arbitrary and variable ploidy levels (including sex chromosomes) 
- Simplified contributor workflows
- Computation of Watterson's Theta and Tajima's D
