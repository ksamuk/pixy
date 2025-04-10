``pixy``<img src="https://raw.githubusercontent.com/ksamuk/pixy/master/docs/images/pixy_logo.png" align="right" width="20%">
====================

[![DOI](https://zenodo.org/badge/181987337.svg)](https://zenodo.org/badge/latestdoi/181987337)
[![Anaconda Version](https://anaconda.org/conda-forge/pixy/badges/version.svg)](https://anaconda.org/conda-forge/pixy)
[![Python Versions](https://img.shields.io/badge/python-3.8_|_3.9_|_3.10_|_3.11-blue)](https://github.com/ksamuk/pixy)

[![CI](https://github.com/ksamuk/pixy/actions/workflows/python_package.yml/badge.svg?branch=master)](https://github.com/ksamuk/pixy/actions/workflows/python_package.yml?query=branch%3Amaster)
[![mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![Poetry](https://img.shields.io/endpoint?url=https://python-poetry.org/badge/v0.json)](https://python-poetry.org/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)

[![Install with Conda](https://img.shields.io/badge/Install%20with-conda-brightgreen.svg)](https://anaconda.org/conda-forge/pixy)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`pixy` is a command-line tool for painlessly estimating average nucleotide diversity within (π) and between (d<sub>xy</sub>) populations from a VCF. In particular, pixy facilitates the use of VCFs containing invariant (monomorphic) sites, which are **essential** for the correct computation of π and d<sub>xy</sub> in the face of missing data (i.e. always).

The [manuscript describing pixy](https://doi.org/10.1111/1755-0998.13326) is now published in Molecular Ecology Resources.

## Authors
Kieran Samuk (UC Riverside) and Katharine Korunes (Duke University) <p>

## Citation
If you use `pixy` in your research, please cite the manuscript below, and the [Zenodo DOI](https://zenodo.org/badge/latestdoi/181987337) of the specific version of pixy used for your project.

**Manuscript**<br>
Korunes, K.L. and Samuk, K. (2021), pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data. Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13326

**Zenodo DOI for various versions of pixy**<br>
Go to https://zenodo.org/badge/latestdoi/181987337 and find the DOI that matches the version used (the current version is shown first).

## Supported Organisms and Data Formats

Currently, pixy only supports computation using biallelic SNPs (and invariant sites) from diploid organisms. VCFs need to be compressed with bgzip and indexed with tabix.

## Documentation

https://pixy.readthedocs.io/

## Installation

`pixy` is currently available for installation on Linux/OSX systems via conda, and [hosted on conda-forge](https://anaconda.org/conda-forge/pixy). To install pixy using conda, you will first need to add conda-forge as a channel (if you haven't already):
```
conda config --add channels conda-forge
```
Then, create and activate a new conda environment for pixy:
```
conda create -n "pixy" python=3.8
conda activate pixy
```

Then install pixy, htslib, and samtools 1.9:
```
conda install -c conda-forge pixy=1.2.5
conda install -c bioconda htslib
conda install -c bioconda samtools=1.9 --force-reinstall -y
```

You can test your pixy installation by running:
```
pixy --help
```
If you have trouble installing pixy in an environment using python 3.9, try rolling back to python 3.8.

For information in installing conda, see here:

anaconda (more features and initial modules): https://docs.anaconda.com/anaconda/install/

miniconda (lighter weight): https://docs.conda.io/en/latest/miniconda.html

## A note on accuracy
We have made every effort to ensure that pixy provides accurate and unbiased results. As described in the paper, we use population genetic simulations, where the true value of parameters is exactly known, to assess the performance of pixy. However, because of the huge biological and methodological parameter space around preparing VCFs, it is not possible to guarantee that pixy will specifically work for your organism of interest. As such, it is ultimately up to the investigator to check that pixy is performing as expected for their use case, e.g. by simulating their data-generation process, including missingness. 

## Contribute to pixy
We are very open to pull requests for new features or bugfixes. If a pull request implements a new substantial feature or fixes a substantial bug, we would be happy to considering including contributors as authors on future manuscripts decscribing new versions of pixy.

## Development Roadmap (Planned Features as of Feb 2024)
- Update to handle GATK missing data formats - COMPLETE (as of version 1.2.11.beta1)
- Simplified alternative to "All-Sites VCF" workflow
- Support for arbitrary and variable ploidy levels (including sex chromosomes)
- Computation of summary statistics from genotype likelihoods
- Simplified contributor workflows
- Computation of Tajima's D
