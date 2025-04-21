************
Installation
************

Via Anaconda
============
pixy is currently available for installation on Linux/OSX systems via conda-forge. You can install it using following conda install commands::

    conda create -n "pixy" python=3.11
    conda activate pixy
    conda install -c conda-forge pixy
    conda install -c bioconda htslib
    conda install -c bioconda samtools=1.21

You can test the installation by running::

    pixy --help 

.. note::
    For information on installing conda:
        anaconda (more features and initial modules): https://docs.anaconda.com/anaconda/install/
        miniconda (lighter weight): https://docs.conda.io/en/latest/miniconda.html

Via Git or Download
===================

Alternatively, you can clone from Github, with manual dependency installation: https://github.com/ksamuk/pixy
