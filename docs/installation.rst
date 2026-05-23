************
Installation
************

pixy is available for Linux and macOS. For a guided walkthrough — including
generating an AllSites VCF, filtering it, and running your first analysis —
see :doc:`guide/pixy_guide`. The instructions below are the minimal
copy-pasteable version for users who already know what they're doing.

Via conda
=========

We recommend installing into a fresh conda environment. If you don't already
have a conda distribution installed, we recommend `Miniforge
<https://github.com/conda-forge/miniforge>`_ (lightweight, fully free,
conda-forge pre-configured, ships with the fast ``mamba`` solver) or
`Miniconda <https://docs.anaconda.com/miniconda/>`_. The full Anaconda
distribution also works but has more restrictive licensing for commercial
or large-organization use.

``pixy`` supports Python 3.10–3.14::

    conda create -n "pixy" python=3.12
    conda activate pixy
    conda install -c conda-forge pixy
    conda install -c bioconda samtools

(The ``samtools`` install pulls in ``htslib`` as a dependency, which provides
the ``bgzip`` and ``tabix`` utilities you'll need to prepare your VCF.)

Test the installation by running::

    pixy --help

.. note::
    If the conda solver hangs on environment creation, install ``mamba`` and
    use it in place of ``conda`` for the install steps (``mamba install ...``).
    Miniforge ships with ``mamba`` already, and recent conda releases support
    ``--solver=libmamba`` for the same speedup.

From source
===========

To work from the latest code (or to contribute), clone the repository and
install in editable mode. ``pixy`` uses `Poetry
<https://python-poetry.org/>`_ for dependency management::

    git clone https://github.com/ksamuk/pixy.git
    cd pixy
    poetry install
    poetry run pixy --help

If you'd rather not install Poetry, you can install the package with ``pip``
into an existing environment::

    git clone https://github.com/ksamuk/pixy.git
    cd pixy
    pip install -e .

In either case you'll still need ``bgzip`` and ``tabix`` on your ``PATH``
(install ``samtools`` from bioconda, or your distribution's package manager,
to get them).
