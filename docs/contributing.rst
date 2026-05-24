************
Contributing
************

Contributions to ``pixy`` are very welcome — bug reports, documentation
fixes, new features, and improved tests are all useful. This page walks
through the setup needed to work on ``pixy`` locally and the commands used
to check, lint, and test the code.

If you are just looking to install ``pixy`` for use rather than development,
see :doc:`installation`. The instructions below set up a full development
environment with all dev dependencies.

.. note::
    **Not sure where to start?** Browse the
    `open issues <https://github.com/ksamuk/pixy/issues>`_ — anything tagged
    ``good first issue`` or ``documentation`` is a fine place to begin. If
    nothing obvious jumps out, fixing a typo, clarifying a docs paragraph,
    or adding a test for an under-covered code path is always welcome and
    requires very little setup.

Prerequisites
=============

Before installing the development version, install:

* `mamba <https://mamba.readthedocs.io/en/latest/index.html>`_ — a fast,
  drop-in replacement for ``conda``. (If you have ``conda`` installed but
  not ``mamba``, ``conda install -n base -c conda-forge mamba`` will add it.
  `Miniforge <https://github.com/conda-forge/miniforge>`_ ships with
  ``mamba`` already.)
* `Poetry <https://python-poetry.org/docs/#installation>`_ — used to manage
  the Python dependencies and the dev tooling.

Installing the development version
==================================

Clone the repository, create the dev environment from the provided
``environment.yaml`` (which pulls a supported Python — 3.10 through 3.14
— plus ``samtools`` and ``htslib`` from bioconda), and then install the
Python dependencies with
Poetry::

    git clone https://github.com/ksamuk/pixy.git
    cd pixy
    mamba env create -f environment.yaml
    mamba activate pixy-dev
    poetry install

You should now be able to run the development version directly::

    poetry run pixy --help

Project layout
==============

A quick map of the repository for orientation:

* ``pixy/`` — the main package.

  * ``__main__.py`` — CLI entry point and argument wiring.
  * ``args_validation.py`` — argument parsing and validation.
  * ``calc.py`` — core statistical calculations (π, d\ :sub:`xy`, F\ :sub:`ST`, θ\ :sub:`W`, Tajima's *D*).
  * ``core.py`` — windowing, parallelization, and the per-chromosome compute loop.
  * ``stats/`` — modular per-statistic implementations.
  * ``models.py`` / ``enums.py`` — shared dataclasses and enums.

* ``tests/`` — pytest suite.

  * ``tests/args_validation/`` — tests for the CLI argument layer.
  * ``tests/main/`` — end-to-end CLI tests with VCF fixtures under
    ``data/`` and golden outputs under ``expected_outputs/``.
  * ``tests/stats/`` — unit tests for the per-statistic models.
  * ``tests/test_calc.py`` — unit tests for the core calculation routines.
  * ``tests/conftest.py`` — shared pytest fixtures.

* ``docs/`` — Sphinx documentation source (this page lives here).
* ``stubs/`` — type stubs for third-party libraries without their own.
* ``environment.yaml`` — conda environment for development.
* ``pyproject.toml`` — Poetry config, dependency pins, ruff/mypy/poe tasks.

Primary development commands
============================

The following commands cover the day-to-day developer workflow. They are
listed here in long form; equivalent shortcut commands using the
``poethepoet`` task runner are described in the next section.

Resolve linting issues::

    poetry run ruff check --fix

Resolve formatting issues::

    poetry run ruff format

Run the unit tests::

    poetry run pytest

To run a single test file, a single test class, or a single test function::

    poetry run pytest tests/test_calc.py
    poetry run pytest tests/test_calc.py::TestPi
    poetry run pytest tests/test_calc.py::TestPi::test_basic_case

To run tests matching a substring (``-k``) or with verbose output (``-v``)::

    poetry run pytest -v -k "fst"

Run static type checks::

    poetry run mypy

Generate a local HTML coverage report after running tests::

    poetry run coverage html

Verify the lockfile is up to date::

    poetry check --lock

Shortcut task commands
======================

``pixy`` uses `poethepoet <https://poethepoet.natn.io/index.html>`_ to bundle
the most common dev commands into short task names. To enable them, install
the Poetry plugin once::

    poetry self add 'poethepoet[poetry_plugin]'

.. note::
    Once the project is migrated to Poetry 2.0, the plugin will install
    automatically from a declaration in ``pyproject.toml``.

Run an individual check::

    poetry task check-lock
    poetry task check-format
    poetry task check-lint
    poetry task check-tests
    poetry task check-typing

Run all checks (the same set CI runs)::

    poetry task check-all

Apply individual fixes::

    poetry task fix-format
    poetry task fix-lint

Apply all fixes::

    poetry task fix-all

Apply all fixes and then run all checks::

    poetry task fix-and-check-all

Working on the documentation
============================

The documentation lives under ``docs/`` and is built with Sphinx. To preview
your changes locally, install the doc dependencies and run a build::

    pip install -r docs/requirements.txt
    sphinx-build -b html docs/ /tmp/pixy_html

Then open ``/tmp/pixy_html/index.html`` in a browser. On Windows, substitute
a writable path such as ``%TEMP%\\pixy_html`` and open the resulting
``index.html`` with ``start``.

Code style
==========

A few light conventions that make review smoother:

* **Formatting and linting** are handled by ``ruff``, configured in
  ``pyproject.toml``. ``ruff format`` is the source of truth — please don't
  apply a different formatter on top.
* **Type hints** are checked by ``mypy``. New functions should have type
  annotations on their parameters and return values; CI will reject PRs
  whose annotations don't pass ``mypy``.
* **Docstrings** use a plain-prose style. A short one-line summary is
  enough for most internal helpers; user-facing functions and classes
  should describe parameters and behavior.
* **Imports** are sorted by ``ruff`` (isort rules) — running
  ``poetry task fix-all`` will tidy this for you.

Continuous integration
======================

Every push and pull request runs the ``Code checks`` workflow defined in
``.github/workflows/python_package.yml``. CI executes ``ruff format --check``,
``ruff check``, ``mypy``, and ``pytest`` on Python **3.9**, **3.10**, and
**3.11**. The fastest way to be confident a PR will pass CI is to run
``poetry task check-all`` locally — it runs the same set on your machine.

Reporting bugs and requesting features
======================================

Bug reports, feature requests, and questions should go to the
`GitHub issue tracker <https://github.com/ksamuk/pixy/issues>`_. We have
issue templates for bug reports, feature requests, installation help, and
general questions — please pick the one that fits, since the prompts there
collect the information we usually need.

When reporting a bug, the most useful information to include is:

* The ``pixy`` version (``pixy --version``).
* Your operating system and Python version.
* The full command that triggered the problem.
* The full error output (please paste in text, not screenshots, when
  possible).
* A small VCF that reproduces the issue, if you can share one. A few
  thousand sites from a single chromosome is plenty — the goal is something
  small enough to attach to the issue.

Pull request workflow
=====================

If you're new to contributing on GitHub, the basic flow is:

1. Fork ``ksamuk/pixy`` on GitHub.
2. Clone your fork locally and create a topic branch::

       git checkout -b fix/missing-data-edge-case

3. Make your changes, then run ``poetry task fix-and-check-all`` to format,
   lint, type-check, and test in one go.
4. Commit with a clear message describing *what* changed and *why* (see
   :ref:`commit-messages` below).
5. Push the branch to your fork and open a pull request against
   ``ksamuk/pixy:main``.

A good pull request:

* **Has a descriptive title** that summarizes the change in plain English.
* **Links any related issues** in the description (``Fixes #123``).
* **Stays focused** on a single concern — one PR per logical change makes
  review much faster.
* **Adds or updates tests** under ``tests/`` for any behavior change.
* **Updates the docs** under ``docs/`` if a user-facing flag or output
  changes.
* **Passes** ``poetry task check-all`` locally before pushing.

Don't worry about getting every detail right on the first push — reviewers
are happy to suggest small adjustments. The most important thing is to open
the PR and start the conversation.

.. _commit-messages:

Commit messages
===============

We don't enforce a strict commit-message format, but PRs are much easier to
review (and the git history is much more useful) when commits follow a few
loose conventions:

* **Use the imperative mood** in the subject line: "Add Tajima's D
  example", not "Added Tajima's D example" or "Adds Tajima's D example".
* **Keep the subject line short** — aim for roughly 50 characters and no
  hard period at the end.
* **Separate subject and body with a blank line.** If your change needs
  explanation, write a body that says *what* changed and, more importantly,
  *why*. Wrap body lines at ~72 characters.
* **Reference issues** in the body (``Refs #123``) or, when the commit
  closes one, in the PR description (``Fixes #123``).
* **One logical change per commit** when feasible. If you find yourself
  writing "and" in a commit subject, it's probably two commits.

A short prefix tag like ``docs:``, ``fix:``, or ``test:`` at the start of
the subject is welcome but optional.

Code of conduct
===============

``pixy`` follows the `Contributor Covenant
<https://www.contributor-covenant.org/version/2/1/code_of_conduct/>`_
(version 2.1). Participation in the project — opening issues, submitting
pull requests, commenting on either — implies that you'll treat other
contributors with respect. The full text is reproduced in
``CODE_OF_CONDUCT.md`` at the repo root, and instances of unacceptable
behavior may be reported to the project maintainer (see :doc:`about`).

License
=======

``pixy`` is released under the MIT License (see ``LICENSE`` in the repo
root). By contributing code, documentation, or other materials to the
project, you agree that your contributions will be made available under the
same license.

Release process (for maintainers)
=================================

For maintainers cutting a new release of ``pixy``:

1. **Update the version.** Bump the version in ``pyproject.toml`` (and in
   the ``release`` string in ``docs/conf.py`` if it has drifted).
2. **Update the changelog.** Add a new section to ``docs/changelog.rst``
   summarizing the user-visible changes, grouped into Features, Bug fixes,
   and Packaging as appropriate.
3. **Verify everything passes.** Run ``poetry task check-all`` and rebuild
   the docs (``sphinx-build -b html docs/ /tmp/pixy_html``) to confirm
   nothing regressed.
4. **Tag the release.** ``git tag -a v<version> -m "pixy <version>"`` and
   ``git push --tags``.
5. **Release on GitHub.** Draft a new Release on the GitHub repo using the
   tag; paste the changelog section as the release notes.
6. **Update the conda-forge feedstock.** Open a PR against
   `conda-forge/pixy-feedstock <https://github.com/conda-forge/pixy-feedstock>`_
   bumping the version and sha256 in ``recipe/meta.yaml``. The
   conda-forge bots usually open this PR automatically a few hours after
   the GitHub release; you only need to do it manually if the bot misses
   it.
7. **Archive on Zenodo.** Once GitHub-Zenodo integration is enabled for the
   repo, every GitHub Release is automatically archived with a new DOI.
   Confirm the new DOI appears at https://zenodo.org/record/4432294 and
   that the version metadata is correct.
