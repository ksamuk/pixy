# Development and Testing

## Installing development version of pixy

Create a fresh conda environment with Python 3.8, install samtools and htslib
1.9 using the provided `environment.yaml`, and then install `pixy` and its
Python dependencies with poetry.

```console
mamba create -n pixy-dev python=3.8
mamba env update -f environment.yaml
poetry install
```

## Primary Development Commands

To check and resolve linting issues in the codebase, run:

```console
poetry run ruff check --fix
```

To check and resolve formatting issues in the codebase, run:

```console
poetry run ruff format
```

To check the unit tests in the codebase, run:

```console
poetry run pytest
```

To check the typing in the codebase, run:

```console
poetry run mypy
```

To generate a code coverage report after testing locally, run:

```console
poetry run coverage html
```

To check the lock file is up to date:

```console
poetry check --lock
```

## Shortcut Task Commands

To be able to run shortcut task commands, first install the Poetry plugin [`poethepoet`](https://poethepoet.natn.io/index.html):

```console
poetry self add 'poethepoet[poetry_plugin]'
```

> [!NOTE]
> Upon the release of Poetry [v2.0.0](https://github.com/orgs/python-poetry/discussions/9793#discussioncomment-11043205), Poetry will automatically support bootstrap installation of [project-specific plugins](https://github.com/python-poetry/poetry/pull/9547) and installation of the task runner will become automatic for this project.
> The `pyproject.toml` syntax will be:
> 
> ```toml
> [tool.poetry]
> requires-poetry = ">=2.0"
> 
> [tool.poetry.requires-plugins]
> poethepoet = ">=0.29"
> ```

###### For Running Individual Checks

```console
poetry task check-lock
poetry task check-format
poetry task check-lint
poetry task check-tests
poetry task check-typing
```

###### For Running All Checks

```console
poetry task check-all
```

###### For Running Individual Fixes

```console
poetry task fix-format
poetry task fix-lint
```

###### For Running All Fixes

```console
poetry task fix-all
```

###### For Running All Fixes and Checks

```console
poetry task fix-and-check-all
```
