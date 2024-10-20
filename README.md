# Genovisio MarCNV

[![Python version](https://img.shields.io/badge/python-3.12+-green.svg)](https://www.python.org/downloads/)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-%23FE5196?logo=conventionalcommits&logoColor=white)](https://conventionalcommits.org)

Automatic evaluation of ACMG criteria for CNV based on .json file with data from MongoDB.

## Installation

In python3.12 you can simply use pip to install marCNV:

```bash
pip install git+https://github.com/geneton-ltd/genovisio_marcnv.git
```

Without python3.12, you can install marcnv using mamba:

```bash
mamba env create -f conda_example.yaml
```

This gives you the following entrypoint:

- `marcnv-classify` - run marCNV for the JSON annotation of a particular CNV region.

## Running

First, you need annotated CNV region, for example `annotation.json`. To get this annotation, see [Annotation package](https://github.com/geneton-ltd/genovisio_annotation). Then to predict, run:

```sh
marcnv-classify annotation.json --output isv.json 2> log.err
```

## Development

Poetry is used to package the application. It is required to run `poetry build` and `poetry install` to recreate the `poetry.lock` containing frozen versions of dependencies.

There are some changes whether you develop on machines with Python 3.12 and on machines without Python 3.12.
If you cannot install python3.12 and are limited to conda environments, skip to the next guide.

### Development with Python 3.12

Install poetry using pipx:

```sh
pipx install poetry
```

Now in the cloned repository, install the package:

```sh
poetry install --with dev
```

Activate the virtual environment where dependencies are installed:

```sh
poetry shell
```

All dependencies are now installed and you can run entrypoint command directly or using python `{script}`.

### Development without python 3.12

Install custom conda environment with python3.12.

Then, install poetry (python packaging management library):

```sh
curl -sSL https://install.python-poetry.org | python3 -
```

Now in the cloned repository, install the package:

```sh
poetry install --with dev
```

All dependencies are now installed and you can run entrypoint command directly or using `python {script}`.

### Adding/removing dependencies

When adding or removing dependencies, you need to define them in `pyproject.toml`.

Then, locked versions need to be redefined:

```sh
poetry lock
```

Install again:

```sh
poetry install
```

### Style and formatting

Pre-commit is used to enforce the common style and linting, defined in .pre-commit-config.yaml.

To set up pre-commit, initially run:

```sh
pre-commit install
```

Then, before each commit, an automatic linting and formatting will run and potentially prompt to review the made changes.

### Conventional PRs

When committing, you must follow the [Conventional Commits spec](https://www.conventionalcommits.org/en/v1.0.0/). Each PR is automatically validated by the GH action.
This means that there the PR title must be like `feat: XY` or `fix: XY` and so on, and the PR should contain at least one commit named like this.

Further, any push (i.e. after merged PR) to the `main` branch creates in a new PR:

- a new release following the [Semantic Versioning](https://semver.org/)
- an automatic changelog as parsed from the commit history
