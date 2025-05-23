# MOCPy

[![PyPI version](https://badge.fury.io/py/mocpy.svg)](https://badge.fury.io/py/MOCPy)
[![Build/Test status](https://github.com/cds-astro/mocpy/actions/workflows/test.yml/badge.svg)](https://github.com/cds-astro/mocpy/actions/workflows/test.yml)
[![Doc](https://img.shields.io/badge/Documentation-link-green.svg)](https://cds-astro.github.io/mocpy/)

<img align="right" width="200px" alt="mocpy's logo" src="./docs/_static/MOCpy-light.svg">

MOCPy is a Python library allowing easy creation and manipulation of MOCs (Multi-Order Coverage maps).

MOC is an [IVOA standard](https://ivoa.net/documents/MOC/20220727/index.html) enabling
description of arbitrary sky, time, (or frequency) coverages.

## What's in MOCpy?

|           | Space              | Time               | Frequency          |
|-----------|--------------------|--------------------|--------------------|
| Space     | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| Time      | :white_check_mark: | :white_check_mark: | :x:                |
| Frequency | :white_check_mark: | :x:                | :white_check_mark: |

There is support for one dimensional MOCs:

- **Space MOCs**. The sky tesselation is based on the HEALPix sky tessellation, it maps
  regions of the sky into hierarchically grouped cells.
- **Time MOCs**. The segmentation of the time axis is linear. This is described in the
  MOC standard since its version 2.0.
- **Frequency MOCs**. :warning: this is not standard yet! The segmentation of the frequency
  axis is logarithmic. See [this IVOA presentation][https://wiki.ivoa.net/internal/IVOA/InterOpMay2023Apps/FMOC_IVOA_Bologne_Fernique.pdf] for a description.

But also these combinations:

- **Space-Time MOCs**
- **Space-Frequency MOCs**

![a MOC plotted with MOCPy](./resources/readme.png)

*Rendered with MOCpy!*

## Other resources

For a command line tool, see the [moc-cli](https://github.com/cds-astro/cds-moc-rust/tree/main/crates/cli).

For more information about MOCPy's Rust core, see the [moc crate](https://crates.io/crates/moc)

## Installation

We strongly recommend to work in an environnement

### Latest stable version

- from pip `pip install mocpy`
- from conda `conda install -c conda-forge mocpy`

### Unreleased latest version

To install `MOCPy` from this repository, make sure you have Python (https://www.python.org/downloads/)
and Rust (https://www.rust-lang.org/tools/install) on your machine. Then, run

```
pip install git+https://github.com/cds-astro/mocpy.git
```

## To run the notebooks

The example notebooks require additional dependencies. They can be installed with

```sh
pip install mocpy[notebooks]
```

## Development installation

Contributions are very welcome. To build and test the library, clone this repository,
make sure that you've got Rust on your machine (https://www.rust-lang.org/tools/install)
and do:

```sh
pip install pre-commit
pre-commit install
pip install maturin
maturin develop --release
```

And the library is ready to be edited. After a change on the Rust side, you'll have to
do `maturin develop --release` again. If you only edit the python side, then a
`pip install .` at the root of the repository will suffice.

We use pre-commit to keep consistent style conventions between all of us. This means
that your commits will not pass if they require modifications. If this happens, do the
required changes, `git add` the files again and attempt to commit as usual.

To run the tests, do:

```sh
pip install .[dev]
python -m pytest
```

## For use in pyodide

Wheels that run in pyodide can be downloaded from each github's release for the module.

## Other python packages for MOCs

- https://github.com/grahambell/pymoc/