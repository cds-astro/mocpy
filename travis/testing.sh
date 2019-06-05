#!/bin/bash
# Abort if any simple command returns a non-zero value.
# See https://stackoverflow.com/questions/821396/aborting-a-shell-script-if-any-command-returns-a-non-zero-value
set -e

# Generate the dynamic library from the cdshealpix Rust crate.
# This will download the crate from crates.io and build it first.
python setup.py build_rust
# Move the dynamic lib to the python package folder
find build/ -name "*.so" -type f -exec cp {} ./mocpy \;
python -m pytest -v mocpy --cov-report=term --cov=mocpy
# Test notebooks
pip install .
# Install other dependencies for running the notebooks
pip install jupyter astroquery
python test_notebooks.py