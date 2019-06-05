#!/bin/bash
# Abort if any simple command returns a non-zero value.
# See https://stackoverflow.com/questions/821396/aborting-a-shell-script-if-any-command-returns-a-non-zero-value
set -e

python -m pytest -v mocpy --cov-report=term --cov=mocpy
# Test notebooks
pip install .
# Install other dependencies for running the notebooks
pip install jupyter astroquery
python test_notebooks.py