#!/bin/bash
# Abort if any simple command returns a non-zero value.
# See https://stackoverflow.com/questions/821396/aborting-a-shell-script-if-any-command-returns-a-non-zero-value
set -e

$PYTHON -m pytest -v mocpy --cov-report=term --cov=mocpy
# Test notebooks
$PIP install .
# Install other dependencies for running the notebooks
$PIP install jupyter astroquery regions
$PYTHON test_notebooks.py