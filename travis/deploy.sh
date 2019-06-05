#!/bin/bash
# Abort if any simple command returns a non-zero value.
# See https://stackoverflow.com/questions/821396/aborting-a-shell-script-if-any-command-returns-a-non-zero-value
set -e

### Build the wheels ###
$PIP install cibuildwheel==0.10.0 setuptools-rust
export CIBW_BEFORE_BUILD="pip install setuptools-rust && curl https://sh.rustup.rs -sSf | sh -s -- --default-toolchain nightly -y"
export CIBW_ENVIRONMENT='PATH="$HOME/.cargo/bin:$PATH"'
cibuildwheel --output-dir dist
### Upload the wheels to PyPI ###
# If the commit is tagged
if [[ $TRAVIS_TAG ]]; then
    $PIP install twine
    $PYTHON -m twine upload --repository-url https://upload.pypi.org/legacy/ dist/*.whl --skip-existing
fi
