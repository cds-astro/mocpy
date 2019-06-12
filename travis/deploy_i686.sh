#!/bin/sh
#
# Deploy mocpy for i686 archs with python3
#
# Usage:
#   deploy_i686.sh
#
# See https://stackoverflow.com/questions/821396/aborting-a-shell-script-if-any-command-returns-a-non-zero-value
set -e

### Build the wheels ###
$PIP install cibuildwheel==0.11.1 setuptools-rust
export CIBW_BEFORE_BUILD="yum install -y libpng-devel &&
 yum install -y wget &&
 wget https://download.savannah.gnu.org/releases/freetype/freetype-2.4.0.tar.gz &&
 tar -zxvf freetype-2.4.0.tar.gz &&
 cd freetype-2.4.0/ &&
 make install &&
 cd .. &&
 pip install setuptools-rust &&
 curl https://sh.rustup.rs -sSf | sh -s -- --default-toolchain nightly -y"
export CIBW_ENVIRONMENT='PATH="$HOME/.cargo/bin:$PATH"'
cibuildwheel --output-dir dist
### Upload the wheels to PyPI ###
# If the commit is tagged
if [[ $TRAVIS_TAG ]]; then
    $PIP install twine
    $PYTHON -m twine upload --repository-url https://upload.pypi.org/legacy/ dist/*.whl --skip-existing
fi