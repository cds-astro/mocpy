#!/bin/bash

set -e

# install python and curl for downloading the rust
# compiler
apt update
apt install -y python3 python3-pip git curl

# configure the remote and permissions
# to push to the cds-astro/gh-pages branch
git config --global user.email "$GH_EMAIL"
git config --global user.name "$GH_NAME"

git remote rm origin
git remote add cds-astro https://"$GH_NAME":"$GH_TOKEN"@github.com/cds-astro/mocpy.git
git fetch cds-astro

git checkout gh-pages
rm -rf *

git checkout master .

curl https://sh.rustup.rs -sSf | sh -s -- --default-toolchain nightly -y
export PATH="$HOME/.cargo/bin:$PATH"

ln -s /usr/bin/python3 /usr/bin/python

python -m pip install -U pip
python -m pip install -r requirements/docs.txt

python setup.py build_rust
find build/ -name "*.so" -type f -exec cp {} ./mocpy \;

cd docs
make html
cd ..
mv docs/_build/html/ /tmp/
rm -rf *
mv /tmp/html/* .

touch .nojekyll
git add --all
git commit -am "doc update"
git push cds-astro gh-pages
