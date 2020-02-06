#!/bin/sh
#
# Build and test mocpy in a 32-bit environment with python3
#
# Usage:
#   testing_py3_ubuntu32.sh
#

# Update packages to the latest available versions
linux32 --32bit i386 sh -c '
    apt update > /dev/null &&
    apt install -y curl pkg-config libfreetype6-dev \
	python3 python3-pip >/dev/null &&
    # Symbolic link so that pyO3 points to the
    # python3 version
    ln -s /usr/bin/python3 /usr/bin/python
' &&

# Run the tests
linux32 --32bit i386 sh -c '
    pip3 install -U pip &&
    # Download the dependencies for compiling cdshealpix
    pip install -r requirements/contributing.txt &&
    pip install setuptools_rust &&
    # Install Rust compiler
    curl --proto "=https" --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- --default-toolchain nightly -y &&
    export PATH="$HOME/.cargo/bin:$PATH" &&
    # Build the rust dynamic library
    python3 setup.py build_rust &&
    # Move the dynamic lib to the python package folder
    find build/ -name "*.so" -type f -exec cp {} ./mocpy \; &&
    python3 -m pytest -v mocpy
'
