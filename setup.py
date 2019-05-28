#!/usr/bin/env python

import sys

from setuptools import setup, find_packages
from setuptools_rust import RustExtension

PYTHON_MAJOR_VERSION = sys.version_info.major

exec(open('mocpy/version.py').read())

setup(name='MOCPy',
    packages=find_packages(),
    version=__version__,
    description='MOC parsing and manipulation in Python',
    author='Thomas Boch, Matthieu Baumann',
    author_email='thomas.boch@astro.unistra.fr',
    license='BSD',
    url='https://github.com/cds-astro/mocpy',
    install_requires=['astropy',
        'astropy_healpix',
        'cdshealpix',
        'matplotlib', # Used in fill and border
        'networkx', # Used in get_boundaries
        # Disable the installation of spherical-geometry as it
        # fails the appveyor tests for python 3 of astroquery. See https://github.com/astropy/astroquery/pull/1343
        # 'spherical-geometry',
        'lark-parser', # Used in from_str for parsing the string given and create the MOC from it
    ],
    provides=['mocpy'],
    long_description="MOCPy is a Python library allowing easy creation \
     and manipulation of MOCs (Multi-Order Coverage maps). \
     MOC is an `IVOA standard <http://ivoa.net/documents/MOC/>` \
     enabling description of arbitrary sky regions. \
     Based on the HEALPix sky tessellation, it maps regions on the sky \
     into hierarchically grouped predefined cells.\n \
     An experimental support for TMOC (temporal MOC) has been added since version 0.4.0.\
     It allows creation, parsing and comparison of TMOCs.",
    rust_extensions=[RustExtension(
        # Package name
        "mocpy.core",
        # The path to the cargo.toml file defining the rust-side wrapper.
        # This file usually contains the name of the project, its version, the author
        # and the dependencies of the crate (in our case the rust wrapper depends on the cdshealpix
        # crate). 
        'Cargo.toml',
        # Specify python version to setuptools_rust
        rustc_flags=['--cfg=Py_{}'.format(PYTHON_MAJOR_VERSION)],
        features=['numpy/python{}'.format(PYTHON_MAJOR_VERSION)],
        # Add the --release option when building the rust code
        debug=False)],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python',
        'License :: OSI Approved :: BSD License',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False,
)
