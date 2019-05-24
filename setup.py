#!/usr/bin/env python

from setuptools import setup, find_packages

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
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python',
        'License :: OSI Approved :: BSD License',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
)
