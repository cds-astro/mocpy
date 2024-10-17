*****
MOCPy
*****
|PyPI version| |Build/Test status| |Notebook Binder| |Doc|


.. figure:: ./docs/_static/MOCpy-light.svg
   :width: 200 px
   :alt: mocpy's logo

MOCPy is a Python library allowing easy creation and manipulation of MOCs (Multi-Order Coverage maps).

MOC is an IVOA standard  enabling description of arbitrary sky regions.
Based on the HEALPix sky tessellation, it maps regions on the sky
into hierarchically grouped predefined cells.

A support for TMOC (temporal MOC) has been added since version 0.4.0.
It allows creation, parsing and comparison of TMOCs.
The support for Frequency MOCs has been added since version 0.13.0

Space & Time coverages (STMOC) are an extension of MOC to add time information.
It is possible to get a TMOC by querying a STMOC with a MOC and/or get a MOC
by querying a STMOC with a TMOC.

Please check the mocpy's `documentation <https://cds-astro.github.io/mocpy/>`__
for more details about installing MOCPy and using it.

For a command line tool, see the `moc-cli <https://github.com/cds-astro/cds-moc-rust/tree/main/crates/cli>`__.

For more information about the MOCPy Rust core, see the `moc crate <https://crates.io/crates/moc>`__.

.. figure:: ./resources/readme.png
   :width: 500 px
   :align: center
   :alt: a moc plotted with mocpy

   *Rendered with MOCpy!*

.. |PyPI version| image:: https://badge.fury.io/py/mocpy.svg
    :target: https://badge.fury.io/py/MOCPy

.. |Build/Test status| image:: https://github.com/cds-astro/mocpy/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cds-astro/mocpy/actions/workflows/test.yml

.. |Notebook Binder| image:: http://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/cds-astro/mocpy/master

.. |Doc| image:: https://img.shields.io/badge/Documentation-link-green.svg
    :target: https://cds-astro.github.io/mocpy/

Installation
------------

We strongly recommend to work in an environnement

Latest stable version
*********************

- from pip ``pip install mocpy``
- from conda ``conda install -c conda-forge mocpy``
- from this repository

Unreleased latest version
*************************

To install ``mocpy`` from this repository, make sure you have Python (https://www.python.org/downloads/)
and Rust (https://www.rust-lang.org/tools/install) on your machine. Then, run

.. code::

   pip install git+https://github.com/cds-astro/mocpy.git

To run the notebooks
********************

The example notebooks require additional dependencies. They can be installed with

.. code::

    pip install mocpy[notebooks]

For use in pyodide
******************

Wheels that run in pyodide can be downloaded from `this repository assets <https://github.com/cds-astro/mocpy/releases/download/v0.12.3/mocpy-0.12.3-cp310-cp310-emscripten_3_1_27_wasm32.whl>`__. This is not fully tested.
