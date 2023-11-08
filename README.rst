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

An experimental support for TMOC (temporal MOC) has been added since version 0.4.0.
It allows creation, parsing and comparison of TMOCs.

Space & Time coverages (STMOC) are an extension of MOC to add time information.
It is possible to get a TMOC by querying a STMOC with a MOC and/or get a MOC
by querying a STMOC with a TMOC.

Please check the mocpy's `documentation <https://cds-astro.github.io/mocpy/>`__
for more details about installing MOCPy and using it.

For a command line tool, see the `moc-cli <https://github.com/cds-astro/cds-moc-rust/tree/main/crates/cli>`__.

For more information about the MOCPy Rust core, see the `moc crate <https://crates.io/crates/moc>`__.

.. figure:: ./resources/readme.png
   :scale: 50 %
   :align: center
   :alt: map to buried treasure

   *Rendered with MOCpy!*

.. |PyPI version| image:: https://badge.fury.io/py/mocpy.svg
    :target: https://badge.fury.io/py/MOCPy

.. |Build/Test status| image:: https://github.com/cds-astro/mocpy/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cds-astro/mocpy/actions/workflows/test.yml

.. |Notebook Binder| image:: http://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/cds-astro/mocpy/master

.. |Doc| image:: https://img.shields.io/badge/Documentation-link-green.svg
    :target: https://cds-astro.github.io/mocpy/

Migrating to version 0.12
-------------------------

Since 0.12.3
************

- ``MOC.MAX_ORDER` and `TimeMOC.MAX_ORDER`` replace the former ``IntervalSet.HPX_MAX_ORDER`` and ``IntervalSet.TIME_MAX_ORDER``
- ``MOC.to_depth29_ranges`` is now a public method replacing the former private ``IntervalSet.nested`` and addition of ``TimeMOC.to_depth61_ranges`` for a time counterpart

Since v0.12.0
*************

- ``MOC.contains_skycoords`` and ``MOC.contains_lonlat`` replace ``MOC.contains`` (``contains`` will be removed in v1.0.0)
- ``TimeMOC.contains_with_timeresolution`` has been added with the previous behaviour of  ``TimeMOC.contains``
- ``from_uniq` removed from `IntervalSet`` and added to ``MOC``
- ``MOC.from_healpix_cells`` now requires the ``max_depth`` argument, the depth of the MOC we want to create
- ``World2ScreenMPL`` has been renamed ``WCS``

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
.. code::

   git clone https://github.com/cds-astro/mocpy.git
   cd mocpy
   pip install .

Note that the point is important.

To run the notebooks
********************

The example notebooks require additional dependencies. They can be installed with

.. code::

    pip install mocpy[notebooks]

For use in pyodide
******************

Wheels that run in pyodide can be downloaded from `this repository assets <https://github.com/cds-astro/mocpy/releases/download/v0.12.3/mocpy-0.12.3-cp310-cp310-emscripten_3_1_27_wasm32.whl>`__. This is not fully tested.
