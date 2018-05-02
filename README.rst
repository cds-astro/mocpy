*****
MOCPy
*****

MOCPy is a Python library allowing easy creation, parsing and manipulation of MOCs (Multi-Order Coverage maps).

============
Introduction
============

MOC is an IVOA standard ( http://ivoa.net/documents/MOC/ ) enabling description 
of arbitrary sky regions. Based on the HEALPix sky tessellation, it maps 
regions on the sky into hierarchically grouped predefined cells.

=======
License
=======

MOCPy is distributed under GPLv3 license.

============
Requirements
============

``numpy``, ``astropy_healpix`` and ``astropy`` are required.
``matplotlib`` is also needed if you want to plot MOC objects.

MOCPy runs under Python 2 and Python 3.

============
Installation
============

``pip install mocpy``

===========
Using MOCPy
===========

-------------
Reading a MOC
-------------

===========
Examples
===========

Different examples on how to use MOCPy can be found in the notebooks directory.

.. image:: http://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/cds-astro/mocpy/tmoc

Launch the Binder to execute live and interact with the notebooks examples.  


=====
Tests
=====

To run the automated tests of ``mocpy`` you can use ``pytest`` and this command:

    python -m pytest -v mocpy

Continuous integration test status:

* .. image:: http://img.shields.io/travis/cds-astro/mocpy.svg?branch=master
    :target: https://travis-ci.org/cds-astro/mocpy
    :alt: Test Status Travis-CI
