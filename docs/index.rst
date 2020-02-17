..  _index:

Welcome to MOCPy's documentation!
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   examples/examples
   api
   contribute

MOCPy is a Python library allowing easy creation, parsing and manipulation
of MOCs (Multi-Order Coverage maps).

Its code is hosted on `GitHub <https://github.com/cds-astro/mocpy/>`__ and
distributed under the BSD-3 license.

What is a MOC ?
---------------

MOC is an `IVOA standard <http://ivoa.net/documents/MOC/>`__ enabling description
of arbitrary sky regions. Based on the HEALPix sky tessellation, it maps
regions on the sky into hierarchically grouped predefined cells.

MOCPy provides the `~mocpy.MOC` and `~mocpy.TimeMOC` classes handling
respectively the manipulation of spatial and temporal MOCs.

As an example, here is the sky coverage of the SDSS sky survey:

.. plot:: examples/plot_SDSS_r.py

As well as its time coverage:

.. plot:: examples/plot_TMOC_SDSS_r.py

- :cite:`2014ivoa.spec.0602F`

References
----------

.. bibliography:: references.bib
    :all:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
