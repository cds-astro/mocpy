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
As of the version 0.6.0, MOCPy is only distributed for Python 3.5/3.6
and 3.7. Python 2 users can still refer to the
`0.5.7 <https://github.com/cds-astro/mocpy/releases/tag/v0.5.7>`__ version.
Its code is hosted on `GitHub <https://github.com/cds-astro/mocpy/>`__ and
distributed under the BSD-3 license.

MOC is an `IVOA standard <http://ivoa.net/documents/MOC/>`__ enabling description
of arbitrary sky regions. Based on the HEALPix sky tessellation, it maps
regions on the sky into hierarchically grouped predefined cells.

MOCPy provides the `~mocpy.MOC` and `~mocpy.TimeMOC` classes handling
respectively the manipulation of spatial and temporal MOCs.

As an example, here is the sky coverage of the SDSS sky survey:

.. plot:: examples/plot_SDSS_r.py

As well as its time coverage:

.. plot:: examples/plot_TMOC_SDSS_r.py

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
