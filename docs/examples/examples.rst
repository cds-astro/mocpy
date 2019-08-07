Examples
========

Notebooks examples can be found on:

* `The official github repo <https://github.com/cds-astro/mocpy/tree/master/notebooks>`__.
* `Binder <https://mybinder.org/v2/gh/cds-astro/mocpy/master>`__. This allows you to run and modify the notebooks.

Spatial coverages
-----------------

Here are some code examples manipulating :py:class:`MOC` objects.

Loading and plotting the MOC of SDSS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example:

* Load the coverage of SDSS from a FITS file.
* Plot the MOC by:

1. Defining a matplotlib figure.
2. Defining an astropy WCS representing the field of view of the plot.
3. Call the :py:meth:`mocpy.moc.MOC.fill` and :py:meth:`mocpy.moc.MOC.border` so that mocpy plot on a matplotlib axis.
4. Set the axis labels, a title, enable the grid and plot the final figure.


.. plot:: examples/plot_SDSS_r.py
    :include-source:

Intersection between GALEX and SDSS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example:

* Load the coverages of SDSS and GALEX from FITS files.
* Compute their intersection
* Compute their union
* Plot the resulting intersection and union on a same matplotlib axis.


.. plot:: examples/logical_operations.py
    :include-source:

Create a MOC from a concave polygon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example shows how to call :py:meth:`mocpy.moc.MOC.from_polygon` or :py:meth:`mocpy.moc.MOC.from_polygon_skycoord`.

.. plot:: examples/polygon_coverage.py
    :include-source:

Get the border(s) of a MOC
~~~~~~~~~~~~~~~~~~~~~~~~~~

This example shows how to call :py:meth:`mocpy.moc.MOC.get_boundaries`. The borders are returned as a list of `~astropy.coordinates.SkyCoord` each defining one border.
In this example:

1. The sky coordinates defining the border of the MOC are projected to the pixel image system.
2. Then, a matplotlib path is defined from the projected vertices.
3. Finally the path is plot on top of the MOC.

.. plot:: examples/compute_borders.py
    :include-source:

Temporal coverages
------------------

A class :py:class:`TimeMOC` describes temporal coverages. 

Please refer to the following notebook `here <https://github.com/cds-astro/mocpy/blob/ranges2D/notebooks/tmoc.ipynb>`__ for how to use it.

Space & Time coverages
----------------------

Space-Time coverages are a new feature of `mocpy` since its version 0.7.0 and are
an attempt initiated by the Virtual Observatory for binding spatial and temporal coverages together.
See its description formalized by the IVOA `here <http://www.ivoa.net/documents/stmoc/20190515/NOTE-stmoc-1.0-20190515.pdf>`__.

Space-Time coverages allows you to:

1. Retrieve the spatial coverage observed by a mission within a set of time frames (i.e. `~astropy.time.Time` ranges).
2. Retrieve the temporal coverage observed by a mission within a spatial coverage.

As we do for spatial or temporal coverages, one can also perform the union, intersection or difference between two Space-Time coverages.

Please refer to the following notebook `here <https://github.com/cds-astro/mocpy/blob/ranges2D/notebooks/Space%20%26%20Time%20coverages.ipynb>`__ for how to compute and query Space-Time coverages.