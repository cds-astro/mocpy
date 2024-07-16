##################
User documentation
##################

***********************************
SMOC (a.k.a MOC): Spatial coverages
***********************************

Gallery of notebooks examples using SMOCs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. nbgallery::
    ../_collections/notebooks/00-MOCpy_introduction
    ../_collections/notebooks/compute_moc_borders
    ../_collections/notebooks/filtering_astropy_table
    ../_collections/notebooks/FITS-image-pixels-intersecting-MOC
    ../_collections/notebooks/from_astropy_table.ipynb
    ../_collections/notebooks/01-Creating_MOCs_from_shapes
    ../_collections/notebooks/02-Creating_MOCs_from_astropy_regions
    ../_collections/notebooks/from_fits_and_intersection
    ../_collections/notebooks/from_image_with_mask
    ../_collections/notebooks/from_vizier_table
    ../_collections/notebooks/query_vizier_table

Here are some code examples manipulating :py:class:`MOC` objects.

Examples use cases
==================


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

Create a MOC from a polygon
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example shows how to call :py:meth:`mocpy.moc.MOC.from_polygon` or :py:meth:`mocpy.moc.MOC.from_polygon_skycoord`.

.. plot:: examples/polygon_coverage.py
    :include-source:

For a more extended description on how to create MOCs from shapes, you can check the example notebooks
:doc:`../_collections/notebooks/01-Creating_MOCs_from_shapes` and
:doc:`../_collections/notebooks/02-Creating_MOCs_from_astropy_regions`.

Get the border(s) of a MOC
~~~~~~~~~~~~~~~~~~~~~~~~~~

This example shows how to call `mocpy.moc.MOC.get_boundaries`. The borders are returned as a list of `~astropy.coordinates.SkyCoord` each defining one border.
In this example:

1. The sky coordinates defining the border of the MOC are projected to the pixel image system.
2. Then, a matplotlib path is defined from the projected vertices.
3. Finally the path is plot on top of the MOC.

.. plot:: examples/compute_borders.py
    :include-source:

Gravitational Waves MOCs
~~~~~~~~~~~~~~~~~~~~~~~~

This example shows the probability confidence regions of gravitational waves.
HEALPix cells are given under the
`uniq pixel notation <http://www.ivoa.net/documents/Notes/MOC/20120412/NOTE-MOC-1.0-20120412.pdf>`__.
Each pixel is associated with a specific probability density value. We convert this into
a probability by multiplying it with the area of each cell.
Then, we can create a MOC from which a GW has x% of chance of being localized in it.
By definition the MOC which has 100% of chance of containing a GW is the full sky MOC.

.. plot:: examples/bayestar.py
    :include-source:

Calculate MOC sky area
~~~~~~~~~~~~~~~~~~~~~~

This example shows how to Calculate the sky area of a MOC instance.

.. plot:: examples/calculate_moc_sky_area.py
    :include-source:

Performing computation on the pixels of an FITS image lying in a MOC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example shows how a MOC can filter pixels from a specific FITS image (i.e. associated with a WCS). These pixels can
then be retrieved from the image for performing some computations on them: e.g. mean, variance analysis thanks to numpy/scikit-learn...

.. plot:: examples/filter_image_pixels.py
    :include-source:

************************
TMOC: Temporal coverages
************************

The :py:class:`TimeMOC` class represents a temporal coverage.

Gallery of notebooks examples using TMOCs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. nbgallery::
    ../_collections/notebooks/tmoc


*****************************
STMOC: Space & Time coverages
*****************************

Space-Time coverages are a new feature of ``mocpy`` since its version 0.7.0 and are bind spatial and temporal coverages together.
The standard description is published by the IVOA `here <http://www.ivoa.net/documents/stmoc/20190515/NOTE-stmoc-1.0-20190515.pdf>`__.

Space-Time coverages allow to:

1. Retrieve the spatial coverage observed by a mission within a set of time frames (i.e. `astropy.time.Time` ranges).
2. Retrieve the temporal coverage observed by a mission within a spatial coverage.

As we do for spatial or temporal coverages, one can also perform the union, intersection or difference between two Space-Time coverages.

Gallery of notebooks examples using STMOCs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. nbgallery::
    ../_collections/notebooks/STMOC from time ranges
    ../_collections/notebooks/Space & Time coverages

*************************
FMOC: Frequency coverages
*************************

Gallery of notebooks examples using FMOCs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. nbgallery::
    ../_collections/notebooks/First_Steps_with_FMOCs

