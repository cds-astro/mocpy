##################
User documentation
##################

***********************************
SMOC (a.k.a MOC): Spatial coverages
***********************************

Introduction to Space-MOCs
==========================

Space MOCs are an ensemble of HEALPix cells of mixed orders. They represent sky regions
in an approximate way. They are designed for efficient calculations between sky regions.

The coordinates for Space-MOCs are always in IRCS at epoch J2000 by definition.

Space-MOCs can represent arbitrary shapes. Common examples are an approximated cone, an
ensemble of approximated cones, or the coverage of a specific survey.

The following notebook is an introduction to manipulation of Space-MOCs with MOCPy:

.. nbgallery::
    ../_collections/notebooks/00-MOCpy_introduction

As you saw, Space-MOCs can be visualized with either astropy+matplotlib or with the
ipyaladin widget.

Space-MOC creation methods
==========================

Reading a MOC from a file or from a server
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are diverse ways to instantiate a MOC. A first approach is to get a MOC from a
FITS file on our disk, or from a distant server that already has a pre-calculated MOC.

.. nbgallery::
    ../_collections/notebooks/from_fits_and_intersection

Calculating a MOC on the fly from a region of the sky
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a simple example of the creation of a MOC that fully covers a polygon defined
by its vertices on the sky. The vertices are linked by great circles.

.. plot:: examples/polygon_coverage.py
    :include-source:

For a more extended description on how to create MOCs from diverse shapes, you can check
the example notebooks :doc:`../_collections/notebooks/01-Creating_MOCs_from_shapes` and
:doc:`../_collections/notebooks/02-Creating_MOCs_from_astropy_regions`.

.. nbgallery::
    ../_collections/notebooks/01-Creating_MOCs_from_shapes
    ../_collections/notebooks/02-Creating_MOCs_from_astropy_regions
    ../_collections/notebooks/from_astropy_table.ipynb
    ../_collections/notebooks/from_image_with_mask

Instantiating a MOC when we already know its HEALPix cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This can be done with methods like :py:meth:`~mocpy.moc.MOC.from_string`,
:py:meth:`~mocpy.moc.MOC.from_healpix_cells`, :py:meth:`~mocpy.moc.MOC.from_depth29_ranges`.
A simple MOC that represents the full sky is the MOC of order 0 that contains all the 12
HEALPix cells of the order 0. Creating it from a string looks like:

.. plot:: examples/all_sky.py
    :include-source:

We could also take only odd numbered cells of order 1. Let's use
:py:meth:`~mocpy.moc.MOC.from_json` that takes a dictionary:

.. plot:: examples/odd_cells.py
    :include-source:

Useful methods
==============

Calculating a Space-MOC sky area
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example shows how to calculate the sky area of a Space-MOC.

.. plot:: examples/calculate_moc_sky_area.py
    :include-source:

Other examples (operations, use-cases)
======================================

Here are some other code examples manipulating :py:class:`MOC` objects.

.. nbgallery::
    ../_collections/notebooks/compute_moc_borders
    ../_collections/notebooks/filtering_astropy_table
    ../_collections/notebooks/FITS-image-pixels-intersecting-MOC
    ../_collections/notebooks/query_vizier_table

Loading and plotting the MOC of SDSS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example:

* Load the coverage of SDSS from a FITS file.
* Plot the MOC by:

1. Defining a matplotlib figure.
2. Defining an astropy WCS representing the field of view of the plot.
3. Call the :py:meth:`~mocpy.moc.MOC.fill` and :py:meth:`~mocpy.moc.MOC.border` so that mocpy plot on a matplotlib axis.
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

Get the border(s) of a MOC
~~~~~~~~~~~~~~~~~~~~~~~~~~

This example shows how to call `~mocpy.moc.MOC.get_boundaries`. The borders are returned as a list of `~astropy.coordinates.SkyCoord` each defining one border.
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

Please check the following notebooks if you want to know more about time MOCs.


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

Please check the following notebooks if you want to know more about space-time MOCs.

.. nbgallery::
    ../_collections/notebooks/STMOC from time ranges
    ../_collections/notebooks/Space & Time coverages

*************************
FMOC: Frequency coverages
*************************

Please check the following notebooks if you want to know more about frequency MOCs.


.. nbgallery::
    ../_collections/notebooks/First_Steps_with_FMOCs

********************************
SFMOC: Space-Frequency coverages
********************************

Creating SF-MOC instances
=========================

1. From a FITS, ASCII, JSON file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `SFMOC` class allows manipulation and generation of
**Space-Frequency Multi-Order Coverages** (SF-MOCs).

New empty SF-MOC `~mocpy.SFMOC.new_empty`:

    >>> from mocpy import SFMOC
    >>> SFMOC.new_empty(max_order_frequency=40, max_order_space=12)
    f40/ s12/

Loading an SF-MOC from a FITS, JSON or ASCII file `~mocpy.SFMOC.load`:

    >>> # we generate a fake ASCII file
    >>> from tempfile import NamedTemporaryFile
    >>> with NamedTemporaryFile(delete_on_close=False) as tf:
    ...     tf.write(b"f20/0-100 s12/0-100")
    ...     tf.close()
    ...     # and we read it
    ...     SFMOC.load(tf.name, format="ascii")
    19
    f14/0
    15/2
    18/24
    20/100
    s9/0
    10/4-5
    11/24
    12/100
    f20/ s12/

This illustrated the load method that we use when we have a MOC in a file, but for
ASCII strings, one can also directly use `~mocpy.SFMOC.from_string`:

    >>> SFMOC.from_string("f15/2 s12/1")
    f15/2
    s12/1
    f15/ s12/

Creating SF-MOCs from physical parameters:

2. From frequency-position "points"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In `~mocpy.SFMOC.from_frequencies_and_positions`, we give a set of frequency-position
"points" and the SF-MOC will be built from the cells containing these punctual data at
the maximum requested order. Use the utility functions
`~mocpy.MOC.spatial_resolution_to_order` to chose the space resolution and
`~mocpy.FrequencyMOC.relative_precision_to_order` to chose the frequency resolution.

    >>> import astropy.units as u
    >>> from mocpy import SFMOC
    >>> SFMOC.from_frequencies_and_positions(
    ...     [0.01, 0.1, 1, 10, 100] * u.Hz,
    ...     lon=[0, 1, 2, 3, 4] * u.deg, lat=[0, 1, 2, 3, 4] * u.deg,
    ...     max_order_frequency=25, max_order_space=8)
    f25/22879928
    s8/311296
    f25/23750246
    s8/311316
    f25/24641536
    s8/311378
    f25/25493504
    s8/311558
    f25/26361856
    s8/311624
    f25/ s8/

3. From frequency ranges associated to punctual positions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With `mocpy.SFMOC.from_frequency_ranges_and_positions` we do not get a single cell in
frequency, bu the cells covering a given frequency band/range:

    >>> import astropy.units as u
    >>> from mocpy import SFMOC
    >>> SFMOC.from_frequency_ranges_and_positions(
    ...     frequencies_min=[0.01, 0.1, 1, 10, 100] * u.Hz,
    ...     frequencies_max=[0.1, 1, 10, 100, 1000] * u.Hz,
    ...     lon=[0, 1, 2, 3, 4] * u.deg, lat=[0, 1, 2, 3, 4] * u.deg,
    ...     max_order_frequency=6, max_order_space=8)
    f6/43-44
    s8/311296
    f6/45
    s8/311296 311316
    f6/46
    s8/311316
    f6/47
    s8/311378
    f6/48
    s8/311378 311558
    f6/49
    s8/311558
    f6/50
    s8/311558 311624
    f6/51
    s8/311624
    f6/ s8/

4. From frequency ranges associated to Space-MOCs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In `~mocpy.SFMOC.from_spatial_coverages`, the positions are not punctual anymore. They
are represented by S-MOCs, which can be built with the `~mocpy.MOC` module. Let's for
example build an SF-MOC from a square observation in two given frequency bands:

    >>> from mocpy import MOC, SFMOC
    >>> import astropy.units as u
    >>> # a MOC representing a field of view
    >>> moc = MOC.from_box(lon=0*u.deg, lat=0*u.deg, a=2*u.deg, b=2*u.deg,
    ...                    angle=0*u.deg, max_depth=10)
    >>> # we use the same moc twice, and give two ranges
    >>> sfmoc_fov = SFMOC.from_spatial_coverages(
    ...    frequencies_min=([200, 300]*u.nm).to(u.Hz, equivalencies=u.spectral()),
    ...    frequencies_max=([100, 200]*u.nm).to(u.Hz, equivalencies=u.spectral()),
    ...    spatial_coverages=[moc, moc],
    ...    max_order_frequency=13)

Utilities
=========

To inspect an `~mocpy.SFMOC` instance, one can read its `~mocpy.SFMOC.max_order`,
`~mocpy.SFMOC.min_frequency` and `~mocpy.SFMOC.max_frequency` properties. We can also
check wether it is empty with `~mocpy.SFMOC.is_empty`. Let's use ``sfmoc_fov`` from the
previous example:

    >>> sfmoc_fov.max_order
    (13, 10)
    >>> sfmoc_fov.is_empty()
    False
    >>> sfmoc_fov.min_frequency
    <Quantity 9.93958512e+14 Hz>
    >>> sfmoc_fov.max_frequency
    <Quantity 3.025856e+15 Hz>


General manipulation
====================

Check wether a list of "punctual" position/frequency values fall within an SF-MOC:

    >>> from mocpy import SFMOC, MOC
    >>> import astropy.units as u
    >>> # we generate an SF-MOC from a cone and a range of frequencies
    >>> moc = MOC.from_cone(0*u.deg, 0*u.deg, radius=10*u.deg, max_depth=10)
    >>> sfmoc_cone = SFMOC.from_spatial_coverages(0.01*u.Hz, 100*u.Hz,
    ...                                           moc, max_order_frequency=40)
    >>> # one inside, one outside
    >>> sfmoc_cone.contains([10, 10000]*u.Hz, [0.1, 20]*u.deg, [0.1, 20]*u.deg)
    array([ True, False])

Fold operation on one of the SF-MOC dimensions, use `mocpy.SFMOC.query_by_frequency` or
`~mocpy.SFMOC.query_by_space`

    >>> from mocpy import FrequencyMOC as FMOC
    >>> fmoc = FMOC.from_frequency_ranges(15, 0.1*u.Hz, 1*u.Hz)
    >>> moc = sfmoc_cone.query_by_frequency(fmoc)
    >>> type(moc)
    <class 'mocpy.moc.moc.MOC'>

Here we recover the Space-MOC corresponding to the parts of the ``sfmoc_cone`` coverage
with frequencies within ``fmoc``.
