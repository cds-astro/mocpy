# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.3] - 2019-01-30
### Added
- A new constructor from_cells. It returns a new MOC instance from a numpy structured array representing the cells of a MOC. The structured array passed must contains 3 attributes: "ipix": np.uint64, "depth": np.uint32, "fully_covered": np.uint8 (a flag bit. For the moment its value has no effect for the newly created MOC).

## [0.5.2] - 2019-01-28
### Added
- A new **from_polygon_skycoord** method added where you can pass an astropy.coordinates.SkyCoord describing the polygon coordinates instead of two lon, lat astropy.units.Quantity. The *max_depth*, and *inside* optional arguments remain.

### Changed
- Remove spherical geom from dependency so that astroquery.cds wheel for windows/py3 can be generated. Spherical Geom is only used in MOC.from_polygon. A message is addressed to the user telling him to install sphrical geom if it is not installed and if he wants to create a MOC from a polygon.

## [0.5.1] - 2019-01-25
### Changed
- `pip install mocpy` now installs all the dependencies of mocpy. See the setup.py file. (requires changed to install_requires).

## [0.5.0] - 2019-01-09
### Added
- Two methods `fill` and `border` taking an astropy.wcs.WCS and a matplotlib axis. `fill` projects the MOC into the WCS and draws it on the MPL axis using pathpatches for each HEALPix cell. `border` draws the border the same way and requires a WCS and an MPL axe too. You can pass to these functions additional keyword arguments that will directly be passed to MPL when plotting (e.g. the color, the linewidth, and alpha component etc...). Check the notebooks to see how to use these new methods.
- You can retrieve the list of skycoords describing the border MOC. Each skycoord from the list refers to one border of the MOC (either an external border or the border of a hole in a connexe MOC). The process takes for the moment a quite amount of time and thus may be optimized in the future. A new GALEX boundaries notebook has been added so that you can check how it works. I recommend to decrease the order of GALEX to 5 before computing the boundaries otherwise it will take some time. This add relies on the issue [#29](https://github.com/cds-astro/mocpy/issues/29) initiated by [@ggreco77](https://github.com/ggreco77).
- A new `from_polygon` method taking the vertices (i.e. (skycoords) or (lon, lat) tuples) responsible for setting up a new MOC from a polygon.
An inside SkyCoord point is requested and says to the algorithm which area (2 possible as on the sphere, an area and its complement are both finite) must be chosen. If no inside sky coord is given, we consider the mean of the vertices of the polygon as belonging to the MOC (This is without ambiguity for convex polygons but it may not work for concave ones). Vertices describing a convex or concave polygon are accepted. Vertices describing a self-intersecting polygon are not accepted.
This method does not rely on astropy_healpix as there is for the moment no method returning the set of HEALPix cells belonging to a polygon and is thus implemented purely in Python.
- A new `serialize` public method allows to serialize the MOC in two possible format, FITS and json. For a FITS serialization the method returns an astropy HDU binary table. For a JSON serialization, the method returns a python dictionary storing order-[list of HEALPix cells] as key-value pairs.

### Changed
- `write` method does not take a `write_to_file` argument. For serialization purpose, there is now a new `serialize` method.
- astropy_healpix.HEALPix.lonlat_to_healpix seems to not accept astropy MaskedColumn anymore. For lon as a MaskedColumn, please pass lon.T * u.deg to mocpy.MOC.from_lonlat. We need to transpose the column and then convert it to an astropy.units.Quantity.
- Notebooks have been updated and all the plots now use the new methods `fill` and `border`.
- A new package `spatial`, invisible from the user, but keeping all the code of spatial   MOCs (plotting methods, border computation, special utils for creating WCS...) has been created. tmocs and core functions are still located in the root of the project.
- Add of a changelog
