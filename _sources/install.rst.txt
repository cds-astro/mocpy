Install
=======

MOCPy relies on the following mandatory external dependencies :

- `cdshealpix <https://cds-astro.github.io/cds-healpix-python/>`__
- `astropy <http://docs.astropy.org/en/stable/>`__
- `matplotlib <https://matplotlib.org/>`__
- `networkx <http://networkx.github.io/>`__ used in `~mocpy.MOC.get_boundaries`
- `lark-parser <https://github.com/lark-parser/lark>`__ used for
  serializing the sky coverages to ASCII format

It has also a few optional dependencies that are mainly used in
deprecated methods and will be removed in the future:

- `astropy-healpix <http://astropy-healpix.readthedocs.io/en/latest/>`__
  is used by the old plotting methods as well as in `~mocpy.MOC.get_boundaries`

To upgrade the ``mocpy`` package to the latest version::

    pip install mocpy -U

To install ``mocpy`` type::

    pip install mocpy
