Install
=======

Stable versions
---------------

MOCPy relies on the following mandatory external dependencies :

- `cdshealpix <https://cds-astro.github.io/cds-healpix-python/>`__
- `astropy <http://docs.astropy.org/en/stable/>`__
- `matplotlib <https://matplotlib.org/>`__
- `networkx <http://networkx.github.io/>`__ used in ``~mocpy.MOC.get_boundaries``

It has an optional dependency to query FITS files by url:

- `requests <https://github.com/psf/requests>`__

To upgrade the ``mocpy`` package to the latest version::

    pip install mocpy -U

To install ``mocpy`` and all mandatory dependencies type::

    pip install mocpy

The documentation contains several examples and notebooks that make use of additional
optional dependencies. If you want to run these notebooks, you can install ``mocpy`` and all
mandatory and optional dependencies for the notebooks with::

    pip install "mocpy[notebooks]"

From source
-----------

To install MOCPy from source, you'll need `Rust <https://www.rust-lang.org/tools/install>`_.

Then you can download the source code from `this link <https://github.com/cds-astro/mocpy/archive/refs/heads/master.zip>`_
or directly from the GitHub repository. Once in the source code folder, you'll have to do

::
    > pip install maturin
    > maturin develop --release
    > pip install .

And you'll have the very last changes running on your machine.