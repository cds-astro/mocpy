:ref:`index`

Contribute
==========

MOCPy is a wrapper around `cds-moc-rust <https://github.com/cds-astro/cds-moc-rust>`__.
This is why this repository contains both python -- in the ``python`` folder -- and rust
files -- in the ``src`` folder.

We use `pyO3 bindings <https://pyo3.rs>`__.

Installation
------------

We recommend setting up a python virtual environment that will act independently from
your current installed interpreter and package manager. We suggest the `standard library
python module venv <https://docs.python.org/3/library/venv.html>`__ or
`uv venv <https://docs.astral.sh/uv/pip/environments/>` if you don't know where to look.

You can now install all the necessary packages for developing and testing MOCpy::

    pip install .[dev,plots,notebooks]
    pre-commit install
    maturin develop --release

Then, if you edit the Rust part of the library you'll have to run again the
``maturin develop --release`` and the ``pip install .`` command. If you only edit the
python part, then ``pip install .`` is enough.

If a ``maturin develop --release`` does not actualize your ``Cargo.toml`` file, try: ::

    rm Cargo.lock && cargo clean

And run the ``maturin`` command again.

These commands also install `git hooks
<https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks>`__` in the repository.
You will see the hooks run linting and formatting scripts at the pre-commit stage. They:

- trim trailing whitespaces automatically in any files
- run ``cargofmt`` to format the ``.rs`` files
- sort requirements  and run ``ruff format`` to format the ``.py`` files
- clear the metadata of the ``.ipynb`` files

If any of these hooks applied changes to the files you wanted to commit, you can
``git add`` them again and commit the new versions.


Running the python tests
------------------------

Once your environment is set up and activated you can run the tests

To run the automated tests and the docstring examples, go to the repo folder and type::

    python -m pytest -v python/mocpy

To run the tests with coverage report locally::

    python -m pytest -v python/mocpy --cov-report=term --cov=python/mocpy

You also can have a html output of the coverage with the flag ``--cov-report=html``.
This will generate an ``htmlcov`` folder where all the static html files can be found.

To be sure that your modifications didn't break the notebooks, do::

    python -m pip install .[notebooks]
    python -m pytest --nbmake "./notebooks"

Adding a changelog entry
------------------------

Add a short description of your modification in ``CHANGELOG.md`` under the
``[unreleased]`` part.


Building the documentation
--------------------------

To see the documentation locally, you'll need to install the additional python dependencies with::

    pip install .[docs]

and the pandoc software (``sudo apt-get install pandoc`` on ubuntu,
``choco install pandoc`` on windows, ``brew install pandoc`` on mac,
``pacman -S haskell-pandoc`` on arch).

- To build the docs from the repo directory::

    cd docs
    make html
    cd ..

You will find the html index file in the :file:`docs/_build/html` folder.
