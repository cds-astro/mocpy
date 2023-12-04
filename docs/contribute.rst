:ref:`index`

Contribute
==========

For contribution purposes the best thing to do is setting up a python
environment that will act independently from your current installed
interpreter and package manager.

Setting up the environment
--------------------------

Here we highlight two ways to set up an environment.
Using conda (installed from miniconda or anaconda) or
using the pip package
`virtualenv <https://python-guide-pt-br.readthedocs.io/fr/latest/dev/virtualenvs.html>`__ .
Choose the one you prefer, or your own one.

Using conda
~~~~~~~~~~~

- Set up the conda environment and activate it::

    conda create -n mocpy-dev python==3.12.*
    source activate mocpy-dev

- Once you are done with your developments you can
  deactivate your conda environment::

    conda deactivate

Using virtualenv
~~~~~~~~~~~~~~~~

- Go to your home location::

    cd ~

- Set up the virtual environment there::

    virtualenv -p /usr/bin/python3 mocpy-env

``Virtualenv`` will create a directory named ``mocpy-env`` in your home.
This directory contains a new independent python interpreter
(in this case, a python3 one, instanciated from /usr/bin/python3)
along with a new empty pip package manager.

- Activate your virtual environment::

    source ~/mocpy-env/bin/activate

``pip list`` will tell you that there is no package installed and ``python``
opens a session with the ``mocpy-env`` interpreter.

- You can now install all the necessary pip packages
  for developing and testing MOCpy::

    pip install .[dev]

- Once you are done with your developments you can deactivate the virtual env::

    deactivate

Pre-commits setup
-----------------

- You'll need to install the pre-commits in your ``.git/hooks`` to check your code locally::

    pre-commit install

It will run linting and formatting tests at each of your commits.

Now build package
-----------------

- You will need to install `maturin <https://github.com/PyO3/maturin>`__, a tool that builds and publishes crates and rust binaries as python packages::

    pip install maturin

- Move to your root's mocpy location and run maturin::

    maturin develop --release

This step will inform you of any issue in the rust part.

- After a new version of mocpy goes out, if a ``maturin develop --release`` does not actualize your
  ``Cargo.toml`` file, you might need to before executing the ``maturin`` command again::

    rm Cargo.lock && cargo clean

Running the python tests
------------------------

Once your environment is set up and activated you can run the tests

- To run the automated tests and the doctring examples, go to the repo folder and type::

    python -m pytest -v python/mocpy

- To run the tests with coverage report locally::

    python -m pytest -v python/mocpy --cov-report=term --cov=python/mocpy

- When contributing to the notebooks::

    python -m pip install .[notebooks]
    python -m pytest --nbmake -n=auto "./notebooks"

You also can have a html output of the coverage with the flag ``--cov-report=html``.
This will generate an ``htmlcov`` folder where all the static html files can be found.


Building the documentation
--------------------------

To see the documentation locally, you'll need to install the additional python dependencies with

    pip install .[docs]

and the pandoc software (``sudo apt-get install pandoc`` on ubuntu,
``choco install pandoc`` on windows, ``brew install pandoc`` on mac,
``pacman -S haskell-pandoc`` on arch).

- To build the docs from the repo directory::

    cd docs
    make html
    cd ..

You will find the html index file in the :file:`docs/_build/html` folder.
