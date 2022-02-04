:ref:`index`

Contribute
==========

For contribution purposes the best thing to do is setting up a python
environment that will act independently from your current installed
interpreter and package manager.

Setting up the environment
--------------------------

There is two ways to set up an environment.
Using conda (installed from miniconda or anaconda) or
using the pip package called
`virtualenv <https://python-guide-pt-br.readthedocs.io/fr/latest/dev/virtualenvs.html>`__ .
Choose the one you prefer.

Using conda
~~~~~~~~~~~

- Set up the conda environment and activate it::

    conda env create -q python=3.7 -f <path_to_mocpy>/requirements/environment.yml -n mocpy-env
    source activate mocpy-env

- Once you are done with your developments you can
  deactivate your conda environment::

    source deactivate

Using virtualenv
~~~~~~~~~~~~~~~~

- Go to your home location::

    cd ~

- Set up the virtual environment there::

    virtualenv -p /usr/bin/python3 mocpy-env

`Virtualenv` will create a directory named `mocpy-env` in your home.
This directory contains a new independent python interpreter
(in this case, a python3 one, instanciated from /usr/bin/python3)
along with a new empty pip package manager.

- Activate your virtual environment::

    source ~/mocpy-env/bin/activate

`pip list` will tell you that there is no package installed and `python`
opens a session with the `mocpy-env` interpreter.

- You will need to install `maturin <https://github.com/PyO3/maturin>`__, a tool that builds and publishes crates and rust binaries as python packages::

    pip install maturin

- Run maturin::

    maturin develop --release

- You can now install all the necessary pip packages
  for developping and testing MOCpy::

    pip install -r <path_to_mocpy_repo>/requirements/contributing.txt

- Once you are done with your developments you can deactivate the virtual env::

    deactivate


Running the tests
-----------------

Once your environment is set up and activated you can run the tests

- To run the automated tests go to the repo folder and type::

    python -m pytest -v mocpy

- To run the tests with coverage report locally::

    python -m pytest -v mocpy --cov-report=term --cov=mocpy mocpy


Building the documentation
--------------------------

- To build the docs from the repo directory::

    cd docs
    make html
    cd ..

You will find the html index file in the :file:`docs/_build/html` folder.
