:ref:`index`

Contribute
==========

TODO: describe how to hack on mocpy.

- Set up the conda environment::

    conda env create -q python=3.5 -f environment.yml
    source activate mocpy-dev

- To run the automated tests type::

    python -m pytest -v mocpy

- To build the docs from the repo directory::

    cd docs
    make html
    cd ..

You will find the html index file in the :file:`docs/_build/html` folder.

- To run the tests with coverage report locally go to the repo folder and type::

    pytest -v mocpy --cov-report=term --cov=mocpy mocpy

