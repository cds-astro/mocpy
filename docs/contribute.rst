:ref:`index`

Contribute
==========

TODO: describe how to hack on mocpy.

- Set up the conda environment::

    conda env create -q python=3.5 -f environment.yml
    source activate mocpy-dev

- To run the automated tests go to the repo folder and type::

    python -m pytest -v mocpy

- To run the tests with coverage report locally::

    pytest -v mocpy --cov-report=term --cov=mocpy mocpy

- If you want to do some profiling on a specific test use the following command::

    pytest mocpy -v -k <your_test_name>  --profile-svg

This will generate for you a nice SVG graph telling you how much time functions in your test take
- To build the docs from the repo directory::

    cd docs
    make html
    cd ..

You will find the html index file in the :file:`docs/_build/html` folder.

