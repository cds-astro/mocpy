import nbformat
from nbconvert.preprocessors import CellExecutionError
from nbconvert.preprocessors import ExecutePreprocessor

from glob import glob

notebook_filenames_l = glob("notebooks/*.ipynb")

for notebook_filename in notebook_filenames_l:
    with open(notebook_filename) as f:
        print("Executing notebook : {0}".format(notebook_filename))
        nb = nbformat.read(f, as_version=4)
        ep = ExecutePreprocessor(timeout=1000, kernel_name='mocpy-env')
        try:
            ep.preprocess(nb, {'metadata': {'path': 'notebooks/'}})
        except CellExecutionError as e:
            print("{0} [FAILED]\n{1}".format(notebook_filename, e))
            # exit with error status code
            exit(1)
        print("{0} [PASSED]".format(notebook_filename))

