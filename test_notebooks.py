import nbformat
from nbconvert.preprocessors import CellExecutionError
from nbconvert.preprocessors import ExecutePreprocessor

from glob import glob

class bcolors:
    OKGREEN = '\033[92m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'


notebook_filenames_l = glob("notebooks/*.ipynb")

for notebook_filename in notebook_filenames_l:
    with open(notebook_filename) as f:
        print("Executing notebook : {0}".format(notebook_filename))
        nb = nbformat.read(f, as_version=4)
        ep = ExecutePreprocessor(timeout=1000, kernel_name='python3')
        try:
            ep.preprocess(nb, {'metadata': {'path': 'notebooks/'}})
        except CellExecutionError as e:
            print(("{0} " + bcolors.FAIL + "[FAILED]" + bcolors.ENDC + "\n{1}").format(notebook_filename, e))
            # exit with error status code
            exit(1)
        print(("{0} " + bcolors.OKGREEN + "[PASSED]" + bcolors.ENDC).format(notebook_filename))

