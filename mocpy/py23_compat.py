"""
py23_compat.py

Python 2 / 3 compatibility layer.

"""
try:
    int = long
    range = xrange
except NameError:
    int = int
    range = range

try:
    from sets import Set as set
except ImportError:
    pass

try:
    from urllib import urlencode
except ImportError:
    from urllib.parse import urlencode

