"""
py23_compat.py

Python 2 / 3 compatibility layer.

"""
try:
    int = long
    range = xrange
except NameError:
    pass

try:
    from sets import Set as set
except ImportError:
    pass



