"""
py23_compat.py

Python 2 / 3 compatibility layer.
"""
import sys

__all__ = [
    'int', 'range', 'urlencode', 'BytesIO',
]

PY2 = (sys.version_info.major == 2)

if PY2:
    int = long
    range = xrange
else:
    int = int
    range = range

if PY2:
    from urllib import urlencode
else:
    from urllib.parse import urlencode

# https://pythonhosted.org/six/#six.BytesIO
if PY2:
    import StringIO

    BytesIO = StringIO.StringIO
else:
    import io

    BytesIO = io.BytesIO
