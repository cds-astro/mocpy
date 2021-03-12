*****
MOCPy
*****
|Travis Status| |AppVeyor status| |Notebook Binder| |Doc|

MOCPy is a Python library allowing easy creation and manipulation of MOCs (Multi-Order Coverage maps).   

MOC is an IVOA standard  enabling description of arbitrary sky regions.  
Based on the HEALPix sky tessellation, it maps regions on the sky
into hierarchically grouped predefined cells.

An experimental support for TMOC (temporal MOC) has been added since version 0.4.0.
It allows creation, parsing and comparison of TMOCs.

Space & Time coverages (STMOC) are an extension of MOC to add time information.
It is possible to get a TMOC by querying a STMOC with a MOC and/or get a MOC 
by querying a STMOC with a TMOC.

Please check the mocpy's `documentation <https://cds-astro.github.io/mocpy/>` __
for more details about installing MOCPy and using it.

.. figure:: ./resources/readme.png
   :scale: 50 %
   :align: center
   :alt: map to buried treasure

   *Rendered with MOCpy!*


.. |Travis Status| image:: http://img.shields.io/travis/cds-astro/mocpy.svg?branch=master
    :target: https://travis-ci.org/cds-astro/mocpy

.. |AppVeyor status| image:: https://ci.appveyor.com/api/projects/status/26xwvddah60lhxrx/branch/master?svg=true
    :target: https://ci.appveyor.com/project/bmatthieu3/mocpy/branch/master

.. |Notebook Binder| image:: http://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/cds-astro/mocpy/master

.. |Doc| image:: https://img.shields.io/badge/Documentation-link-green.svg
    :target: https://cds-astro.github.io/mocpy/
