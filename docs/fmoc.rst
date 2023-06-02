Working with FMOCs
==================

The cells in this documentation page execute one after the other, just like in a notebook.

.. jupyter-execute::

    from mocpy import FrequencyMOC
    import astropy.units as u
    import matplotlib.pyplot as plt

Frequency MOCs are coverages of the electromagnetic axis.

Basic FMOC creation
-------------------

They can be created from a list of frequencies expressed in Hertz, with `mocpy.FrequencyMOC.from_frequencies`.

.. jupyter-execute::

    fmoc_from_frequencies = FrequencyMOC.from_frequencies(10, [1, 0.1, 0.01] * u.Hz)

Other ways to instantiate a frequency MOC are `mocpy.FrequencyMOC.new_empty`



Plotting FMOCs
--------------

This very simple fmoc can be plotted with `plot_frequencies`

.. jupyter-execute::

    fix, ax = plt.subplots(figsize=(15, 1))
    fmoc_from_frequencies.plot_frequencies(ax, color="lightblue")

or with `plot_wavelengths`