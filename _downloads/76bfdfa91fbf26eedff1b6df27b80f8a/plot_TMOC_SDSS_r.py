from mocpy import TimeMOC

filename = "./../../resources/TMOC/HST_SDSSg/TMoc.fits"
tmoc = TimeMOC.from_fits(filename)
tmoc.plot(title="Time coverage of the SDSS-g survey")
