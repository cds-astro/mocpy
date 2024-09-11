from mocpy import MOC

all_sky = MOC.from_str("0/0-11")  # the 12 cells of HEALPix at order 0
all_sky.display_preview()  # the inside of the MOC is represented in red with display_preview
