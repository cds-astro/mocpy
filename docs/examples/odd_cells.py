from mocpy import MOC

moc = MOC.from_json({"1": [i for i in range(12 * 4) if i % 2 == 1]})
moc.display_preview()
