from molbuilder import *

logger = createLogger("MainScript")
logger.info("Starting the script")

mol = Molecule()
mol.parseMol("pentane.sdf")
bond_length = mol.get_bond_length(0, 1)

frame1 = mol.create_frame(1, 0, 2)
frame2 = mol.create_frame(0, 1, 7)

fragment1 = mol.create_fragment(1, 0)
fragment2 = mol.create_fragment(0, 1)

frame1 = fragment1.find_frame(frame1)
frame2 = fragment2.find_frame(frame2)

frame2.shift(bond_length)
frame2.join(frame1)
fragment2.writeMol("test.mol")
