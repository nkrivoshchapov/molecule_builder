from molbuilder import *

logger = createLogger("MainScript")
logger.info("Starting the script")

# TODO atom indexing from 1

mol = Molecule()
mol.parseMol("pentane.sdf")
bond_length = mol.get_bond_length(0, 1)

frame1 = mol.create_frame(1, 0, 2)
frameB1 = mol.create_frame(3, 4, 2)
frameB2 = mol.create_frame(4, 3, 16)
fragment1 = mol.create_fragment(1, 0)
frame1 = fragment1.find_frame(frame1)

frame2 = mol.create_frame(0, 1, 7)
logger.info("Frame 1 = " + repr(frame1.matrix))
logger.info("Frame 2 = " + repr(frame2.matrix))
fragment2 = mol.create_fragment(0, 1)
frame2 = fragment2.find_frame(frame2)
frame2.shift(bond_length)
fragment2.fix_indicies()
frame2.join(frame1)
fragment2.writeMol("test.mol")


logger.info("Frame 1 = " + repr(frame1.matrix))
logger.info("Frame 2 = " + repr(frame2.matrix))
# mol = Molecule()
# mol.parseMol("pentane.sdf")
# bl1 = mol.get_bond_length(0, 1)
# bl2 = mol.get_bond_length(3, 4)
# frameA1 = mol.create_frame(1, 0, 2)
# frameA2 = mol.create_frame(0, 1, 7)
# frameB1 = mol.create_frame(3, 4, 2)
# frameB2 = mol.create_frame(4, 3, 16)
# middle_frag = mol.create_fragment(1, 0)
# middle_frag = middle_frag.create_fragment(3, 4)
# fragA = mol.create_fragment(0, 1)
# fragB = mol.create_fragment(4, 3)
#
# frame1 = fragA.find_frame(frameA2)
# frame2 = middle_frag.find_frame(frameA1)
# frame1.shift(bl1)
# logger.info("First join")
# frame1.join(frame2)
#
# fragA.writeMol("test.mol")

# frame1 = fragA.find_frame(frameB1)
# frame2 = fragB.find_frame(frameB2)
# frame1.shift(bl2)
# logger.info("Second join")
# frame2.join(frame1)
#
# fragB.writeMol("test.mol")
