from molbuilder import *

logger = createLogger("MainScript")
logger.info("Starting the script")

# The order in which `find_frame` and `join` procedures are called matters here

mol = Molecule()
mol.parseMol("pentane.sdf")
bl1 = mol.get_bond_length(0, 1)
bl2 = mol.get_bond_length(3, 4)
frameA1 = mol.create_frame(1, 0, 2)
frameA2 = mol.create_frame(0, 1, 7)
frameB1 = mol.create_frame(3, 4, 2)
frameB2 = mol.create_frame(4, 3, 16)

middle_frag = mol.create_fragment(1, 0)
middle_frag = middle_frag.create_fragment(3, 4)
fragA = mol.create_fragment(0, 1)
fragB = mol.create_fragment(4, 3)

frameA1 = middle_frag.find_frame(frameA1)
frameA2 = fragA.find_frame(frameA2)
frameB2 = fragB.find_frame(frameB2)

frameA2.shift(bl1)
frameA2.rotate(15)
logger.info("First join")
frameA2.join(frameA1)

frameB1 = fragA.find_frame(frameB1) # Can't call sooner because `join` adds necessary atoms in fragA

frameB1.shift(bl2)
frameB1.rotate(-17)
logger.info("Second join")
frameB2.join(frameB1)

fragB.writeMol("test.mol")
