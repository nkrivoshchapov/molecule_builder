from molbuilder import *

logger = createLogger("MainScript")
logger.info("Starting the script")

# TODO atom indexing from 1

mol1 = Molecule()
mol1.parseMol("pentane.sdf")
mol2 = Molecule(mol1)                       # mol2 is a copy of mol1
bond_length = mol1.get_bond_length(0, 1)    # get bond length
frame1 = mol1.create_frame(1, 0, 2)
fragment1 = mol1.create_fragment(1, 0)
fragment1.writeMol("test1.mol")

frame2 = mol2.create_frame(0, 1, 7)
fragment2 = mol2.create_fragment(0, 1)
fragment2.associated_frames[0].shift(bond_length)
fragment2.associated_frames[0].join(fragment1.associated_frames[0])

fragment1.merge(fragment2)
fragment1.writeMol("test2.mol")
