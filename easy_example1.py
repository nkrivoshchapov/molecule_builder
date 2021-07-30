from molbuilder import *

logger = createLogger("MainScript")
logger.info("Starting the script")

main_mol = Molecule()
main_mol.parseMol("artemisinin.sdf")
bond_length = main_mol.get_bond_length(17, 19)

main_frame = main_mol.create_frame(17, 19, 16)
main_fragment = main_mol.create_fragment(17, 19)
main_frame = main_fragment.find_frame(main_frame)

substituent_mol = Molecule()
substituent_mol.parseMol("toluene.sdf")

substituent_frame = substituent_mol.create_frame(6, 12, 1)
substituent_fragment = substituent_mol.create_fragment(6, 12)
substituent_frame = substituent_fragment.find_frame(substituent_frame)
substituent_frame.shift(bond_length)
substituent_frame.rotate(60)
main_frame.join(substituent_frame)
main_fragment.writeMol("test.mol")
