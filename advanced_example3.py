from molbuilder import *

logger = createLogger("MainScript")
logger.info("Starting the script")

pentane_mol = Molecule()
pentane_mol.parseMol("pentane.sdf")
cc_length = pentane_mol.get_bond_length(0, 1)

# Its possible to use bivalent fragments to create chains. Such fragments cannot be reused, so the list `carbene_frags`
# of identical fragments is created.

main_mol = Molecule()
main_mol.parseMol("methane.sdf")
methyl_frame = main_mol.create_frame(1, 0, 2)
carbene_frame = main_mol.create_frame(1, 3, 2)

methylA_frag = main_mol.create_fragment(1, 0)
methylB_frag = main_mol.create_fragment(1, 0)
N = 100
def get_torsion():
    return 40

carbene_frags = []
for i in range(N):
    carbene_frags.append(methylA_frag.create_fragment(1, 3))

methylA_frame = methylA_frag.find_frame(methyl_frame)
methylB_frame = methylB_frag.find_frame(methyl_frame)

carbene_framesA = []
carbene_framesB = []
for i in range(N):
    carbene_framesA.append(carbene_frags[i].find_frame(methyl_frame))
    carbene_framesB.append(carbene_frags[i].find_frame(carbene_frame))

carbene_framesA[0].shift(cc_length)
carbene_framesA[0].rotate(get_torsion())
carbene_framesA[0].join(methylA_frame)

for i in range(N-1):
    carbene_framesA[i + 1].shift(cc_length)
    carbene_framesA[i + 1].rotate(get_torsion())
    carbene_framesA[i + 1].join(carbene_framesB[i])

carbene_framesB[N - 1].shift(cc_length)
carbene_framesB[N - 1].rotate(get_torsion())
carbene_framesB[N - 1].join(methylB_frame)

carbene_frags[N - 1].writeMol("test.mol")
