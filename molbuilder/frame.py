import numpy as np
from numpy.linalg import norm, inv
from .logging import createLogger
from copy import deepcopy

logger = createLogger("Frame")

class Frame:
    def __init__(self, origin = None):
        if origin is not None:
            self.center = origin.center
            self.matrix = deepcopy(origin.matrix)

    def construct(self, mol, at0, at1, at2):
        self.mol = mol
        self.center = at0
        av = mol.G.nodes[at1]['xyz'] - mol.G.nodes[at0]['xyz']
        av /= norm(av)
        bv = mol.G.nodes[at2]['xyz'] - mol.G.nodes[at0]['xyz']
        cv = np.cross(av, bv)
        cv /= norm(cv)
        bv = np.cross(cv, av)
        # logger.info("av = " + repr(av))
        # logger.info("center = " + repr(mol.G.nodes[at0]['xyz']))
        self.matrix = np.zeros((4, 4))
        self.matrix[:3, :3] = np.array([av, bv, cv]).transpose()
        self.matrix[:3, 3] = deepcopy(mol.G.nodes[at0]['xyz']).transpose()
        self.matrix[3, 3] = 1
        # logger.info(repr(self.matrix))
        mol.associated_frames.append(self)

    def shift(self, dist):
        shift_matrix = np.identity(4)
        shift_matrix[0, 3] = dist
        self.matrix = self.matrix @ shift_matrix

    def join(self, other):
        myframe = deepcopy(self.matrix)
        myframe[:, 0] = -myframe[:, 0]
        basis_change = other.matrix @ inv(myframe)
        for i in set(self.mol.G.nodes()):
            vec = np.array([self.mol.G.nodes[i]['xyz'][0],
                            self.mol.G.nodes[i]['xyz'][1],
                            self.mol.G.nodes[i]['xyz'][2],
                            1])
            vec = basis_change @ vec
            self.mol.G.nodes[i]['xyz'][:3] = vec[:3]
