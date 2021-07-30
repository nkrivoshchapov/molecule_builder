import numpy as np
from numpy.linalg import norm, inv, det
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
        self.matrix = np.zeros((4, 4))
        self.matrix[:3, :3] = np.array([av, bv, cv]).transpose()
        self.matrix[:3, 3] = deepcopy(mol.G.nodes[at0]['xyz']).transpose()
        self.matrix[3, 3] = 1
        mol.associated_frames.append(self)

        self.mol.fix_indicies()
        for i, frame in enumerate(self.mol.associated_frames):
            logger.info("Old frame : " + repr(frame.matrix))
            logger.info("Old center : " + repr(self.mol.G.nodes[frame.center]['xyz']))

    def shift(self, dist):
        shift_matrix = np.identity(4)
        shift_matrix[0, 3] = dist
        self.matrix = self.matrix @ shift_matrix

    def join(self, other):
        logger.info("Running join")
        myframe = deepcopy(self.matrix)
        myframe[:, 0] = -myframe[:, 0]
        myframe[:, 2] = -myframe[:, 2]
        basis_change =  other.matrix @ inv(myframe)
        logger.warning("Det = %f" % det(basis_change))
        for i in set(self.mol.G.nodes()):
            vec = np.array([self.mol.G.nodes[i]['xyz'][0],
                            self.mol.G.nodes[i]['xyz'][1],
                            self.mol.G.nodes[i]['xyz'][2],
                            1])
            vec = basis_change @ vec
            self.mol.G.nodes[i]['xyz'][:3] = vec[:3]

        self.mol.fix_indicies()
        other.mol.fix_indicies(start=self.mol.G.number_of_nodes())
        for i, frame in enumerate(self.mol.associated_frames):
            logger.info("Old frame : " + repr(frame.matrix))
            logger.info("Old center : " + repr(self.mol.G.nodes[frame.center]['xyz']))

        for node in set(other.mol.G.nodes()):
            self.mol.G.add_node(node)
            self.mol.G.nodes[node]['xyz'] = deepcopy(other.mol.G.nodes[node]['xyz'])
            self.mol.G.nodes[node]['atom_symbol'] = deepcopy(other.mol.G.nodes[node]['atom_symbol'])

        for edge in set(other.mol.G.edges()):
            self.mol.G.add_edge(*edge)
            self.mol.G[edge[0]][edge[1]]['bondtype'] = deepcopy(other.mol.G[edge[0]][edge[1]]['bondtype'])

        for i, frame in enumerate(self.mol.associated_frames):
            self.mol.associated_frames[i].matrix = basis_change @ frame.matrix
            logger.info("New frame : " + repr(self.mol.associated_frames[i].matrix))
            logger.info("New center : " + repr(self.mol.G.nodes[self.mol.associated_frames[i].center]['xyz']))

        for frame in other.mol.associated_frames:
            newframe = Frame(frame)
            newframe.mol = self.mol
            self.mol.associated_frames.append(newframe)

        # Create bond between frame centers
        self.mol.G.add_edge(self.center, other.center)
        self.mol.G[self.center][other.center]['bondtype'] = 1
        self.mol.fix_indicies()
