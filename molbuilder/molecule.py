import networkx as nx
import numpy as np
from numpy.linalg import norm
from .logging import createLogger
from copy import deepcopy
from .frame import Frame

logger = createLogger("Molecule")

class Molecule:
    def __init__(self, origin = None):
        if origin is not None: # A simple copy constructor
            self.G = deepcopy(origin.G)
            self.associated_frames = []
            for f in origin.associated_frames:
                newframe = Frame(f)
                newframe.mol = self
                self.associated_frames.append(newframe)
        else: # Create the main structures
            self.G = nx.Graph()
            self.associated_frames = []

    def parseMol(self, filename):
        lines = open(filename, "r").readlines()
        natoms = int(lines[3][0:3])
        nbonds = int(lines[3][3:6])
        logger.debug("NAtoms = %d   NBonds = %d" % (natoms, nbonds))
        for i in range(4, 4 + natoms):
            parts = lines[i].replace("\n", "").split()
            self.G.add_node(i-4)
            self.G.nodes[i - 4]['xyz'] = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
            self.G.nodes[i - 4]['atom_symbol'] = parts[3]
        for i in range(4 + natoms, 4 + natoms + nbonds):
            at1 = int(lines[i][0:3])
            at2 = int(lines[i][3:7])
            self.G.add_edge(at1 - 1, at2 - 1)
            self.G[at1 - 1][at2 - 1]['bondtype'] = int(lines[i][7:10])

    def fix_indicies(self, start = 0):
        self.G = nx.convert_node_labels_to_integers(self.G, first_label = start)

    def get_bond_length(self, at1, at2):
        return norm(self.G.nodes[at1]['xyz'] - self.G.nodes[at2]['xyz'])

    def writeMol(self, filename):
        self.fix_indicies()
        lines = ["", "", ""]
        lines.append("%3d%3d  0  0  0  0  0  0  0  0999 V2000" % (self.G.number_of_nodes(), self.G.number_of_edges()))
        #   -5.5250    1.6470    1.0014 C   0  0  0  0  0  0  0  0  0  0  0  0
        for i in list(self.G.nodes()):
            lines.append("%10.4f%10.4f%10.4f%3s  0  0  0  0  0  0  0  0  0  0  0  0" % (
                self.G.nodes[i]['xyz'][0],
                self.G.nodes[i]['xyz'][1],
                self.G.nodes[i]['xyz'][2],
                self.G.nodes[i]['atom_symbol']))

        for edge in list(self.G.edges()):
            lines.append("%3s%3s%3s  0" % (edge[0] + 1,
                                           edge[1] + 1,
                                           self.G[edge[0]][edge[1]]['bondtype']))
        lines.append("M  END\n")

        file = open(filename, "w")
        file.write("\n".join(lines))
        file.close()

    def create_fragment(self, central_atom, other_atom):
        newmol = Molecule()
        newmol.G = deepcopy(self.G)
        newmol.G.remove_edge(central_atom, other_atom)
        for c in nx.connected_components(newmol.G):
            if other_atom in set(c):
                removepart = set(c)
                break
        for i in removepart:
            newmol.G.remove_node(i)

        for i, frame in enumerate(self.associated_frames):
            if frame.center not in removepart:
                newframe = Frame(frame)
                newframe.mol = newmol
                newmol.associated_frames.append(newframe)
        return newmol

    def create_frame(self, at0, at1, at2):
        newframe = Frame()
        newframe.construct(self, at0, at1, at2)
        return newframe

    def merge(self, other):
        self.fix_indicies()
        other.fix_indicies(start = self.G.number_of_nodes())
        for node in set(other.G.nodes()):
            self.G.add_node(node)
            self.G.nodes[node]['xyz'] = deepcopy(other.G.nodes[node]['xyz'])
            self.G.nodes[node]['atom_symbol'] = deepcopy(other.G.nodes[node]['atom_symbol'])

        for edge in set(other.G.edges()):
            self.G.add_edge(*edge)
            self.G[edge[0]][edge[1]]['bondtype'] = deepcopy(other.G[edge[0]][edge[1]]['bondtype'])
