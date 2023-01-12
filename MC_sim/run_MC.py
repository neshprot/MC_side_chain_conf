from utils import *
from graph import Graph
from MC_code import MonteCarlo

if __name__ == '__main__':
    mol = read_pdb('6GUX_t.pdb')
    bonds, rot_bonds = amino_acid(mol)

    graph = Graph(bonds)

    start_energy = 100

    print(MonteCarlo(mol, graph, rot_bonds, start_energy, attempts=100))
