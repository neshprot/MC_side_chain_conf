import configparser

from utils import *
from graph import Graph
from MC_code import MonteCarlo


# PARSING CONFIG
config = configparser.ConfigParser()
config.read('config.ini')

# config constants
pdb_file = config['PDB']['File']
value1 = float(config['PDB']['ENERGY'])
attempts = float(config['PARAMS']['attempts'])
stop_step = float(config['PARAMS']['stop_step'])

# config files
result_file_name = config['COMPUTING']['ResultFileName']

if __name__ == '__main__':
    mol = read_pdb('6GUX_t.pdb')
    bonds, rot_bonds = amino_acid(mol)

    graph = Graph(bonds)

    start_energy = 100

    rotations, best_energy = MonteCarlo(mol, graph, rot_bonds, start_energy, attempts, stop_step)

    write_result(result_file_name, rotations, best_energy)
