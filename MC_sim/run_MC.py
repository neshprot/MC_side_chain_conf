import configparser
import json

from utils import *
from graph import Graph
from MC_code import MonteCarlo

# PARSING CONFIG
config = configparser.ConfigParser()
config.read('config.ini')

# config constants
pdb_file = config['PDB']['File']
value1 = float(config['PDB']['ENERGY'])
psf_file = config['PSF']['File']
inp_file = config['INP']['File']
attempts = float(config['PARAMS']['attempts'])
stop_step = float(config['PARAMS']['stop_step'])
rotating_resid = json.loads(config.get('ROTATING RESID', 'numbers'))

# config files
result_file_name = config['COMPUTING']['ResultFileName']

if __name__ == '__main__':
    const_dict = read_inp(inp_file)
    mol = read_pdb(pdb_file, const_dict)
    bonds, rot_bonds, rot_bonds_CA = amino_acid(mol, rotating_resid)
    read_psf(psf_file, mol)

    graph = Graph(bonds)

    rotations, best_energy = MonteCarlo(mol, graph, rot_bonds, attempts, stop_step, rotating_resid, rot_bonds_CA)

    write_result(result_file_name, rotations, best_energy)
