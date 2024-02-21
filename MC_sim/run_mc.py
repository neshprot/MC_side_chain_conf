"""
read config file and run script
"""
import configparser
from json import loads


from utils import read_inp, read_pdb, amino_acid, read_psf, read_results, write_result, write_pdb, rotate_trp_tors_angle
from graph import Graph
from mc_code import montecarlo

# PARSING CONFIG
config = configparser.ConfigParser()
config.read('/home/mikhail/work/Protein_Model/MC_side_chain_conf/MC_sim/config.ini')


# config constants
pdb_file = config['PDB']['file']
pdb_wild_file = config['PDB']['Wild_file']
value1 = float(config['PDB']['ENERGY'])
psf_file = config['PSF']['File']
inp_file = config['INP']['File']
attempts = float(config['PARAMS']['attempts'])
stop_step = float(config['PARAMS']['stop_step'])
rotating_resid = loads(config.get('ROTATING RESID', 'numbers'))
trp_resid = loads(config.get('ROTATING RESID', 'trp_numbers'))

# config files
result_file_name = config['COMPUTING']['ResultFilename']

if __name__ == '__main__':

    const_dict = read_inp(inp_file)
    mol = read_pdb(pdb_file, const_dict)
    wild_mol = read_pdb(pdb_wild_file, const_dict)
    mol = rotate_trp_tors_angle(wild_mol, mol, trp_resid)

    bonds, rot_bonds = amino_acid(mol, rotating_resid)
        
    read_psf(psf_file, mol)

    graph = Graph(bonds)

    read_results('result1', mol, graph)


    rotations, best_energy = montecarlo(mol, graph, rot_bonds, attempts, stop_step, rotating_resid)

    write_result(result_file_name, rotations)
