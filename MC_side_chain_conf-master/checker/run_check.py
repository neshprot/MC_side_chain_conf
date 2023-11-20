import numpy as np
import math
from pathlib import Path
import os.path


class Atom:
    def __init__(self):
        self.datatype = 'ATOM'  # "ATOM"/"HETATM"
        self.name = ''  # Atom name
        self.altloc = ''  # Alternate location indicator.
        self.resname = ''  # Residue name
        self.chainid = 'U'  # Chain identifier
        self.resseq = 0  # Residue sequence number
        self.rescode = ''  # Code for insertions of residues
        self.coordin = np.array([0, 0, 0])  # (X,Y,Z) orthogonal
        self.occup = 0.0  # occupancy
        self.tempfac = 0.0  # Temperature factor
        self.element = 'Xx'  # element symbol
        self.charge = 0.0  # Atom charge
        self.partialcharge = 0.0 # Atom partial charge
        self.bonded = []    # atoms that is bonded to the current one
        self.epsilon = 0.0  # epsilon LJ
        self.rmin = 0.0     # rmin/2 LJ

def read_pdb(fname, input_dir):
    '''
    not only read pdb, but also initialize constants from inp file
    '''
    myself = Path(__file__).resolve()
    res = myself.parents[1]

    with open(f"{res}\{input_dir}\{fname}", "r") as file:
        molecule = dict()
        for line in file.readlines():
            data_type = line[0:6].strip()
            if data_type not in ['ATOM', 'HETATM']:
                continue
            atom = Atom()
            atom.datatype = line[0:6].strip()
            num = int(line[6:11])
            atom.name = line[12:16].strip()
            atom.altloc = line[16].strip()
            atom.resname = line[17:20].strip()
            atom.chainid = line[21].strip()
            atom.resseq = int(line[22:26])
            atom.rescode = line[26].strip()
            atom.coordin = np.array(list(map(float, [line[30:38], line[38:46], line[46:54]])))
            atom.occup = 0.0  # float(line[54:60])
            atom.tempfac = 0.0  # float(line[60:66])
            atom.element = atom.name[0]  # line[76:78].strip()

            molecule[num] = atom

    return molecule

def post_proc(ini_mol, start_mol, end_mol, rotating_resid):
    value = 0
    counter = 0
    for atom in start_mol:
        if start_mol[atom].resseq not in rotating_resid:
            continue
        else:
            ini_dist = sum(ini_mol[atom].coordin[i]**2 - start_mol[atom].coordin[i]**2 for i in range(3))
            fin_dist = sum(end_mol[atom].coordin[i] ** 2 - start_mol[atom].coordin[i] ** 2 for i in range(3))
            value += 1/(1+math.exp(-abs(ini_dist-fin_dist))) - 1/2
            counter += 1
    result = 1 - value/counter
    return result

if __name__ == '__main__':
    rotating_resid = [15, 219, 88, 215, 89, 85, 207, 208, 197, 188, 189, 86, 60, 93, 211, 212, 141, 137, 192, 118, 196, 92]
    ini_mol = read_pdb('6GUX_t.pdb', 'input_files')
    start_mol = read_pdb('d95n_0_autopsf_tmpfile_v1.pdb', 'input_files')
    end_mol = read_pdb('out_1000_2000_0.00001_80_20.pdb', 'results')

    #0.5032092470676031 - если не изменять начальную

    result = post_proc(ini_mol, start_mol, end_mol, rotating_resid)
    print(result)
