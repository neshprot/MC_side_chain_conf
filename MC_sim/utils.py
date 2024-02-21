"""
useful scripts for code
"""
import math
from pathlib import Path
import os.path
import numpy as np
from graph import Graph

vdwr = {
    "H": 1.10,
    "O": 1.52,
    "N": 1.55,
    "C": 1.70,
    "S": 1.80
}

resname_3to1 = {
    "GLY": "G",
    "LEU": "L",
    "TYR": "Y",
    "SER": "S",
    "GLU": "E",
    "GLN": "Q",
    "ASP": "D", 
    "ASN": "N",
    "PHE": "F",
    "ALA": "A",
    "LYS": "K",
    "ARG": "R",
    "HIS": "H",
    "CYS": "C",
    "VAL": "V",
    "PRO": "P",
    "TRP": "W",
    "ILE": "I",
    "MET": "M",
    "THR": "T",

    # Not standard
    "ASH": "D",
    "HSP": "H",
    "HSE": "H",
    "GLH": "E",
}

atoma_radius = {
    "H": 0.53,
    "N": 0.75,
    "O": 0.48,
    "C": 0.7,
    "S": 1.27
}


class TorsionAngle:
    def __init__(self, residue_number):
        self.residue_number = residue_number
        self.torsion_angles = {}
        self.rotamer_dict = {}


    def update_torsion_angles(self, bond, new_torsion):
        if bond in self.torsion_angles:
            self.torsion_angles[bond].append(new_torsion)
        else:
            self.torsion_angles[bond] = [new_torsion]
    def update_rotamer(self, rotamers):
        self.rotamer_dict = rotamers


    def print_angels_bond(self, bond):
        return self.torsion_angles[bond]
    

    def print_angles(self):
        return self.torsion_angles
    def rotamers(self):
        return self.rotamer_dict


class Atom:
    """
    describing atom
    """
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


def read_inp(fname):
    """
    reading input file
    :param fname: file name
    :return: LJ constants such epsilon and rmin
    """
    myself = Path(__file__).resolve()
    res = myself.parents[1]
    name_part = ''
    const_dict = {}
    with open(fr"{res}/input_files/{fname}", "r", encoding="utf-8") as file:
        for line in file.readlines():
            newname = line[0:9].strip()
            if newname == 'NONbonded':
                name_part = newname
            if name_part == 'NONbonded' and line[7:15] == '0.000000':
                split_line = line.split()
                name = split_line[0]
                epsilon = float(split_line[2])
                rmin = float(split_line[3])
                const_dict.update({name: (epsilon, rmin)})
    return const_dict


def read_pdb(fname, const_dict):
    '''
    not only read pdb, but also initialize constants from inp file
    '''
    myself = Path(__file__).resolve()
    res = myself.parents[1]
    with open(fr"{res}/input_files/{fname}", "r", encoding="utf-8") as file:
        molecule = {}
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
            atom.epsilon = const_dict.get(atom.name[0], (0.0, 0.0))[0]
            atom.rmin = const_dict.get(atom.name[0], (0.0, 0.0))[1]

            molecule[num] = atom

    return molecule


def read_psf(fname, mol):
    """
    read psf file
    :param fname: file name 
    :param mol: class for protein
    :return: initialisation of atom's params
    """
    myself = Path(__file__).resolve()
    res = myself.parents[1]
    name = ''
    with open(fr"{res}/input_files/{fname}", "r", encoding="utf-8") as file:
        def natom():
            unused_zero = line[69:70].strip()
            if unused_zero == '0':
                num = int(line[0:8])
                print(num)
                mol[num].partialcharge = float(line[35:44])

        def nbond():
            atom_list = line.split()
            first_atom = atom_list[0::2]
            second_atom = atom_list[1::2]
            for i, j in zip(first_atom, second_atom):
                mol[float(i)].bonded.append(float(j))
                mol[float(j)].bonded.append(float(i))

        def ntheta():
            pass

        def section_name(name):
            return {'!NATOM': natom,
                    '!NBOND': nbond,
                    '!NTHETA': ntheta,
                    }.get(name, lambda: 0)()

        for line in file.readlines():
            newname = line[9:15].strip()
            if newname in ['!NATOM', '!NBOND', '!NTHET']:
                name = newname
                continue
            
            section_name(name)


def write_pdb(molecule, fname):
    """
    write pdb file with protein
    :param molecule: protein
    :param fname: file name
    :return: pdb file
    """
    with open(fname, "w", encoding="utf-8") as file:
        for (idx, atom) in molecule.items():
            line = '{a0:<6}{a1:>5}{s}{a2:>4}{a3:>1}{a4:>3}{s}{a5:>1}' \
                   '{a6:>4}{a7:<1}{s:>3}{a8[0]:>8.3f}{a8[1]:>8.3f}{a8[2]:>8.3f}' \
                   '{a9:>6.2f}{a10:>6.2f}{s:>11}{a11:<2}\n'.format(
                a0=atom.datatype, a1=idx, a2=atom.name, a3=atom.altloc,
                a4=atom.resname, a5=atom.chainid, a6=atom.resseq,
                a7=atom.rescode, a8=atom.coordin, a9=atom.occup,
                a10=atom.tempfac, a11=atom.element, s=' '
            )
            file.write(line)


def amino_acid(mol, rotating_resid):
    """graph representation of side bonds"""
    bonds = {}
    rot_bonds = []

    def ala(pos):
        # CA, C, CB, HB1, HB2, HB3
        bonds.update({pos+6: [pos], pos: [pos+2], pos+2: [pos+3, pos+4, pos+5],
                      pos+3: [pos+3], pos+4: [pos+4],
                      pos+5: [pos+5]})
        nrot_bonds = [(pos, pos+2)]
        return nrot_bonds

    def ser(pos):
        # CA, C, CB, HB1, HB2, OG, HG1
        bonds.update({pos + 7: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5],
                      pos + 3: [pos + 3], pos + 4: [pos + 4],
                      pos + 5: [pos + 6], pos + 6: [pos + 6]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 5)]
        return nrot_bonds

    def thr(pos):
        # CA, CB, HB, OG1, HG1, CG2, HG21, HG22, HG23
        bonds.update(
            {pos+10: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 6],
             pos + 3: [pos + 3], pos + 4: [pos + 5],
             pos + 5: [pos + 5], pos + 6: [pos + 7, pos + 8, pos + 9],
             pos + 7: [pos + 7], pos+8: [pos+8], pos+9: [pos+9]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 4), (pos + 2, pos + 6)]
        return nrot_bonds

    def cys(pos):
        # CA, CB, HB1, HB2, SG, HG1
        bonds.update(
            {pos + 7: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5],
             pos + 3: [pos + 3], pos + 4: [pos + 4],
             pos + 5: [pos + 6], pos + 6: [pos + 6]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 5)]
        return nrot_bonds

    def val(pos):
        # CA, CB, HB, CG1, HG11, HG12, HG13, CG2, HG21, HG22, HG23
        bonds.update(
            {pos+12: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 8],
             pos + 3: [pos + 3], pos + 4: [pos + 5, pos+6, pos+7],
             pos + 5: [pos + 5], pos + 6: [pos+6],
             pos + 7: [pos+7], pos+8: [pos+9, pos+10, pos+11],
             pos+9: [pos+9], pos+10: [pos+10],
             pos+11: [pos+11]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 4), (pos + 2, pos + 8)]
        return nrot_bonds

    def leu(pos):
        # CA, CB, HB1, HB2, CG, HG, CD1, HD11, HD12, HD13, CD2, HD21, HD22, HD23
        bonds.update(
            {pos+15: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5],
             pos + 3: [pos + 3], pos + 4: [pos + 4],
             pos + 5: [pos + 6, pos + 7, pos + 11], pos + 6: [pos + 6],
             pos + 7: [pos + 8, pos + 9, pos + 10], pos + 8: [pos + 8],
             pos + 9: [pos + 9], pos + 10: [pos + 10],
             pos + 11: [pos + 12, pos+13, pos+14],
             pos+12: [pos+12], pos+13: [pos+13], pos+14: [pos+14]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 5), (pos + 5, pos + 7), (pos + 5, pos + 11)]
        return nrot_bonds

    def ile(pos):
        # CA, CB, HB, CG2, HG21, HG22, HG23, CG1, HG11, HG12, CD, HD1, HD2, HD3
        bonds.update(
            {pos + 15: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 8],
             pos + 3: [pos + 3], pos + 4: [pos + 5, pos + 6, pos + 7],
             pos + 5: [pos + 5], pos + 6: [pos + 6],
             pos + 7: [pos + 7], pos + 8: [pos + 9, pos + 10, pos + 11],
             pos + 9: [pos + 9], pos + 10: [pos + 10],
             pos + 11: [pos + 12, pos + 13, pos + 14], pos + 12: [pos + 12],
             pos + 13: [pos + 13], pos + 14: [pos + 14]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 8)]
        return nrot_bonds

    def met(pos):
        # CA, CB, HB1, HB2, CG, HG1, HG2, SD, CE, HE1, HE2, HE3
        bonds.update(
            {pos + 13: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5],
             pos + 3: [pos + 3], pos + 4: [pos + 4],
             pos + 5: [pos + 6, pos + 7, pos + 8], pos + 6: [pos + 6],
             pos + 7: [pos + 7], pos + 8: [pos + 9],
             pos + 9: [pos + 10, pos + 11, pos+12], pos + 10: [pos + 10],
             pos + 11: [pos + 11], pos + 12: [pos+12]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 5), (pos + 5, pos + 8), (pos + 8, pos + 9)]
        return nrot_bonds

    def phe(pos):
        # CA, CB, HB1, HB2, CG, CD1, HD1, CE1, HE1, CZ, HZ, CD2, HD2, CE2, HE2
        bonds.update(
            {pos+16: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5], pos + 3: [pos + 3],
             pos + 4: [pos + 4], pos + 5: [pos + 6, pos + 8],
             pos + 6: [pos + 7, pos + 12], pos + 7: [pos + 7],
             pos + 8: [pos + 9, pos+14], pos + 9: [pos + 9],
             pos + 10: [pos + 11], pos + 11: [pos + 11],
             pos+12: [pos+10, pos+13], pos+13: [pos+13],
             pos+14: [pos+10, pos+15], pos+15: [pos + 15]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 5)]
        return nrot_bonds

    def tyr(pos):
        # CA, CB, HB1, HB2, CG, CD1, HD1, CE1, HE1, CZ, OH, HH, CD2, HD2, CE2, HE2
        bonds.update(
            {pos+17: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5],
             pos + 3: [pos + 3], pos + 4: [pos + 4],
             pos + 5: [pos + 6, pos + 8], pos + 6: [pos + 7, pos + 13],
             pos + 7: [pos + 7], pos + 8: [pos + 9, pos+15],
             pos + 9: [pos + 10], pos + 10: [pos + 11],
             pos + 11: [pos + 12], pos+12: [pos+12],
             pos+13: [pos+10, pos+14], pos+14: [pos+14],
             pos+15: [pos + 10, pos + 16], pos+16: [pos+16]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 5), (pos + 10, pos + 11)]
        return nrot_bonds

    def trp(pos):
        # CA, CB, HB1, HB2, CG, CD1, HD1, NE1, HE1,
        # CE2, CD2, CE3, HE3, CZ3, HZ3, CZ2, HZ2, CH2, HH2
        bonds.update(
            {pos + 20: [pos],pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5],
             pos + 3: [pos + 3], pos + 4: [pos + 4],
             pos + 5: [pos + 6], pos + 6: [pos + 7, pos + 18],
             pos + 7: [pos + 8, pos + 9], pos + 8: [pos + 8],
             pos + 9: [pos + 10, pos + 11], pos + 10: [pos + 10],
             pos + 11: [pos + 12, pos + 18], pos+12: [pos+13, pos+14],
             pos+13: [pos+13], pos+14: [pos+15, pos+16],
             pos+15: [pos + 15], pos+16: [pos+17, pos+18],
             pos+17: [pos+17], pos+18: [pos+19],
             pos+19: [pos+19]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 5)]
        return nrot_bonds

    def asp(pos):
        # CA, CB, HB1, HB2, CG, OD1, OD2
        bonds.update(
            {pos + 8: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5],
             pos + 3: [pos + 3], pos + 4: [pos + 4],
             pos + 5: [pos + 6, pos + 7], pos + 6: [pos + 6],
             pos + 7: [pos + 7]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 5)]
        return nrot_bonds

    def glu(pos):
        # CA, CB, HB1, HB2, CG, HG1, HG2, CD, OE1, OE2
        bonds.update(
            {pos + 11: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5],
             pos + 3: [pos + 3], pos + 4: [pos + 4],
             pos + 5: [pos + 6, pos + 7, pos + 8], pos + 6: [pos + 6],
             pos + 7: [pos + 7], pos + 8: [pos + 9, pos + 10],
             pos + 9: [pos + 9], pos + 10: [pos + 10]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 5), (pos + 5, pos + 8)]
        return nrot_bonds

    def asn(pos):
        # CA, CB, HB1, HB2, CG, OD1, ND2, HD21, HD22
        bonds.update(
            {pos + 10: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5],
             pos + 3: [pos + 3], pos + 4: [pos + 4],
             pos + 5: [pos + 6, pos + 7], pos + 6: [pos + 6],
             pos + 7: [pos + 8, pos + 9], pos + 8: [pos + 8],
             pos + 9: [pos + 9]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 5), (pos + 5, pos + 7)]
        return nrot_bonds

    def lys(pos):
        # CA, CB, HB1, HB2, CG, HG1, HG2, CD, HD1, HD2, CE, HE1, HE2, NZ, NZ1, NZ2, NZ3
        bonds.update(
            {pos + 18: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5],
             pos + 3: [pos + 3], pos + 4: [pos + 4],
             pos + 5: [pos + 6, pos + 7, pos + 8], pos + 6: [pos + 6],
             pos + 7: [pos + 7], pos + 8: [pos + 9, pos + 10, pos + 11],
             pos + 9: [pos + 9], pos + 10: [pos + 10],
             pos + 11: [pos + 12, pos + 13, pos + 14], pos+12: [pos+12],
             pos + 13: [pos + 13], pos + 14: [pos+15, pos + 16, pos + 17],
             pos + 15: [pos + 15], pos + 16: [pos + 16],
             pos + 17: [pos + 17]})
        nrot_bonds = [(pos, pos + 2),   (pos + 2, pos + 5), (pos + 5, pos + 8),
                      (pos + 8, pos + 11), (pos + 11, pos + 14)]
        return nrot_bonds

    def arg(pos):
        # CA, CB, HB1, HB2, CG, HG1, HG2, CD, HD1, HD2, NE, HE, CZ, NH1, HH11, HH12, NH2, HH21, HH22
        bonds.update(
            {pos + 20: [pos], pos: [pos + 2], pos + 2: [pos + 3, pos + 4, pos + 5], pos + 3: [pos + 3],
             pos + 4: [pos + 5], pos + 5: [pos + 6, pos + 7, pos + 8],
             pos + 6: [pos + 6], pos + 7: [pos + 7],
             pos + 8: [pos + 9, pos + 10, pos + 11], pos + 9: [pos + 9],
             pos + 10: [pos + 10], pos + 11: [pos + 12, pos + 13],
             pos + 12: [pos + 12], pos + 13: [pos + 14, pos + 17],
             pos + 14: [pos + 15, pos + 16], pos + 15: [pos + 15],
             pos + 16: [pos + 16], pos + 17: [pos + 18, pos + 19],
             pos + 18: [pos + 18], pos + 19: [pos + 19]})
        nrot_bonds = [(pos, pos + 2), (pos + 2, pos + 5), (pos + 5, pos + 8),
                      (pos + 8, pos + 11), (pos + 11, pos + 13),
                     (pos + 13, pos + 14), (pos + 13, pos + 17)]
        return nrot_bonds

    def side_chains(resname, num):
        return {'ALA': lambda: ala(num),
                'SER': lambda: ser(num),
                'THR': lambda: thr(num),
                'CYS': lambda: cys(num),
                'VAL': lambda: val(num),
                'LEU': lambda: leu(num),
                'ILE': lambda: ile(num),
                'MET': lambda: met(num),
                'PHE': lambda: phe(num),
                'TYR': lambda: tyr(num),
                'TRP': lambda: trp(num),
                'ASP': lambda: asp(num),
                'GLU': lambda: glu(num),
                'ASN': lambda: asn(num),
                'LYS': lambda: lys(num),
                'ARG': lambda: arg(num)}.get(resname, lambda: [])()
    rot_bonds = {}
    for i, atom in enumerate(mol.values()):
        if atom.name == 'CA' and atom.resseq in rotating_resid:
            nrot_bonds = side_chains(atom.resname, i+1)
            rot_bonds[atom.resseq] = nrot_bonds
    return bonds, rot_bonds


def rotation(origin_point1, origin_point2, point, angle):
    """
    rotation around an arbitrary axis

    :param origin_point1: the first point of the vector around which the rotation occurs
    :param origin_point2: the second point of the vector around which the rotation occurs
    :param point: the point that rotates around the vector
    :param angle: angle of rotation of the point around the vector
    :return : new point coordinates
    """
    vector = origin_point2 - origin_point1
    d_vec = math.sqrt(sum(x**2 for x in vector))
    norm_vector = vector/d_vec
    dist = math.sqrt(norm_vector[1]**2+norm_vector[2]**2)

    def count_matrix(origin_point, norm_vector, angle):
        trans_matrix = np.array([[1.0, 0.0, 0.0, 0.0],
                                [0.0, 1.0, 0.0, 0.0],
                               [0.0, 0.0, 1.0, 0.0],
                               [-origin_point[0], -origin_point[1], -origin_point[2], 1.0]])
        rot_x_matrix = np.array([[1.0, 0.0, 0.0, 0.0],
                                   [0.0, norm_vector[2]/dist, norm_vector[1]/dist, 0.0],
                                   [0.0, -norm_vector[1]/dist, norm_vector[2]/dist, 0.0],
                                   [0.0, 0.0, 0.0, 1.0]])
        rot_y_matrix = np.array([[dist, 0.0, norm_vector[0], 0.0],
                                   [0.0, 1.0, 0.0, 0.0],
                                   [-norm_vector[0], 0.0, dist, 0.0],
                                   [0.0, 0.0, 0.0, 1.0]])
        alpha = np.radians(angle)
        rot_a_matrix = np.array([[np.cos(alpha), np.sin(alpha), 0.0, 0.0],
                                   [-np.sin(alpha), np.cos(alpha), 0.0, 0.0],
                                   [0.0, 0.0, 1.0, 0.0],
                                   [0.0, 0.0, 0.0, 1.0]])
        return trans_matrix, rot_x_matrix, rot_y_matrix, rot_a_matrix

    trans_matrix, rot_x_matrix, rot_y_matrix, rot_a_matrix = \
        count_matrix(origin_point1, norm_vector, angle)
    matrix = np.dot(np.dot(trans_matrix, rot_x_matrix), rot_y_matrix)

    change_point = np.append(point, 1)
    inv_matrix = np.copy(matrix)
    inv_matrix = np.transpose(inv_matrix)
    inv_matrix[:, 3] = [0.0, 0.0, 0.0, 1.0]
    inv_matrix[3] = list(origin_point1) + [1]
    final_matrix = np.dot(np.dot(matrix, rot_a_matrix), inv_matrix)
    new_point = np.dot(change_point, final_matrix)
    return np.delete(new_point, 3)


def rotate(mol, coords_for_rot, bond_num1, bond_num2, angle):
    """
    rotation of all points of the side chain connected to the rotation axis

    :param mol: protein
    :param coords_for_rot: coordinates of the points that will be rotated
    :param bond_num1: the number of the first atom in the vector around which the rotation occurs
    :param bond_num2: the number of the second atom in the vector around which the rotation occurs
    :param angle: point rotation angle
    :return: updates the coordinates of atoms in a protein
    """
    for i in coords_for_rot:
        new_coord = rotation(mol[bond_num1].coordin, mol[bond_num2].coordin, mol[i].coordin, angle)
        mol[i].coordin = new_coord


def write_result(fname, rotations):
    '''
    :param fname: name of the file
    :param rotations: list of all rotations with nums of atoms and angle
    :param best_energy: energy of the fragment of side chain
    :return: create file with good rotations
    '''
    with open(fname, "w", encoding="utf-8") as file:
        for i in rotations:
            file.write(f"Rotation around an axis: {i[0]} {i[1]} by an angle: {i[2]}\n")


def read_results(fname, mol, graph):
    """
    function read previous results and rotate protein

    :param fname: dile with previous results
    :param mol: protein
    :param graph: graph class
    :return: rotate protein
    """
    if os.path.isfile(fname):
        with open(fname, "r", encoding="utf-8") as file:
            for line in file.readlines():
                lst = line.split()
                rotate(mol, graph.bfs(float(lst[5])), float(lst[4]), float(lst[5]), float(lst[9]))


def post_proc(ini_mol, start_mol, end_mol, rotating_resid):
    """
    assessment of model fit, need upgrade
    :param ini_mol:
    :param start_mol:
    :param end_mol:
    :param rotating_resid:
    :return:
    """
    value = 0
    counter = 0
    for atom in start_mol:
        if start_mol[atom].resseq not in rotating_resid:
            continue
        ini_dist = sum(ini_mol[atom].coordin[i]**2 - start_mol[atom].coordin[i]**2
                       for i in range(3))
        fin_dist = sum(end_mol[atom].coordin[i] ** 2 - start_mol[atom].coordin[i] ** 2
                       for i in range(3))
        value += 1/(1+math.exp(-abs(ini_dist-fin_dist)))-1/2
        counter += 1
    result = 1 - value/counter
    return result


def get_torsion_angle(point1, point2, point3, point4):
    rad2deg = 180.0 / math.pi
    vector1 = np.cross(np.array(point2) - np.array(point1), np.array(point2) - np.array(point3))
    vector2 = np.cross(np.array(point3) - np.array(point2), np.array(point3) - np.array(point4))  
    norm_vector1 = vector1 / np.linalg.norm(vector1)
    norm_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(norm_vector1, norm_vector2)
    sign = 2 * float(np.dot(norm_vector1, (np.array(point3) - np.array(point4)) / np.linalg.norm(np.array(point3) - np.array(point4))) < 0) - 1
    torsion = -np.arccos(dot_product) * rad2deg * sign
    torsion = (torsion + 360) % 360
    return torsion


def get_torsions(mol, bonds):
    torsions = []
    for position, atoms in bonds.items():
        n_position_bonds = len(bonds[position])
        for a in range(n_position_bonds):
            k = bonds[position][a]
            n_kbonds = len(bonds[k])
            for b in range(n_kbonds):
                i = bonds[k][b]
                n_ibonds = len(bonds[i])
                for c in range(n_ibonds):
                    l = bonds[i][c]
                    if (l == position or l == i):
                        continue
                    t1234 = get_torsion_angle(mol[position].coordin, mol[k].coordin, mol[i].coordin, mol[l].coordin)
                    if mol[k].name[0] == 'C' and mol[i].name[0] == 'C':
                        torsions.append([position, k, i, l, t1234])
    return torsions


def print_torsions(mol, torsions):
    n_torsions = len(torsions)
    print('%i torsion(s) found (degrees) %i in  resseq' % (n_torsions, mol[torsions[0][0]].resseq))
    for q in range(n_torsions):
        n1, n2, n3, n4 = torsions[q][0:4]
        t1234 = torsions[q][4]
        nstr = '%i-%i-%i-%i' % (n1, n2, n3, n4)
        tstr = '(%s-%s-%s-%s) ' % (mol[n1].name, mol[n2].name, mol[n3].name, mol[n4].name)
        print(' %-15s  %-13s  %8.3f\n' % (nstr, tstr, t1234), end='')
    print('\n', end='')


def rotate_trp_tors_angle(mol_1, mol_2, rotating_resid):
    for trp_number in rotating_resid:
        trp_dif = []
        bonds_1, rot_bonds_1 = amino_acid(mol_1, [trp_number])
        bonds_2, rot_bonds_2 = amino_acid(mol_2, [trp_number])
        graph_2 = Graph(bonds_2)
        torsions_1 = [get_torsions(mol_1, bonds_1)[0], get_torsions(mol_1, bonds_1)[1]]
        torsions_2 = [get_torsions(mol_2, bonds_2)[0], get_torsions(mol_2, bonds_2)[1]]
        trp_dif += [torsions_2[0][4] - torsions_1[0][4], torsions_2[1][4] - torsions_1[1][4]]
        rot_1 = trp_dif[0]
        rot_2 = trp_dif[1]
        rotate(mol_2, graph_2.bfs(torsions_2[0][2]), torsions_2[0][1], torsions_2[0][2], (-1) * np.sign(rot_1) * abs(rot_1))
        rotate(mol_2, graph_2.bfs(torsions_2[1][2]), torsions_2[1][1], torsions_2[1][2], (-1) * np.sign(rot_2) * abs(rot_2))
        #print_torsions(mol_1, get_torsions(mol_1, bonds_1))
        #print_torsions(mol_2, [get_torsions(mol_2, bonds_2)[0], get_torsions(mol_2, bonds_2)[1]])
    return mol_2
