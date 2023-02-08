import numpy as np
import math
from pathlib import Path

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


class Atom:
    def __init__(self):
        self.DataType = 'ATOM'  # "ATOM"/"HETATM"
        self.Name = ''  # Atom name
        self.AltLoc = ''  # Alternate location indicator.
        self.ResName = ''  # Residue name
        self.ChainId = 'U'  # Chain identifier
        self.ResSeq = 0  # Residue sequence number
        self.ResCode = ''  # Code for insertions of residues
        self.Coordin = np.array([0, 0, 0])  # (X,Y,Z) orthogonal
        self.Occup = 0.0  # Occupancy
        self.TempFac = 0.0  # Temperature factor
        self.Element = 'Xx'  # Element symbol
        self.Charge = 0.0  # Atom charge
        self.PartialCharge = 0.0 # Atom partial charge
        self.Bonded = []    # atoms that is bonded to the current one
        self.Epsilon = 0.0  # epsilon LJ
        self.Rmin = 0.0     # Rmin/2 LJ


def read_inp(fname):
    myself = Path(__file__).resolve()
    res = myself.parents[1]
    name = ''
    const_dict = dict()
    with open(f"{res}\input_files\{fname}", "r") as file:
        for line in file.readlines():
            newname = line[0:9].strip()
            if newname == 'NONBONDED':
                name = newname
            if name == 'NONBONDED' and line[7:15] == '0.000000':
                split_line = line.split()
                Name = split_line[0]
                Epsilon = float(split_line[2])
                Rmin = float(split_line[3])
                const_dict.update({Name: (Epsilon, Rmin)})
    return const_dict


def read_pdb(fname, const_dict):
    '''
    not only read pdb, but also initialize constants from inp file
    '''
    myself = Path(__file__).resolve()
    res = myself.parents[1]
    with open(f"{res}\input_files\{fname}", "r") as file:
        molecule = dict()
        for line in file.readlines():
            data_type = line[0:6].strip()
            if data_type not in ['ATOM', 'HETATM']:
                continue
            atom = Atom()
            atom.DataType = line[0:6].strip()
            num = int(line[6:11])
            atom.Name = line[12:16].strip()
            atom.AltLoc = line[16].strip()
            atom.ResName = line[17:20].strip()
            atom.ChainId = line[21].strip()
            atom.ResSeq = int(line[22:26])
            atom.ResCode = line[26].strip()
            atom.Coordin = np.array(list(map(float, [line[30:38], line[38:46], line[46:54]])))
            atom.Occup = 0.0  # float(line[54:60])
            atom.Tempfac = 0.0  # float(line[60:66])
            atom.Element = atom.Name[0]  # line[76:78].strip()
            atom.Epsilon = const_dict.get(atom.Name, (0.0, 0.0))[0]
            atom.Rmin = const_dict.get(atom.Name, (0.0, 0.0))[1]

            molecule[num] = atom

    return molecule


def read_psf(fname, mol):
    myself = Path(__file__).resolve()
    res = myself.parents[1]
    name = ''
    with open(f"{res}\input_files\{fname}", "r") as file:
        def NATOM():
            unused_zero = line[69:70].strip()
            if unused_zero == '0':
                num = int(line[0:8])
                mol[num].PartialCharge = float(line[35:44])

        def NBOND():
            a = line.split()
            a1 = a[0::2]
            a2 = a[1::2]
            for i, j in zip(a1, a2):
                mol[float(i)].Bonded.append(float(j))
                mol[float(j)].Bonded.append(float(i))

        def NTHETA():
            pass

        def section_name(name):
            return {'!NATOM': lambda: NATOM(),
                    '!NBOND': lambda: NBOND(),
                    '!NTHETA': lambda: NTHETA(),
                    }.get(name, lambda: 0)()

        for line in file.readlines():
            newname = line[9:15].strip()
            if newname in ['!NATOM', '!NBOND', '!NTHET']:
                name = newname
                continue
            section_name(name)


def write_pdb(molecule, fname):
    with open(fname, "w") as f:
        for (idx, atom) in molecule.items():
            line = '{a0:<6}{a1:>5}{s}{a2:>4}{a3:>1}{a4:>3}{s}{a5:>1}' \
                   '{a6:>4}{a7:<1}{s:>3}{a8[0]:>8.3f}{a8[1]:>8.3f}{a8[2]:>8.3f}' \
                   '{a9:>6.2f}{a10:>6.2f}{s:>11}{a11:<2}\n'.format(
                a0=atom.DataType, a1=idx, a2=atom.Name, a3=atom.AltLoc,
                a4=atom.ResName, a5=atom.ChainId, a6=atom.ResSeq,
                a7=atom.ResCode, a8=atom.Coordin, a9=atom.Occup,
                a10=atom.TempFac, a11=atom.Element, s=' '
            )
            f.write(line)
            

def amino_acid(mol, rotating_resid):
    """graph representation of side bonds"""
    bonds = dict()
    rot_bonds = []

    def GLY(x):
        return []

    def ALA(x):
        # CA, CB, HB1, HB2, HB3
        bonds.update({x: [x+2], x+2: [x+3, x+4, x+5], x+3: [x+3], x+4: [x+4], x+5: [x+5]})
        nrot_bonds = [(x, x+2)]
        return nrot_bonds

    def SER(x):
        # CA, CB, HB1, HB2, OG, HG1
        bonds.update({x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 4], x + 5: [x + 6],
                      x + 6: [x + 6]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5)]
        return nrot_bonds

    def THR(x):
        # CA, CB, HB, OG1, HG1, CG2, HG21, HG22, HG23
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 6], x + 3: [x + 3], x + 4: [x + 5], x + 5: [x + 5],
             x + 6: [x + 7, x + 8, x + 9], x + 7: [x + 7], x+8: [x+8], x+9: [x+9]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 4), (x + 2, x + 6)]
        return nrot_bonds

    def CYS(x):
        # CA, CB, HB1, HB2, SG, HG1
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 4], x + 5: [x + 6],
             x + 6: [x + 6]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5)]
        return nrot_bonds

    def VAL(x):
        # CA, CB, HB, CG1, HG11, HG12, HG13, CG2, HG21, HG22, HG23
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 8], x + 3: [x + 3], x + 4: [x + 5, x+6, x+7],
             x + 5: [x + 5], x + 6: [x+6], x + 7: [x+7], x+8: [x+9, x+10, x+11], x+9: [x+9], x+10: [x+10],
             x+11: [x+11]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 4), (x + 2, x + 8)]
        return nrot_bonds

    def LEU(x):
        # CA, CB, HB1, HB2, CG, HG, CD1, HD11, HD12, HD13, CD2, HD21, HD22, HD23
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 4],
             x + 5: [x + 6, x + 7, x + 11], x + 6: [x + 6], x + 7: [x + 8, x + 9, x + 10], x + 8: [x + 8],
             x + 9: [x + 9], x + 10: [x + 10], x + 11: [x + 12, x+13, x+14],
             x+12: [x+12], x+13: [x+13], x+14: [x+14]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5), (x + 5, x + 7), (x + 5, x + 11)]
        return nrot_bonds

    def ILE(x):
        # CA, CB, HB, CG2, HG21, HG22, HG23, CG1, HG11, HG12, CD, HD1, HD2, HD3
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 8], x + 3: [x + 3], x + 4: [x + 5, x + 6, x + 7],
             x + 5: [x + 5], x + 6: [x + 6], x + 7: [x + 7], x + 8: [x + 9, x + 10, x + 11],
             x + 9: [x + 9], x + 10: [x + 10], x + 11: [x + 12, x + 13, x + 14],
             x + 12: [x + 12], x + 13: [x + 13], x + 14: [x + 14]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 4), (x + 2, x + 8), (x + 8, x + 11)]
        return nrot_bonds

    def MET(x):
        # CA, CB, HB1, HB2, CG, HG1, HG2, SD, CE, HE1, HE2, HE3
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 4],
             x + 5: [x + 6, x + 7, x + 8], x + 6: [x + 6], x + 7: [x + 7], x + 8: [x + 9],
             x + 9: [x + 10, x + 11, x+12], x + 10: [x + 10], x + 11: [x + 11],
             x + 12: [x+12]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5), (x + 5, x + 8), (x + 8, x + 9)]
        return nrot_bonds

    def PRO(x):
        return []

    def PHE(x):
        # CA, CB, HB1, HB2, CG, CD1, HD1, CE1, HE1, CZ, HZ, CD2, HD2, CE2, HE2
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 4],
             x + 5: [x + 6, x + 8], x + 6: [x + 7, x + 12], x + 7: [x + 7], x + 8: [x + 9, x+14],
             x + 9: [x + 9], x + 10: [x + 11], x + 11: [x + 11],
             x+12: [x+10, x+13], x+13: [x+13], x+14: [x+10, x+15], x+15: [x + 15]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5)]
        return nrot_bonds

    def TYR(x):
        # CA, CB, HB1, HB2, CG, CD1, HD1, CE1, HE1, CZ, OH, HH, CD2, HD2, CE2, HE2
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 4],
             x + 5: [x + 6, x + 8], x + 6: [x + 7, x + 13], x + 7: [x + 7], x + 8: [x + 9, x+15],
             x + 9: [x + 10], x + 10: [x + 11], x + 11: [x + 12],
             x+12: [x+12], x+13: [x+10, x+14], x+14: [x+14], x+15: [x + 10, x + 16], x+16: [x+16]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5), (x + 2, x + 6), (x + 10, x + 11)]
        return nrot_bonds

    def TRP(x):
        # CA, CB, HB1, HB2, CG, CD1, HD1, NE1, HE1, CE2, CD2, CE3, HE3, CZ3, HZ3, CZ2, HZ2, CH2, HH2
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 4],
             x + 5: [x + 6], x + 6: [x + 7, x + 18], x + 7: [x + 8, x + 9], x + 8: [x + 8],
             x + 9: [x + 10, x + 11], x + 10: [x + 10], x + 11: [x + 12, x + 18],
             x+12: [x+13, x+14], x+13: [x+13], x+14: [x+15, x+16], x+15: [x + 15], x+16: [x+17, x+18], x+17: [x+17],
             x+18: [x+19], x+19: [x+19]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5)]
        return nrot_bonds

    def ASP(x):
        # CA, CB, HB1, HB2, CG, OD1, OD2
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 4],
             x + 5: [x + 6, x + 7], x + 6: [x + 6], x + 7: [x + 7]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5)]
        return nrot_bonds

    def GLU(x):
        # CA, CB, HB1, HB2, CG, HG1, HG2, CD, OE1, OE2
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 4],
             x + 5: [x + 6, x + 7, x + 8], x + 6: [x + 6], x + 7: [x + 7], x + 8: [x + 9, x + 10],
             x + 9: [x + 9], x + 10: [x + 10]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5), (x + 5, x + 8)]
        return nrot_bonds

    def ASN(x):
        # CA, CB, HB1, HB2, CG, OD1, ND2, HD21, HD22
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 4],
             x + 5: [x + 6, x + 7], x + 6: [x + 6], x + 7: [x + 8, x + 9], x + 8: [x + 8], x + 9: [x + 9]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5), (x + 5, x + 7)]
        return nrot_bonds

    def GLN(x):
        return []

    def HIS(x):
        return []

    def LYS(x):
        # CA, CB, HB1, HB2, CG, HG1, HG2, CD, HD1, HD2, CE, HE1, HE2, NZ, NZ1, NZ2, NZ3
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 4],
             x + 5: [x + 6, x + 7, x + 8], x + 6: [x + 6], x + 7: [x + 7], x + 8: [x + 9, x + 10, x + 11],
             x + 9: [x + 9], x + 10: [x + 10], x + 11: [x + 12, x + 13, x + 14],
             x+12: [x+12], x + 13: [x + 13], x + 14: [x+15, x + 16, x + 17], x + 15: [x + 15], x + 16: [x + 16],
             x + 17: [x + 17]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5), (x + 5, x + 8), (x + 8, x + 11), (x + 11, x + 14)]
        return nrot_bonds

    def ARG(x):
        # CA, CB, HB1, HB2, CG, HG1, HG2, CD, HD1, HD2, NE, HE, CZ, NH1, HH11, HH12, NH2, HH21, HH22
        bonds.update(
            {x: [x + 2], x + 2: [x + 3, x + 4, x + 5], x + 3: [x + 3], x + 4: [x + 5],
             x + 5: [x + 6, x + 7, x + 8], x + 6: [x + 6], x + 7: [x + 7], x + 8: [x + 9, x + 10, x + 11],
             x + 9: [x + 9], x + 10: [x + 10], x + 11: [x + 12, x + 13],
             x + 12: [x + 12], x + 13: [x + 14, x + 17], x + 14: [x + 15, x + 16], x + 15: [x + 15], x + 16: [x + 16],
             x + 17: [x + 18, x + 19], x + 18: [x + 18], x + 19: [x + 19]})
        nrot_bonds = [(x, x + 2), (x + 2, x + 5), (x + 5, x + 8), (x + 8, x + 11), (x + 11, x + 13),
                     (x + 13, x + 14), (x + 13, x + 17)]
        return nrot_bonds

    def side_chains(ResName, num):
        return {'GLY': lambda: GLY(num),
                'ALA': lambda: ALA(num),
                'SER': lambda: SER(num),
                'THR': lambda: THR(num),
                'CYS': lambda: CYS(num),
                'VAL': lambda: VAL(num),
                'LEU': lambda: LEU(num),
                'ILE': lambda: ILE(num),
                'MET': lambda: MET(num),
                'PHE': lambda: PHE(num),
                'PRO': lambda: PRO(num),
                'TYR': lambda: TYR(num),
                'TRP': lambda: TRP(num),
                'ASP': lambda: ASP(num),
                'GLU': lambda: GLU(num),
                'ASN': lambda: ASN(num),
                'GLN': lambda: GLN(num),
                'HIS': lambda: HIS(num),
                'LYS': lambda: LYS(num),
                'ARG': lambda: ARG(num)}.get(ResName, lambda: [])()

    for i, a in enumerate(mol.values()):
        if a.Name == 'CA' and a.ResSeq in rotating_resid:
            nrot_bonds = side_chains(a.ResName, i+1)
            rot_bonds += nrot_bonds
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
    d = math.sqrt(norm_vector[1]**2+norm_vector[2]**2)
    
    def count_matrix(origin_point, norm_vector, angle):
        trans_matrix = np.array([[1.0, 0.0, 0.0, 0.0],
                                [0.0, 1.0, 0.0, 0.0],
                               [0.0, 0.0, 1.0, 0.0],
                               [-origin_point[0], -origin_point[1], -origin_point[2], 1.0]])
        rot_x_matrix = np.array([[1.0, 0.0, 0.0, 0.0],
                                   [0.0, norm_vector[2]/d, norm_vector[1]/d, 0.0],
                                   [0.0, -norm_vector[1]/d, norm_vector[2]/d, 0.0],
                                   [0.0, 0.0, 0.0, 1.0]])
        rot_y_matrix = np.array([[d, 0.0, norm_vector[0], 0.0],
                                   [0.0, 1.0, 0.0, 0.0],
                                   [-norm_vector[0], 0.0, d, 0.0],
                                   [0.0, 0.0, 0.0, 1.0]])
        a = np.radians(angle)
        rot_a_matrix = np.array([[np.cos(a), np.sin(a), 0.0, 0.0],
                                   [-np.sin(a), np.cos(a), 0.0, 0.0],
                                   [0.0, 0.0, 1.0, 0.0],
                                   [0.0, 0.0, 0.0, 1.0]])
        return trans_matrix, rot_x_matrix, rot_y_matrix, rot_a_matrix

    trans_matrix, rot_x_matrix, rot_y_matrix, rot_a_matrix = count_matrix(origin_point1, norm_vector, angle)
    matrix = np.dot(np.dot(trans_matrix, rot_x_matrix), rot_y_matrix)

    change_point = np.append(point, 1)
    inv_matrix = np.copy(matrix)
    inv_matrix = np.transpose(inv_matrix)
    inv_matrix[:, 3] = [0.0, 0.0, 0.0, 1.0]
    inv_matrix[3] = list(origin_point1) + [1]
    final_matrix = np.dot(np.dot(matrix, rot_a_matrix), inv_matrix)
    new_point = np.dot(change_point, final_matrix)
    return np.delete(new_point, 3)


def write_result(fname, rotations, best_energy):
    with open(fname, "w") as f:
        for i in rotations:
            f.write(f"Rotation around an axis: {i[0]}, {i[1]} by an angle: {i[2]}\n")

        f.write(f"The best energy: {best_energy}")

