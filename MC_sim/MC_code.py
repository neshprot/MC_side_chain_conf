from utils import *
from logger import FileLogger
import random


def Energy(mol, coords_for_rot):
    '''
    Energy function, based on CHARMM

    :param mol: protein
    :param coords_for_rot: list of coords that were rotated
    :return: energy of the side chain that were rotated
    '''
    energy = 0

    def Coulomb(mol1, mol2, distance):
        first_charge = mol[mol1].PartialCharge
        second_charge = mol[mol2].PartialCharge
        result = first_charge*second_charge/(4*math.pi*distance)
        return result

    def LJ(mol1, mol2, distance):
        epsilon12 = math.sqrt(mol[mol1].Epsilon*mol[mol2].Epsilon)
        rmin12 = mol[mol1].Rmin + mol[mol2].Rmin
        result = epsilon12*((rmin12/distance)**12 - 2*(rmin12/distance)**6)
        return result

    for mol1 in coords_for_rot:
        for mol2 in mol:
            if mol1 not in mol[mol2].Bonded:
                distance = sum((x - y) ** 2 for x, y in zip(mol[mol1].Coordin, mol[mol2].Coordin))
                if mol1 == mol2:
                    continue
                elif 0 <= distance <= 1:
                    return 1000
                elif 0 < distance <= 100:
                    energy += Coulomb(mol1, mol2, math.sqrt(distance)) + LJ(mol1, mol2, math.sqrt(distance))
    return energy


def MonteCarlo(mol, graph, rot_bonds, attempts, stop_step, rotating_resid, rot_bonds_CA):
    '''
    main block or MC algorythm

    :param mol: protein
    :param graph: class that describes bonds
    :param rot_bonds: bonds for rotation
    :param attempts: amount of attempts to find good rotation
    :param stop_step: amount of maximum rotations
    :param rotating_resid: id of residue for rotation
    :param rot_bonds_CA: bonds with CA atom
    :return: list of good rotations and best energy
    '''

    logger = FileLogger("logout")
    rotations = []
    iterations = 0
    step = 0
    energy_dict = dict()
    res_dict = dict()
    best_energy = 'None'

    def start_energy():
        ini_rot_bonds = rot_bonds
        for atom in mol:
            if mol[atom].ResSeq not in rotating_resid:
                continue
            if any(atom in i for i in ini_rot_bonds) and mol[atom].Name == 'CA':
                atom_lst = graph.bfs(atom)
                res_energy = Energy(mol, atom_lst)
                energy_dict.update({mol[atom].ResSeq: res_energy})
                res_dict.update({mol[atom].ResSeq: atom_lst})

    start_energy()

    while iterations <= attempts and step < stop_step:
        '''        
        if step <= stop_step/6:
            bond = random.choice(rot_bonds_CA)
        else:
            bond = random.choice(rot_bonds)
        '''
        bond = random.choice(rot_bonds)
        angle = random.uniform(-180, 180)
        coords_for_rot = graph.bfs(bond[1])
        initial_coords = [mol[i].Coordin for i in coords_for_rot]
        rotate(mol, coords_for_rot, bond[0], bond[1], angle)
        energy = Energy(mol, res_dict[mol[coords_for_rot[0]].ResSeq])
        step += 1
        if energy >= energy_dict[mol[coords_for_rot[0]].ResSeq]:
            iterations += 1
            for i, num in enumerate(coords_for_rot):
                mol[num].Coordin = initial_coords[i]

            logger(f"Bad current rotation:\n"
                   f"First num: {bond[0]}, second num: {bond[1]}, resid: {mol[coords_for_rot[0]].ResSeq}, angle: {angle}, energy: {energy}\n\n")
        else:
            rotations += [(bond[0], bond[1], angle)]
            iterations = 0
            energy_dict[mol[coords_for_rot[0]].ResSeq] = energy

            logger(f"Good current rotation:\n"
                   f"First num: {bond[0]}, second num: {bond[1]}, resid: {mol[coords_for_rot[0]].ResSeq}, angle: {angle}, energy: {energy}\n\n")

    write_pdb(mol, 'out.pdb')
    logger(f"The best energy: {best_energy}\n"
           f"Step/Stop {step}/{stop_step}\n")
    logger("\n")
    return rotations, best_energy
