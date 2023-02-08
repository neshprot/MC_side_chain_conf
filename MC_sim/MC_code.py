from utils import *
from logger import FileLogger
import random


def Energy(mol, rotating_resid):
    energy = 0

    def Coulomb(mol1, mol2, distance):
        first_charge = mol[mol1].PartialCharge
        second_charge = mol[mol2].PartialCharge
        result = first_charge*second_charge/distance
        return result

    def LJ(mol1, mol2, distance):
        epsilon12 = math.sqrt(mol[mol1].Epsilon*mol[mol2].Epsilon)
        rmin12 = mol[mol1].Rmin + mol[mol2].Rmin
        result = epsilon12*((rmin12/distance)**12 - 2*(rmin12/distance)**6)
        return result

    for mol1 in mol:
        for mol2 in rotating_resid:
            if mol1 != mol2:

                if mol1 not in mol[mol2].Bonded:
                    distance = sum((x - y) ** 2 for x, y in zip(mol[mol1].Coordin, mol[mol2].Coordin))
                    energy += Coulomb(mol1, mol2, distance) + LJ(mol1, mol2, distance)
    return energy


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
        new_coord = rotation(mol[bond_num1].Coordin, mol[bond_num2].Coordin, mol[i].Coordin, angle)
        mol[i].Coordin = new_coord


def MonteCarlo(mol, graph, rot_bonds, start_energy, attempts, stop_step, rotating_resid):

    logger = FileLogger("logout")
    rotations = []
    iterations = 0
    step = 0
    best_energy = start_energy

    while iterations <= attempts and step < stop_step:
        bond = random.choice(rot_bonds)
        angle = random.uniform(-10, 10)
        coords_for_rot = graph.bfs(bond[1])
        initial_coords = [mol[i].Coordin for i in coords_for_rot]
        rotate(mol, coords_for_rot, bond[0], bond[1], angle)
        energy = Energy(mol, rotating_resid)
        step += 1
        if energy >= best_energy:
            iterations += 1
            for i, num in enumerate(coords_for_rot):
                mol[num].Coordin = initial_coords[i]

            logger(f"Bad current rotation:\n"
                   f"Firt num: {bond[0]}, seecond num: {bond[1]}, angle: {angle}, energy: {energy}\n\n")
        else:
            rotations += [(bond[0], bond[1], angle)]
            iterations = 0
            best_energy = energy

            logger(f"Good current rotation:\n"
                   f"Firt num: {bond[0]}, seecond num: {bond[1]}, angle: {angle}, energy: {energy}\n\n")

    write_pdb(mol, 'out.pdb')
    logger(f"The best energy: {best_energy}\n"
           f"Step/Stop {step}/{stop_step}\n")
    logger("\n")
    return rotations, best_energy
