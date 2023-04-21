from utils import *
from logger import FileLogger
import random

import time


def Energy(mol, coords_for_rot, eps_amen):
    '''
    Energy function, based on CHARMM

    :param mol: protein
    :param coords_for_rot: list of coords that were rotated
    :return: energy of the side chain that were rotated
    '''
    energy = 0
    energy_C = 0
    energy_LJ = 0

    def Coulomb(mol1, mol2, distance):
        first_charge = mol[mol1].PartialCharge
        second_charge = mol[mol2].PartialCharge
        result = first_charge*second_charge/distance*332
        return result

    def LJ(mol1, mol2, distance):
        epsilon12 = math.sqrt(mol[mol1].Epsilon*mol[mol2].Epsilon)
        rmin12 = mol[mol1].Rmin + mol[mol2].Rmin
        #result = eps_amen*epsilon12*((rmin12/distance)**12 - 2*(rmin12/distance)**6)
        result = epsilon12*((rmin12/distance)**12 - 2*(rmin12/distance)**6)
        return result

    for mol1 in coords_for_rot:
        for mol2 in mol:
            if mol1 not in mol[mol2].Bonded:
                # distance between two atoms
                distance = sum((x - y) ** 2 for x, y in zip(mol[mol1].Coordin, mol[mol2].Coordin))

                if mol1 == mol2:
                    continue
                elif distance <= 100:
                    energy_C += Coulomb(mol1, mol2, math.sqrt(distance))
                    energy_LJ += LJ(mol1, mol2, math.sqrt(distance))
                    #energy += LJ(mol1, mol2, math.sqrt(distance), eps_amen) + Coulomb(mol1, mol2, math.sqrt(distance))
    return energy_C, eps_amen*energy_LJ


def MonteCarlo(mol, graph, rot_bonds, attempts, stop_step, rotating_resid):
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
    k = 1.987204259e-3  # Boltzmann constant
    T1 = 200  # first temperature
    T2 = 2000   # second temperature
    T = T1
    eps_amen = 0.01   # amendment of epsilon
    eps = 1     # standard amendment of epsilon
    best_energy = 'None'
    T_count = 0
    T_step = 0
    eps_step = 0

    # evaluate start energy
    def start_energy():
        ini_rot_bonds = rot_bonds
        for atom in mol:
            if mol[atom].ResSeq not in rotating_resid:
                continue
            if any(atom in i for i in ini_rot_bonds) and mol[atom].Name == 'CA':
                atom_lst = graph.bfs(atom)
                #res_energy = Energy(mol, atom_lst, eps_amen=1)
                energy_C, energy_LJ = Energy(mol, atom_lst, eps_amen=1)
                res_energy = energy_C + energy_LJ
                energy_dict.update({mol[atom].ResSeq: res_energy})
                res_dict.update({mol[atom].ResSeq: atom_lst})

    start_energy()

    rot_bonds11 = (2184, 2187)
    rot_bonds12 = (2182, 2184)

    rot_bonds31 = (1361, 1363)
    rot_bonds32 = (1363, 1366)

    rot_bonds21 = (3007, 3010)
    rot_bonds22 = (3005, 3007)

    variety = [40, 35, 20, 5]
    counter = 0

    # main block with consecutive rotations
    while iterations <= attempts and step < stop_step:

        f_time = time.time()

        #bond = random.choice(rot_bonds)

        bond4 = random.choice(rot_bonds)
        pos = random.choice([0, 1])
        if counter == 1000:
            variety = variety[-1:] + variety[:-1]
            counter = 0
        counter += 1
        bond = random.choices([(rot_bonds11, rot_bonds12), (rot_bonds21, rot_bonds22), (rot_bonds31, rot_bonds32), (bond4, bond4)],
                             weights=variety)[0][pos]

        angle_sign = random.choice([-1, 1])

        # firstly we want to find place in space where amino acid have lowest energy, after that we clarify postion
        if step < stop_step*0.5:
            first_prob = 30
            second_prob = 70
        else:
            first_prob = 80
            second_prob = 20

        angle1 = random.uniform(0, 10)
        angle2 = random.uniform(90, 180)
        angle = random.choices([angle1, angle2], weights=[first_prob, second_prob])[0] * angle_sign

        step += 1
        T_step += 1
        eps_step += 1

        # temperature change
        if T_step == 100:
            T_step = 0
            if T_count == 0:
                T_count = 1
                T = T2
            else:
                T_count = 0
                T = T1

        # epsilon change
        if eps_step == 50:
            eps_step = 0
            eps = eps_amen

        coords_for_rot = graph.bfs(bond[1])
        initial_coords = [mol[i].Coordin for i in coords_for_rot]
        rotate(mol, coords_for_rot, bond[0], bond[1], angle)
        #new_energy = Energy(mol, res_dict[mol[coords_for_rot[0]].ResSeq], eps)
        energy_C, energy_LJ = Energy(mol, res_dict[mol[coords_for_rot[0]].ResSeq], eps)
        new_energy = energy_C + energy_LJ
        old_energy = energy_dict[mol[coords_for_rot[0]].ResSeq]

        # whether the system has moved to a higher energy step
        def transition_probability(new_energy, old_energy):
            if new_energy <= old_energy:
                return True
            else:
                prob_up = math.exp(-(new_energy - old_energy)/(k*T))
                if random.uniform(0.0, 1.0) <= prob_up:
                    return True
                else:
                    return False

        choice = transition_probability(new_energy, old_energy)
        if choice:
            rotations += [(bond[0], bond[1], angle)]
            iterations = 0
            energy_dict[mol[coords_for_rot[0]].ResSeq] = energy_C + energy_LJ/eps
            write_pdb(mol, f'for_frames/out_{mol[coords_for_rot[0]].ResSeq}_{step}.pdb')
            logger(f"Step: {step}. Good current rotation:\n"
                   f"First num: {bond[0]}, second num: {bond[1]}, resid: {mol[coords_for_rot[0]].ResSeq},"
                   f" angle: {angle}, energy: {new_energy}\n\n")
        else:
            iterations += 1
            for i, num in enumerate(coords_for_rot):
                mol[num].Coordin = initial_coords[i]

            logger(f"Step: {step}. Bad current rotation:\n"
                   f"First num: {bond[0]}, second num: {bond[1]}, resid: {mol[coords_for_rot[0]].ResSeq},"
                   f" angle: {angle}, energy: {new_energy}\n\n")
            write_pdb(mol, 'out.pdb')

        eps = 1

        print(time.time() - f_time)

    logger(f"The best energy: {best_energy}\n"
           f"Step/Stop {step}/{stop_step}\n")
    logger("\n")
    return rotations, best_energy
