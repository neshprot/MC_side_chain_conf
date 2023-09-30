"""
Main script for MC sim
"""
from utils import math, rotate, write_pdb
from logger import FileLogger
import random


def montecarlo(mol, graph, rot_bonds, attempts, stop_step, rotating_resid):
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
    iterations = 0  # const
    step = 0    # const
    energy_dict = {}
    res_dict = {}
    k = 1.987204259e-3  # Boltzmann constant
    t1 = 100  # first temperature
    t2 = 4000   # second temperature 1000
    t = t1
    eps_amen = 0.4   # amendment of epsilon 0.4
    eps = 1     # standard amendment of epsilon
    best_energy = 'None'
    T_count = 0     # const
    T_step = 0      # const
    T_lim = 100   # every T_count steps temperature changes
    eps_step = 0    # const
    eps_lim = 100   # every eps_lim steps epsilon changes

    ILE_prob = 1

    TRP_queue_update = 50
    TRP_queue_step = 0  # const
    angle_TRP = 0   # const
    TRP_lst = []    # const

    charge_dict = dict()

    all_nums = []

    start_was = False


    def Energy(mol, coords_for_rot, eps_amen):
        '''
        Energy function, based on CHARMM

        :param mol: protein
        :param coords_for_rot: list of coords that were rotated
        :return: energy of the side chain that were rotated
        '''

        energy_C = 0
        energy_LJ = 0

        def Coulomb(mol1, mol2, distance):
            '''
            Coulomb energy function between two different atoms in protein

            :param mol1: num of the first atom
            :param mol2: num of the se
            :param distance: distance between two atoms (angstrom)
            :return:
            '''
            first_charge = mol[mol1].partialcharge
            second_charge = mol[mol2].partialcharge
            result = first_charge * second_charge / distance * 332
            return result

        def LJ(mol1, mol2, distance):
            '''
            Lennard-Jones energy function between two different atoms in protein
            :param mol1: num of the first atom
            :param mol2: num of the se
            :param distance: distance between two atoms (angstrom)
            :return:
            '''
            epsilon12 = math.sqrt(mol[mol1].epsilon * mol[mol2].epsilon)
            rmin12 = mol[mol1].rmin + mol[mol2].rmin
            result = epsilon12 * ((rmin12 / distance) ** 12 - 2 * (rmin12 / distance) ** 6)
            return result

        # count the potencial between each atom of the aminoacid
        # and the atoms at the distance of 10 angstroms
        for mol1 in coords_for_rot:
            for mol2 in mol:
                if mol1 not in mol[mol2].bonded:

                    # square distance between two atoms
                    distance = sum((x - y) ** 2 for x, y in 
                                   zip(mol[mol1].coordin, mol[mol2].coordin))

                    if mol1 == mol2:
                        continue
                    elif distance <= 100:
                        sqrt_distance = math.sqrt(distance)  # distance between two atoms
                        energy_C_bet = Coulomb(mol1, mol2, sqrt_distance)
                        energy_LJ_bet = LJ(mol1, mol2, sqrt_distance)
                        energy_C += energy_C_bet
                        energy_LJ += energy_LJ_bet

                        if start_was and (mol2 in all_nums) and (mol1, mol2) in charge_dict:
                            old_energy = energy_dict[mol[mol2].resseq] + \
                                         energy_C_bet + energy_LJ_bet - charge_dict[
                                (mol1, mol2)]
                            energy_dict.update({mol[mol2].resseq: old_energy})

                        charge_dict.update({(mol1, mol2): energy_C_bet + energy_LJ_bet})
                        charge_dict.update({(mol2, mol1): energy_C_bet + energy_LJ_bet})
        return energy_C, eps_amen * energy_LJ

    # evaluate start energy
    def start_energy():
        summ = 0
        TRP_count = 0
        ini_rot_bonds = rot_bonds
        for atom in mol:
            if mol[atom].resseq not in rotating_resid:
                continue
            if any(atom in i for i in ini_rot_bonds) and mol[atom].name == 'CA':
                atom_lst = graph.bfs(atom)
                energy_C, energy_LJ = Energy(mol, atom_lst, eps_amen=1)
                res_energy = energy_C + energy_LJ
                # record the energy of each aminoacid
                energy_dict.update({mol[atom].resseq: res_energy})
                # record all atoms in a particular amino acid
                res_dict.update({mol[atom].resseq: atom_lst})
                summ += res_energy
                all_nums.append(atom)

    #special function for ILE
    def ILE_rot(bond, mol):

        ile_angle1 = random.choice([30, 60, 90, 180]) * random.choice([-1, 1])
        ile_angle2 = random.choice([30, 60, 90, 180]) * random.choice([-1, 1])

        coords_for_rot1 = graph.bfs(bond[1])
        coords_for_rot2 = graph.bfs(bond[1]+6)

        rotate(mol, coords_for_rot1, bond[0], bond[1], ile_angle1)
        rotate(mol, coords_for_rot2, bond[1], bond[1] + 6, ile_angle2)
        return [ile_angle1, ile_angle2], coords_for_rot, initial_coords

    start_energy()

    start_was = True

    # main block with consecutive rotations
    while iterations <= attempts and step < stop_step:

        bond = random.choice(rot_bonds)     # selecting a random bond for rotation

        first_prob = 1
        second_prob = 0

        angle1 = random.randint(0, 10)
        angle2 = random.randint(90, 180)
        angle_sign = random.choice([-1, 1])
        angle = random.choices([angle1, angle2], weights=[first_prob, second_prob])[0] * angle_sign

        step += 1
        T_step += 1
        eps_step += 1

        # temperature changes every T_lim steps
        if T_step == T_lim:
            T_step = 0
            if T_count == 0:
                T_count = 1
                t = t2
            else:
                T_count = 0
                t = t1

        # epsilon change every eps_lim steps
        if eps_step == eps_lim:
            eps_step = 0
            eps = eps_amen

        coords_for_rot = graph.bfs(bond[1])
        initial_coords = [mol[i].coordin for i in coords_for_rot]

        rotate(mol, coords_for_rot, bond[0], bond[1], angle)

        # evaluate new energy and compare it to old one
        energy_C, energy_LJ = Energy(mol, res_dict[mol[coords_for_rot[0]].resseq], eps)
        new_energy = energy_C + energy_LJ
        old_energy = energy_dict[mol[coords_for_rot[0]].resseq]

        # whether the system has moved to a higher energy step
        def transition_probability(new_energy, old_energy):
            if new_energy <= old_energy:
                return True
            else:
                prob_up = math.exp(-(new_energy - old_energy)/(k*t))
                if random.uniform(0.0, 1.0) <= prob_up:
                    return True
                else:
                    return False

        choice = transition_probability(new_energy, old_energy)

        if choice:
            iterations = 0
            energy_dict[mol[coords_for_rot[0]].resseq] = energy_C + energy_LJ/eps
            # use next line if you want to write every step in frame
            #write_pdb(mol, f'for_frames/out_{mol[coords_for_rot[0]].resseq}_{step}.pdb')
            logger(f"Step: {step}. Good current rotation:\n"
                   f"First num: {bond[0]}, second num: {bond[1]},"
                   f" resid: {mol[coords_for_rot[0]].resseq},"
                   f" angle: {angle}, energy: {new_energy}\n\n")
        else:
            iterations += 1
            for i, num in enumerate(coords_for_rot):
                mol[num].coordin = initial_coords[i]

            logger(f"Step: {step}. Bad current rotation:\n"
                   f"First num: {bond[0]}, second num: {bond[1]},"
                   f" resid: {mol[coords_for_rot[0]].resseq},"
                   f" angle: {angle}, energy: {new_energy}\n\n")

        eps = 1

    write_pdb(mol, 'out.pdb')
    logger(f"The best energy: {best_energy}\n"
           f"Step/Stop {step}/{stop_step}\n")
    logger("\n")
    return rotations, best_energy
