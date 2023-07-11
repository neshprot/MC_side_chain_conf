from utils import *
from logger import FileLogger
import random



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
    T1 = 100  # first temperature
    T2 = 1000   # second temperature
    T = T1
    eps_amen = 0.8   # amendment of epsilon
    eps = 1     # standard amendment of epsilon
    best_energy = 'None'
    T_count = 0
    T_step = 0
    T_lim = 100   # every T_count steps temperature changes
    eps_step = 0
    eps_lim = 100   # every eps_lim steps epsilon changes

    TRP_prob = 0.9
    TRP_queue_update = 50
    TRP_queue_step = 0
    angle_TRP = 0
    TRP_lst = []

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
            first_charge = mol[mol1].PartialCharge
            second_charge = mol[mol2].PartialCharge
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
            epsilon12 = math.sqrt(mol[mol1].Epsilon * mol[mol2].Epsilon)
            rmin12 = mol[mol1].Rmin + mol[mol2].Rmin
            result = epsilon12 * ((rmin12 / distance) ** 12 - 2 * (rmin12 / distance) ** 6)
            return result

        # count the potencial between each atom of the aminoacid and the atoms at the distance of 10 angstroms
        for mol1 in coords_for_rot:
            for mol2 in mol:
                if mol1 not in mol[mol2].Bonded:

                    # square distance between two atoms
                    distance = sum((x - y) ** 2 for x, y in zip(mol[mol1].Coordin, mol[mol2].Coordin))

                    if mol1 == mol2:
                        continue
                    elif distance <= 100:
                        sqrt_distance = math.sqrt(distance)  # distance between two atoms
                        energy_C_bet = Coulomb(mol1, mol2, sqrt_distance)
                        energy_LJ_bet = LJ(mol1, mol2, sqrt_distance)
                        energy_C += energy_C_bet
                        energy_LJ += energy_LJ_bet

                        if start_was and (mol2 in all_nums) and (mol1, mol2) in charge_dict:
                            old_energy = energy_dict[mol[mol2].ResSeq] + energy_C_bet + energy_LJ_bet - charge_dict[
                                (mol1, mol2)]
                            energy_dict.update({mol[mol2].ResSeq: old_energy})

                        charge_dict.update({(mol1, mol2): energy_C_bet + energy_LJ_bet})
                        charge_dict.update({(mol2, mol1): energy_C_bet + energy_LJ_bet})

        return energy_C, eps_amen * energy_LJ

    # special function for TRP
    def TRP_rot(bond, mol):
        coords_for_rot1 = graph.bfs(bond[1])
        initial_coords = [mol[i].Coordin for i in coords_for_rot1]

        '''
        angle1 = random.choice([45, 90, 135, 180])*random.choice([1, -1])     # rotation angle around CA and CB
        angle2 = random.choice([45, 90, 135, 180])*random.choice([1, -1])     # rotation angle around CB and CG
        ex_angle1 = random.uniform(-10.0, 10.0)     # small addition to the large angle1
        ex_angle2 = random.uniform(-10.0, 10.0)     # small addition to the large angle2
        trp_angle1 = angle1 + ex_angle1
        trp_angle2 = angle2 + ex_angle2
        '''
        first_prob1 = 80
        second_prob1 = 20

        #angle11 = random.uniform(0, 50)
        angle11 = random.randint(0, 50)
        #angle12 = random.uniform(90, 180)
        angle12 = random.randint(90, 180)
        angle_sign1 = random.choice([-1, 1])
        trp_angle1 = random.choices([angle11, angle12], weights=[first_prob1, second_prob1])[0] * angle_sign1

        first_prob2 = 20
        second_prob2 = 80

        #angle21 = random.uniform(0, 50)
        angle21 = random.randint(0, 50)
        #angle22 = random.uniform(90, 180)
        angle22 = random.randint(0, 50)
        angle_sign2 = random.choice([-1, 1])
        trp_angle2 = random.choices([angle21, angle22], weights=[first_prob2, second_prob2])[0] * angle_sign2

        rotate(mol, coords_for_rot1, bond[0], bond[1], trp_angle1)
        coords_for_rot2 = graph.bfs((bond[1] + 3))
        rotate(mol, coords_for_rot2, bond[1], bond[1]+3, trp_angle2)
        return trp_angle1, coords_for_rot1, initial_coords, trp_angle2

    # evaluate start energy
    def start_energy():
        summ = 0
        TRP_count = 0
        ini_rot_bonds = rot_bonds
        for atom in mol:
            if mol[atom].ResSeq not in rotating_resid:
                continue
            if any(atom in i for i in ini_rot_bonds) and mol[atom].Name == 'CA':
                atom_lst = graph.bfs(atom)
                energy_C, energy_LJ = Energy(mol, atom_lst, eps_amen=1)
                res_energy = energy_C + energy_LJ
                # record the energy of each aminoacid
                energy_dict.update({mol[atom].ResSeq: res_energy})
                # record all atoms in a particular amino acid
                res_dict.update({mol[atom].ResSeq: atom_lst})
                summ += res_energy
                all_nums.append(atom)
                if mol[atom].ResName == 'TRP':
                    TRP_count += 1
        return TRP_count

    '''
    #I want to rotate only TRP in 6GUX_t
    # 141
    coords_for_rot = graph.bfs(2165)
    rotate(mol, coords_for_rot, 2163, 2165, 74)
    rotate(mol, coords_for_rot, 2165, 2168, 195)

    # 192
    coords_for_rot = graph.bfs(2988)
    rotate(mol, coords_for_rot, 2986, 2988, 180)
    rotate(mol, coords_for_rot, 2988, 2991, 90)

    # 15
    coords_for_rot = graph.bfs(207)
    rotate(mol, coords_for_rot, 205, 207, 200)
    rotate(mol, coords_for_rot, 207, 210, 175)
    
    # 89
    coords_for_rot = graph.bfs(207)
    rotate(mol, coords_for_rot, 1360, 1362, 35)
    rotate(mol, coords_for_rot, 1362, 1365, 179)
    
    # 141
    coords_for_rot = graph.bfs(2184)
    rotate(mol, coords_for_rot, 2182, 2184, -45)
    rotate(mol, coords_for_rot, 2184, 2187, 210)
    write_pdb(mol, 'ini.pdb')
    '''


    TRP_count = start_energy()

    start_was = True

    # main block with consecutive rotations
    while iterations <= attempts and step < stop_step:

        bond = random.choice(rot_bonds)     # selecting a random bond for rotation

        '''
        # firstly we want to find place in space where amino acid have lowest energy, after that we clarify postion
        if step < stop_step*0.2:
            first_prob = 20     # probability that angle will be between 0 and 10 degrees
            second_prob = 80    # probability that angle will be between 90 and 180 degrees
        else:
            first_prob = 80
            second_prob = 20
        '''
        first_prob = 30
        second_prob = 70

        #angle1 = random.uniform(0, 10)
        angle1 = random.randint(0, 10)
        #angle2 = random.uniform(90, 180)
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
                T = T2
            else:
                T_count = 0
                T = T1

        # epsilon change every eps_lim steps
        if eps_step == eps_lim:
            eps_step = 0
            eps = eps_amen

        '''
        # rotate amino acid and memorize initial position
        coords_for_rot = graph.bfs(bond[1])
        initial_coords = [mol[i].Coordin for i in coords_for_rot]
        rotate(mol, coords_for_rot, bond[0], bond[1], angle)
        '''

        # do the tryptophan function (if we use TRP function? we have to remove previous block with rotation)
        coords_for_rot = graph.bfs(bond[1])
        initial_coords = [mol[i].Coordin for i in coords_for_rot]

        if mol[bond[0]].ResSeq not in TRP_lst:
            if TRP_prob >= random.uniform(0.0, 1.0) and mol[bond[0]].ResName == 'TRP' and mol[bond[0]].Name == 'CA':
                angle, coords_for_rot, initial_coords, angle_TRP = TRP_rot(bond, mol)
            else:
                rotate(mol, coords_for_rot, bond[0], bond[1], angle)
            smth_happen = True
        else:
            smth_happen = False
            if len(TRP_lst) == TRP_count:
                TRP_lst.pop(0)

        if angle == 0 and angle_TRP == 0:
            smth_happen = False

        if smth_happen:

            # evaluate new energy and compare it to old one
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

            if choice and (mol[bond[0]].ResSeq not in TRP_lst) and smth_happen:

                if new_energy <= 400:
                    write_pdb(mol, 'out_300.pdb')
                if mol[bond[0]].ResName == 'TRP':
                    if abs(angle) >= 40 or abs(angle_TRP) >= 40:
                        TRP_lst.append(mol[bond[0]].ResSeq)


                rotations += [(bond[0], bond[1], angle)]
                iterations = 0
                energy_dict[mol[coords_for_rot[0]].ResSeq] = energy_C + energy_LJ/eps
                # use next line if you want to write every step in frame
                #write_pdb(mol, f'for_frames/out_{mol[coords_for_rot[0]].ResSeq}_{step}.pdb')
                logger(f"Step: {step}. Good current rotation:\n"
                       f"First num: {bond[0]}, second num: {bond[1]}, resid: {mol[coords_for_rot[0]].ResSeq},"
                       f" angle: {angle}, trp_angle2: {angle_TRP}, energy: {new_energy}\n\n")
            elif smth_happen:
                iterations += 1
                for i, num in enumerate(coords_for_rot):
                    mol[num].Coordin = initial_coords[i]

                logger(f"Step: {step}. Bad current rotation:\n"
                       f"First num: {bond[0]}, second num: {bond[1]}, resid: {mol[coords_for_rot[0]].ResSeq},"
                       f" angle: {angle}, trp_angle2: {angle_TRP}, energy: {new_energy}\n\n")


            eps = 1
            angle_TRP = 0
            if mol[bond[0]].ResName == 'TRP' and len(TRP_lst) != 0:
                TRP_queue_step += 1
            if TRP_queue_step == TRP_queue_update:
                if len(TRP_lst) != 0:
                    TRP_lst.pop(0)
                TRP_queue_step = 0

        else:
            step -= 1

    write_pdb(mol, 'out.pdb')
    logger(f"The best energy: {best_energy}\n"
           f"Step/Stop {step}/{stop_step}\n")
    logger("\n")
    return rotations, best_energy
