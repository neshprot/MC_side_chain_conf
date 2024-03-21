"""
Main script for MC sim
"""
import time

from utils import math, rotate, write_pdb, TorsionAngle, get_torsions, amino_acid,\
    read_inp, read_pdb, read_psf, rotate_trp_tors_angle
from logger import FileLogger
import random
import numpy as np
import os
import configparser
from graph import Graph
import subprocess
import numpy




def montecarlo(mol, graph, rot_bonds, rotating_resid, trp_resid, config):
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
    step = 1  # const
    energy_dict = {}
    res_dict = {}
    k = 1.987204259e-3  # Boltzmann constant
    t1 = float(config['PARAMS']['first_temperature'])  # first temperature
    t2 = float(config['PARAMS']['second_temperature'])  # second temperature 1000
    count_mk =  int(config['PARAMS']['count_iteretions_before_md'])  
    t = t1
    eps_amen = 0.4  # amendment of epsilon 0.4
    eps = 1  # standard amendment of epsilon
    best_energy = 'None'
    T_count = 0  # const
    T_step = 0  # const
    T_lim = 100  # every T_count steps temperature changes
    eps_step = 0  # const
    eps_lim = 100  # every eps_lim steps epsilon changes
    probability_for_rotamer = 0.2 #probability for use rotamer map

    files_dir = config['DIRS']['files_dir'] #all files dir
    rotamer_dir = config['DIRS']['rotamers_dir'] #rotamers dir

    md_start = config['COMANDS']['md_start'] #Namd comand for Molecular Dinamic
    new_pdb_after_md = config['COMANDS']['new_pdb_after_md'] #VMD comand for create pdb after MD

    attempts = float(config['PARAMS']['attempts'])
    stop_step = float(config['PARAMS']['stop_step'])

    ILE_prob = 1

    TRP_queue_update = 50
    TRP_queue_step = 0  # const
    angle_TRP = 0  # const
    TRP_lst = []  # const

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
        # print(f'{mol[coords_for_rot[0]].resseq} - {energy_C+energy_LJ}')
        return energy_C, eps_amen * energy_LJ


    # evaluate start energy
    def start_energy():
        summ = 0
        TRP_count = 0
        ini_rot_bonds = rot_bonds
        for atom in mol:
            if mol[atom].resseq not in rotating_resid:
                continue
            res = []
            for x in ini_rot_bonds.values():
                res.extend(x if isinstance(x, list) else [x])
            if any(atom in i for i in res) and mol[atom].name == 'CA':
                atom_lst = graph.bfs(atom)
                energy_C, energy_LJ = Energy(mol, atom_lst, eps_amen=1)
                res_energy = energy_C + energy_LJ
                # record the energy of each aminoacid
                energy_dict.update({mol[atom].resseq: res_energy})
                # record all atoms in a particular amino acid
                res_dict.update({mol[atom].resseq: atom_lst})
                summ += res_energy
                all_nums.append(atom)

    # special function for ILE
    def ILE_rot(bond, mol):

        ile_angle1 = random.choice([30, 60, 90, 180]) * random.choice([-1, 1])
        ile_angle2 = random.choice([30, 60, 90, 180]) * random.choice([-1, 1])

        coords_for_rot1 = graph.bfs(bond[1])
        coords_for_rot2 = graph.bfs(bond[1] + 6)

        rotate(mol, coords_for_rot1, bond[0], bond[1], ile_angle1)
        rotate(mol, coords_for_rot2, bond[1], bond[1] + 6, ile_angle2)
        return [ile_angle1, ile_angle2], coords_for_rot, initial_coords
    

    def rotamer_open(rotamer_map, file_path):
        #special funtion for open rotamers map

        with open(file_path, 'r') as file:
            next(file)  # Пропускаем заголовок
            for line in file:
                data = line.strip().split(',')
                amino = data[0]
                rotamer = int(data[1])
                probability = float(data[2])
                chi_values = [(float(val) + 360) % 360 if val != "NA" else None for val in data[3:]]

                if amino in rotamer_map:
                    rotamer_map[amino][rotamer] = {
                        "Probability": probability,
                        "Chi_Values": chi_values
                    }
                else:
                    rotamer_map[amino] = {
                        rotamer: {
                            "Probability": probability,
                            "Chi_Values": chi_values
                        }
                    }
        return rotamer_map

    start_energy()

    start_was = True
    list_residue_classes = {}
    for res in rotating_resid:
        residue = TorsionAngle(res)
        list_residue_classes[res] = residue
    # main block with consecutive rotations
    

    def transition_probability(new_energy, old_energy):
        if new_energy <= old_energy:
            return True
        else:
            prob_up = math.exp(-(new_energy - old_energy) / (k * t))
            if numpy.random.random()<= prob_up:
                return True
            else:
                return False
            
    def md():
        # Список файлов для удаления
            files_to_delete = [
                "d95n_0_working_min-0.coor",
                "d95n_0_working_min-0.dcd",
                "d95n_0_working_min-0.restart.coor",
                "d95n_0_working_min-0.restart.coor.old",
                "d95n_0_working_min-0.restart.vel",
                "d95n_0_working_min-0.restart.vel.old",
                "d95n_0_working_min-0.restart.xsc",
                "d95n_0_working_min-0.restart.xsc.old",
                "d95n_0_working_min-0.vel",
                "d95n_0_working_min-0.xsc",
                "d95n_0_working_min-0.xst",
                "FFTW_NAMD_2.14_Linux-x86_64-multicore.txt",
                "m1.log"
            ]

            # Путь к директории с файлами
            input_files_dir = files_dir

            # Смена директории
            os.chdir(input_files_dir)
            write_pdb(mol, 'out.pdb')
            time.sleep(5)

            # Вывод содержимого директории

            # Запуск NAMD и запись лога
            os.system(md_start)

            time.sleep(10)

            # Запуск VMD для обработки данных
            os.system(new_pdb_after_md)

            time.sleep(10)

            # Удаление файлов
            for file in files_to_delete:
                try:
                    os.remove(os.path.join(input_files_dir, file))
                    print(f"{file} удален успешно")
                except FileNotFoundError:
                    print(f"{file} не найден")
                except Exception as e:
                    print(f"Ошибка при удалении {file}: {e}")

            # Возврат к исходной директории
            os.chdir("..")
            time.sleep(5)

    while iterations <= attempts and step < stop_step:
        if step % count_mk == 0:
            md()

            # config constants
            pdb_file = config['PDB']['AFTER_MD']
            inp_file = config['INP']['File']
            psf_file = config['PSF']['File']
            pdb_wild_file = config['PDB']['Wild_file']
            

            const_dict = read_inp(inp_file)
            mol = read_pdb(pdb_file, const_dict)
            wild_mol = read_pdb(pdb_wild_file, const_dict)
            mol = rotate_trp_tors_angle(wild_mol, mol, trp_resid)
            bonds, rot_bonds = amino_acid(mol, rotating_resid)
            
            read_psf(psf_file, mol)
            energy_dict = {}
            res_dict = {}
            charge_dict = dict()

            all_nums = []

            start_was = False
            
            graph = Graph(bonds)

            start_energy()
            start_was = True

        random_residue = random.choice(rotating_resid)
        rot = rot_bonds[random_residue]
        if not rot:
            continue
        rotamer_map = {}
        rotamer_map = rotamer_open(rotamer_map, rotamer_dir)

        if random.random() < probability_for_rotamer and mol[rot[0][0]].resname in rotamer_map.keys() and mol[
            rot[0][0]].resname != 'TRP':
            #adjusting the rotameter map for rotation
            list_residue_classes[random_residue].update_rotamer(rotamer_map[mol[rot[0][0]].resname])
            dmap = list_residue_classes[random_residue].rotamers()
            rotamers = list(dmap.keys())
            probabilities = [dmap[rotamer]["Probability"] for rotamer in rotamers]
            chosen_rotamer = random.choices(rotamers, weights=probabilities, k=1)[0]
            torsions = dmap[chosen_rotamer]['Chi_Values']
            bonds, _ = amino_acid(mol, [random_residue])
            my_tors = get_torsions(mol, bonds)
            
            result = [i for i in my_tors if (i[1], i[2]) in rot and i[0]]
            torsions = [torsions[a] for a in range(len(torsions)) if torsions[a] != None and len(result) >= a + 1]
            result = result[0:len(torsions)]
            dif = []
            for i in range(len(torsions)):
                dif += [result[i][4] - torsions[i]]
            rot_dict = {}
            for torsion in range(len(dif)):
                bond = (result[torsion][1], result[torsion][2])
                rot_dict[bond] = (-1) * np.sign(dif[torsion]) * abs(dif[torsion])

        else:
            bond = random.choice(rot)  # selecting a random bond for rotation
            first_prob = 0.9
            second_prob = 0.1
            angle1 = random.randint(0, 10)
            angle2 = random.randint(90, 180)
            angle_sign = random.choice([-1, 1])
            angle = random.choices([angle1, angle2], weights=[first_prob, second_prob])[0] * angle_sign
            rot_dict = {bond: angle}

        for bond, angle in rot_dict.items():
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
            list_residue_classes[random_residue].update_torsion_angles(bond, angle)
            rotate(mol, coords_for_rot, bond[0], bond[1], angle)
            # evaluate new energy and compare it to old one
            energy_C, energy_LJ = Energy(mol, res_dict[mol[coords_for_rot[0]].resseq], eps)
            new_energy = energy_C + energy_LJ
            old_energy = energy_dict[mol[coords_for_rot[0]].resseq]
            # whether the system has moved to a higher energy step
            choice = transition_probability(new_energy, old_energy)
            if choice:
                rotations += [(bond[0], bond[1], angle)]
                iterations = 0
                energy_dict[mol[coords_for_rot[0]].resseq] = energy_C + energy_LJ / eps
                # use next line if you want to write every step in frame
                # write_pdb(mol, f'for_frames/out_{mol[coords_for_rot[0]].resseq}_{step}.pdb')
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
