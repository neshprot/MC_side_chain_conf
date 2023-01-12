from utils import *

def Energy(mol):
    return random.randint(99, 102)


def rotate(mol, coords_for_rot, bond_num1, bond_num2, angle):
    for i in coords_for_rot:
        new_coord = rotation(mol[bond_num1].Coordin, mol[bond_num2].Coordin, mol[i].Coordin, angle)
        mol[i].Coordin = new_coord


def MonteCarlo(mol, graph, rot_bonds, start_energy, attempts):
    rotations = []
    iterations = 0
    count = 0

    while iterations <= attempts and count < 100:
        bond = random.choice(rot_bonds)
        angle = random.uniform(-10, 10)
        coords_for_rot = graph.bfs(bond[1])
        initial_coords = [mol[i].Coordin for i in coords_for_rot]
        rotate(mol, coords_for_rot, bond[0], bond[1], angle)
        energy = Energy(mol)
        count += 1
        if energy > start_energy:
            iterations += 1
            for i, num in enumerate(coords_for_rot):
                mol[num].Coordin = initial_coords[i]
        else:
            rotations += [(bond[0], bond[1], angle)]
            iterations = 0
    return rotations
