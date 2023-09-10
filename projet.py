import sys
import random
import copy
import numpy as np
import Bio
from Bio.PDB import PDBList
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d

RADIUS_H2O = 1.4
NB_POINTS = 100


"""
Script for calculating a protein's solvent accessibility 
air from information contained in its PDB file.  
Generates: 
- a text file with the solvent accessibility of each atom, 
- two graphs to visualize the points created around the atoms and one 
  to visualize the selection of neighboring atoms to an atom, 
- two graphs with statistics on solvent accessibility by atom 
  category and by amino acid category. 
"""



def parse_pdb_file(path_pdb_file):
    """Function to parse PDB files.
    
    For each ATOM lines extract :
        - residu number
        - 
    
    And add Van der Wadd Radius
    For obtain : 
    ['411', 'CB', 'ALA', -27.302, 9.269, -5.772, 1.7]
    
    Parameters
    ----------
    path_pdb_file : str
        the path to the PDB file 
    

    Returns
    -------
    int
        Le produit des deux nombres.
    """

    info = []

    with open(path_pdb_file, "r") as pdb_file:

        for line in pdb_file:
        
            if line.startswith("ATOM"):

                splt_line = line.split()

                if splt_line[1] != "H":

                    new_line = [splt_line[1], splt_line[2], splt_line[3],
                                float(splt_line[6]), float(splt_line[7]), float(splt_line[8])]

                    # add Van der Wadd radius

                    if new_line[1].find("N") != -1:
                        new_line.append(1.55)

                    if new_line[1].find("O") != -1:
                        new_line.append(1.52)

                    if new_line[1].find("S") != -1:
                        new_line.append(1.8)

                    if new_line[1].find("C") != -1:
                        new_line.append(1.7)

                    info.append(new_line)

    return info


# function to generated points

def create_points(x, y, z, radius):

    r = radius + RADIUS_H2O
    theta = np.random.uniform(0, np.pi, NB_POINTS)
    phi = np.random.uniform(0, 2*np.pi, NB_POINTS)

    points_x = x + r * np.sin(theta) * np.cos(phi)
    points_y = y + r * np.sin(theta) * np.sin(phi)
    points_z = z + r * np.cos(theta)

    points = np.array([points_x, points_y, points_z]).T

    return points


def create_points_graphic(info_pdb, num_atom):

    x = info_pdb[num_atom][3]
    y = info_pdb[num_atom][4]
    z = info_pdb[num_atom][5]

    r = info_pdb[num_atom][6] + RADIUS_H2O
    theta = np.random.uniform(0, np.pi, NB_POINTS)
    phi = np.random.uniform(0, 2*np.pi, NB_POINTS)

    points_x = x + r * np.sin(theta) * np.cos(phi)
    points_y = y + r * np.sin(theta) * np.sin(phi)
    points_z = z + r * np.cos(theta)

    plt.figure()
    axes = plt.axes(projection="3d")
    axes.scatter(x, y, z, color='red') # displays the central atom in red
    axes.scatter(points_x, points_y, points_z) # displays points in red
    plt.show()
    
    points = np.array([points_x, points_y, points_z]).T

    return points


def all_atoms_points(info_pdb):

    list_all_points = {}

    for atom in info_pdb:

        points = create_points(atom[3], atom[4], atom[5], atom[6])
        list_all_points[atom[0]] = points

    return list_all_points


# Retourne la liste des voisins d'un atome à partir de ses coordonnées
# et d'un  seuil (treeshold : distance en Angstrom)

def neighbors(num_atom, pdb, threeshold):

    info_pdb = copy.deepcopy(pdb)

    x = info_pdb[num_atom-1][3]
    y = info_pdb[num_atom-1][4]
    z = info_pdb[num_atom-1][5]

    del (info_pdb[num_atom-1])
    list_neighbors = []

    for atom in info_pdb:

        x_nb = atom[3]
        y_nb = atom[4]
        z_nb = atom[5]

        d = np.sqrt((x - x_nb)**2 + (y - y_nb)**2 + (z - z_nb)**2)

        if d < threeshold:
            list_neighbors.append(atom)

    return list_neighbors


def all_neighbors(info_pdb, threeshold):

    list_all_neighbors = []

    for i in range(len(info_pdb)):
        list_all_neighbors.append(neighbors(i+1, info, threeshold))

    return list_all_neighbors


def plot_protein(info, num_atom, distance):

    list_x_pdb = []
    list_y_pdb = []
    list_z_pdb = []

    for atom in info:

        list_x_pdb.append(atom[3])
        list_y_pdb.append(atom[4])
        list_z_pdb.append(atom[5])

    array_x_pdb = np.asarray(list_x_pdb)
    array_y_pdb = np.asarray(list_y_pdb)
    array_z_pdb = np.asarray(list_z_pdb)

    list_atoms_neighbors = neighbors(num_atom, info, distance)

    list_x_neighbors = []
    list_y_neighbors = []
    list_z_neighbors = []

    for atom in list_atoms_neighbors:

        list_x_neighbors.append(atom[3])
        list_y_neighbors.append(atom[4])
        list_z_neighbors.append(atom[5])

    array_x_neighbors = np.asarray(list_x_neighbors)
    array_y_neighbors = np.asarray(list_y_neighbors)
    array_z_neighbors = np.asarray(list_z_neighbors)

    plt.figure()
    axes = plt.axes(projection="3d")
    axes.scatter(array_x_pdb, array_y_pdb, array_z_pdb, color="blue")
    axes.scatter(info[num_atom][3], info[num_atom][4], info[num_atom][5], s=300, color="green")
    axes.scatter(array_x_neighbors, array_y_neighbors, array_z_neighbors, s=200, color='red')
    plt.show()


# fonction pour déterminer le pourcentage de sondes accessibles au solvant
# pour chaque atome de la protéine

def access_solvant(info, list_all_neighbors, list_all_atoms_points):

    sum_acc_area = 0
    i = 0

    for atom in info:

        sondes = list_all_atoms_points[atom[0]]
        sondes_inacessibles = copy.deepcopy(sondes)
        neighbors = list_all_neighbors[i]

        i += 1
        nb_sondes = 0

        for neighbor in neighbors:

            x_n = neighbor[3]
            y_n = neighbor[4]
            z_n = neighbor[5]

            point_n = np.array([x_n, y_n, z_n])

            # calcule les distances euclidiennes entre l'atome voisin et les sondes de l'atome
            distances = np.linalg.norm(point_n - sondes, axis=1)

            points_inacessibles = distances <= 1.4 + atom[6]
            indices_points_inacessibles = np.where(points_inacessibles)[0]
            sondes_inacessibles[indices_points_inacessibles] = -1

        # divise par 3 car 3 coordonnées x,y,z  pour chaque sonde
        nb_sondes_acessibles = (np.size(sondes)/3) - np.count_nonzero(sondes_inacessibles == -1)/3
        prct_acc_atome = ((nb_sondes_acessibles / (np.size(sondes)/3)) * 100)

        area_acc_atome = (prct_acc_atome/100) * 4 * np.pi*((1.4 + atom[6])**2)
        sum_acc_area = sum_acc_area + area_acc_atome

        print(f"{atom} \t")
        print(f"Number of points accessible to the solvent : {nb_sondes_acessibles}% \t")
        print(f"Atom accessibility percentage : {round(prct_acc_atome,2)}% \t")
        print(f"Solvent accessibility area of ​​the atom : {round(area_acc_atome,2)} Å² \n")

        # add solvent accessibility area of ​​each atom for statistics 
        atom.append(area_acc_atome)

    sum_acc_area = round(sum_acc_area, 2)
    print(f"Final result :\t")
    print(f"Total area exposed to solvent : {sum_acc_area} Å² \n")


def stat_by_atom(info_pdb):

    dict_atoms = {"C": 0, "N": 0, "O": 0, "S": 0}

    for atom in info_pdb:

        a = atom[1][0]

        if a in dict_atoms:

            dict_atoms[a] = dict_atoms[a] + atom[7]

    list_atoms = list(dict_atoms.keys())
    list_values = list(dict_atoms.values())

    fig = plt.figure()
    width = 0.5
    plt.bar(list_atoms, list_values, width, color='b')
    # plt.savefig('Barplot_solvant_access_area_by_atom.png')
    plt.show()

    return dict_atoms


def stat_by_residus(info_pdb):

    dict_aa = {"ALA": 0, "ARG": 0, "ASN": 0, "ASP": 0,
               "CYS": 0, "GLY": 0, "HIS": 0, "ILE": 0,
               "LEU": 0, "LYS": 0, "MET": 0, "PHE": 0,
               "PRO": 0, "SER": 0, "THR": 0, "TRP": 0,
               "TYR": 0, "VAL": 0}

    for atom in info_pdb:

        if atom[2] in dict_aa:

            dict_aa[atom[2]] = round(dict_aa[atom[2]] + atom[7], 2)

    list_keys = list(dict_aa.keys())
    list_values = list(dict_aa.values())

    fig = plt.figure()
    width = 0.5
    plt.bar(list_keys, list_values, width, color='b')
    # plt.savefig('Barplot_solvant_access_area_by_atom.png')
    plt.show()

    return dict_aa


########## Main ##############

# Control user arguments 

if len(sys.argv) != 2:
    sys.exit("ERREUR : il faut exactement un argument.")
    
pdb_id = sys.argv[1]

# Download PDB file 
pdbl = PDBList()
pdbl.retrieve_pdb_file(pdb_id, pdir='./', file_format='pdb')

# Calculs
info = parse_pdb_file(f"pdb{pdb_id}.ent")
list_all_atoms_points = all_atoms_points(info)
list_all_neighbors = all_neighbors(info, 15)
access_solvant(info, list_all_neighbors, list_all_atoms_points)

# Plots
create_points_graphic(info, 1)
plot_protein(info, 1, 10)

# Statistics
stat_by_atom(info)
stat_by_residus(info)

