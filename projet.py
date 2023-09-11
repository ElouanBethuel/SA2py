
"""
Script for compute the solvent accessible surface
area from information contained in Protein Data Bank PDB file.

outputs:

- The solvent accessible surface area in Å².
- A text file with the solvent accessibility of each atom.
- A graph to visualize the points created around a atom.
- A pymol file showing the selection of neighboring atoms.
- A graph showing solvent accessibility by atom category.
- A graph showing solvent accessibility by amino acid category.
- A pymol file to visualize the solvent accessibility. 
    
Usage:
======

    python projet.py pdb_id
    
    argument1:  PDB ID of the protein studied
    arugment2:
"""


__authors__ = ("Elouan Bethuel", "M2BI Université Paris Cité")
__contact__ = ("elouan.bethuel@etu-u.paris.fr")
__date__ = "05-09-2023"


import os
import sys
import time
import random
import copy
import pymol
from pymol import cmd
import numpy as np
import Bio
from Bio.PDB import PDBList
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d


RADIUS_H2O = 1.4
NB_POINTS = 100


def parse_pdb_file(path_pdb_file):
    """Function to parse PDB files.
    
    For each ATOM lines extract:
        - Atom serial number
        - Atom name
        - Residue name
        - Cartesian coordinates
    
    And add Van der Wadd Radius
    Obtained a new line like this: 
     
    ['411', 'CB', 'ALA', -27.302, 9.269, -5.772, 1.7]
    
    Parameters
    ----------
    path_pdb_file : str, the path to the PDB file 
    
    Returns
    -------
    info: list of list
    """

    info = []

    with open(path_pdb_file, "r") as pdb_file:

        for line in pdb_file:
        
            if line.startswith("ATOM"):

                splt_line = line.split()

                if splt_line[1] != "H":

                    new_line = [splt_line[1], splt_line[2], 
                               splt_line[3], float(splt_line[6]),
                               float(splt_line[7]), float(splt_line[8])]

                    # Add Van der Waals radius

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



def create_points(x, y, z, radius):
    """Function to generated points around an atom.
    
    This function generates a set of points to form 
    a sphere around an atom. The sphere corresponds 
    to the atom's solvation sphere. Its radius is equal 
    to the sum of the atom's Van der Waals radius 
    and the radius of a water molecule (1.4 Ångström). 
    
    Parameters
    ----------
    x, y, z : int, the atom's cartesian coordinates
    radius : float, the atom's Van der Waals radius
    
    Returns
    -------
    points: a list of list
    """

    r = radius + RADIUS_H2O
    theta = np.random.uniform(0, np.pi, NB_POINTS)
    phi = np.random.uniform(0, 2*np.pi, NB_POINTS)

    points_x = x + r * np.sin(theta) * np.cos(phi)
    points_y = y + r * np.sin(theta) * np.sin(phi)
    points_z = z + r * np.cos(theta)

    points = np.array([points_x, points_y, points_z]).T

    return points


def create_points_graphic(info_pdb, num_atom):
    """Function to generated points around an atom.

    This function generates a set of points 
    to form a sphere around an atom. 
    The sphere corresponds to the atom's 
    solvation sphere. Its radius is equal 
    to the sum of the atom's Van der Waals radius 
    and the radius of a water molecule (1.4 Ångström).
    
    Parameters
    ----------
    x, y, z : int, the atom's cartesian coordinates
    radius : float, the atom's Van der Waals radius

    Returns
    -------
    points: a list of list, list of cartesian coordinates of each point
    """

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
    plt.title("Viewing points on the solvation sphere")
    axes = plt.axes(projection="3d")
    axes.scatter(x, y, z, color='red') # displays the central atom in red
    axes.scatter(points_x, points_y, points_z) # displays points in red
    plt.savefig(pdb_id + "_points.png")
    plt.close()
    
    points = np.array([points_x, points_y, points_z]).T

    return points



def all_atoms_points(info_pdb):
    """Functions to generate points for each atom.

    This function uses the 'create-points' function 
    to generate points for all atoms contained in the PDB file. 

    Parameters
    ----------
    info_pdb : list of list, 
    
    Returns
    -------
    list_all_points: 
    """

    list_all_points = {}

    for atom in info_pdb:

        points = create_points(atom[3], atom[4], atom[5], atom[6])
        list_all_points[atom[0]] = points

    return list_all_points



def neighbors(num_atom, pdb, threshold):
    """Function to identify an atom's neighbors.
    
    This function identifies all neighbors 
    located at a certain distance from the atom.
    Using the atom's Cartesian coordinates, 
    calculate the distance between this atom 
    and all other atoms in the PDB file. 
    If the distance is less than the threshold, 
    the atom is considered as a neighbor. 

    Parameters
    ----------
    num_atom : int,
    pdb : list of list, 
    threshold : int, 

    Returns
    -------
    list_neighbors: 
    """

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

        if d < threshold:
            list_neighbors.append(atom)

    return list_neighbors


def all_neighbors(info_pdb, threshold):
    """Function to identify neighbors of all atoms. 
    
    This function uses the "neighbors" function 
    to identify the neighbors of each atom in the PDB file.

    Parameters
    ----------
    info_pdb : list of list, 
    threshold : int, 

    Returns
    -------
    list_all_neighbors: list of list, 
    """

    list_all_neighbors = []

    for i in range(len(info_pdb)):
        list_all_neighbors.append(neighbors(i+1, info, threshold))

    return list_all_neighbors



def plot_pymol_surface(info):
	
	cmd.load("pdb" + pdb_id + ".ent", "proteine")

	list_color_pdb = []
    
	for atom in info:
		list_color_pdb.append(atom[7])

	array_color_pdb = np.asarray(list_color_pdb)

	def colormap(valeur):
		if valeur <= 1.0 :
		    return 'blue'
		elif valeur <= 5.0 :
		    return 'green'
		elif valeur <= 15.0 :
		    return 'yellow'
		elif valeur <= 30.0 :
		    return 'orange'
		else:
		    return 'red'

	for idx, value in enumerate(array_color_pdb):
		
		color = colormap(value)
		cmd.color(color, f'proteine and id {idx + 1}')
	
	cmd.show('surface', 'proteine')
	cmd.save(pdb_id + "_surface.pse", "proteine")
	

def plot_pymol_prot_n(info, num_atom, distance):
	
	cmd.load("pdb" + pdb_id + ".ent", "proteine")
	list_color_pdb = []
	list_atoms_neighbors = neighbors(num_atom, info, distance)
	
	for id_atom, atom in enumerate(info):
		if atom in list_atoms_neighbors:
			cmd.color("red", f'proteine and id {id_atom + 1}')
		else:
			cmd.color("blue", f'proteine and id {id_atom + 1}')
			
	#cmd.show('sphere', 'proteine')
	cmd.save(pdb_id + "_neighbors.pse", "proteine")


def access_solvant(info, list_all_neighbors, list_all_atoms_points):
    """Fonction pour déterminer le pourcentage de sondes accessibles au solvant
       pour chaque atome de la protéine
    """
    
    file_name = pdb_id + "_ouput.txt"
    with open(file_name, "w") as filout:
    
        prot_acc_area = 0
        i = 0 

        for atom in info:

            points = list_all_atoms_points[atom[0]]
            occulted_pts = copy.deepcopy(points)
            neighbors = list_all_neighbors[i]
            i = i+1

            for neighbor in neighbors:

                x_n = neighbor[3]
                y_n = neighbor[4]
                z_n = neighbor[5]

                point_n = np.array([x_n, y_n, z_n])

                # calculates Euclidean distances between the neighboring atom and the atom's points
                distances = np.linalg.norm(point_n - points, axis=1)

                occulted = distances <= RADIUS_H2O + atom[6]
                indices_occulted = np.where(occulted)[0]
                occulted_pts[indices_occulted] = -1

            # divise par 3 car 3 coordonnées x,y,z  pour chaque sonde
            accessible_pts = (np.size(points)/3) - np.count_nonzero(occulted_pts == -1)/3
            prct_acc_atom = ((accessible_pts / (np.size(points)/3)) * NB_POINTS)

            acc_area_atom = (prct_acc_atom/NB_POINTS) * 4 * np.pi*((RADIUS_H2O + atom[6])**2)
            prot_acc_area = prot_acc_area + acc_area_atom

            filout.write(f"{atom}\n")
            filout.write(f"Number of points accessible to the solvent : {accessible_pts}%\n")
            filout.write(f"Atom accessibility percentage : {round(prct_acc_atom,2)}%\n")
            filout.write(f"Solvent accessibility area of ​​the atom : {round(acc_area_atom,2)} Å²\n\n")

            # add solvent accessibility area of ​​each atom to pdb info
            atom.append(acc_area_atom)

        prot_acc_area = round(prot_acc_area, 2)
        
        filout.write("\n")
        filout.write("=================\n")
        filout.write(f"Final result :\n")
        filout.write(f"Total surface of protein exposed to solvent : {prot_acc_area} Å² \n")
        
        print("\n")
        print("=================")
        print(f"Final result :\t")
        print(f"Total surface of protein exposed to solvent : {prot_acc_area} Å² \n")


def stat_by_atom(info_pdb):

    dict_atoms = {"C": 0, "N": 0, "O": 0, "S": 0}

    for atom in info_pdb:

        a = atom[1][0]

        if a in dict_atoms:

            dict_atoms[a] = dict_atoms[a] + atom[7]

    list_atoms = list(dict_atoms.keys())
    list_values = list(dict_atoms.values())

    fig = plt.figure()
    plt.title("Barplot solvent accessible surface area by atoms")
    width = 0.5
    plt.bar(list_atoms, list_values, width, color='b')
    plt.savefig(pdb_id + "_barplot_atoms.png")
    plt.close()

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
    plt.title("Barplot solvent accessible surface area by residue")
    width = 0.5
    plt.bar(list_keys, list_values, width, color='b')
    plt.savefig(pdb_id + "_barplot_aa.png")
    plt.close()

    return dict_aa


########## Main ##############

if __name__ == "__main__":

    # Checking user arguments
    if len(sys.argv) != 2:
        sys.exit("ERROR: exactly one argument is required")
    pdb_id = sys.argv[1]
    if len(pdb_id) < 4:
        sys.exit("ERROR: PDB ID must have at least 4 letters")
    print("\n")
    print(f"Calculates the solvent accessible surface area from the {pdb_id} PDB file :\n")
    time.sleep(2)

    # Download PDB file
    print(f"Downloading the PDB file\t")
    time.sleep(1.2)
    print(f"Loading the PDB file\n")
    time.sleep(1.2)
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, pdir='./', file_format='pdb')

    # Compute the solvent accessible surface area
    print("\n")
    print("Calculation ongoing, please wait a few seconds...\n")
    info = parse_pdb_file(f"pdb{pdb_id}.ent")
    list_all_atoms_points = all_atoms_points(info)
    list_all_neighbors = all_neighbors(info, 15)
    access_solvant(info, list_all_neighbors, list_all_atoms_points)

    # Plots
    create_points_graphic(info, 1)
    plot_pymol_surface(info)
    plot_pymol_prot_n(info, 1, 10)

    # Statistics
    stat_by_atom(info)
    stat_by_residus(info)

    # Create and delete folders
    if not os.path.exists("pdb_folder"):
        os.mkdir("pdb_folder")
    if not os.path.exists("outputs"):
        os.mkdir("outputs")
    if os.path.exists("obsolete"):
        os.rmdir("obsolete")

    # Move files
    path = os.getcwd()
    os.rename(path + "/pdb" + pdb_id + ".ent", path + "/pdb_folder/pdb" + pdb_id + ".ent")
    os.rename(path + "/" + pdb_id + "_ouput.txt", path + "/outputs/" + pdb_id + "_ouput.txt")
    os.rename(path + "/" + pdb_id + "_barplot_aa.png", path + "/outputs/" + pdb_id + "_barplot_aa.png")
    os.rename(path + "/" + pdb_id + "_barplot_atoms.png", path + "/outputs/" + pdb_id + "_barplot_atoms.png")
    os.rename(path + "/" + pdb_id + "_points.png", path + "/outputs/" + pdb_id + "_points.png")
    os.rename(path + "/" + pdb_id + "_surface.pse", path + "/outputs/" + pdb_id + "_surface.pse")
    os.rename(path + "/" + pdb_id + "_neighbors.pse", path + "/outputs/" + pdb_id + "_neighbors.pse")

