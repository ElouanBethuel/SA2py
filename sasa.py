import random
import copy
import pymol
from pymol import cmd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d

RADIUS_H2O = 1.4
NB_POINTS = 20


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
                                splt_line[3], float(splt_line[5]),
                                float(splt_line[6]),
                                float(splt_line[7]),
                                float(splt_line[8])]

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

        points = create_points(atom[4], atom[5], atom[6], atom[7])
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

    x = info_pdb[num_atom-1][4]
    y = info_pdb[num_atom-1][5]
    z = info_pdb[num_atom-1][6]

    del (info_pdb[num_atom-1])
    list_neighbors = []

    for atom in info_pdb:

        x_nb = atom[4]
        y_nb = atom[5]
        z_nb = atom[6]

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
        list_all_neighbors.append(neighbors(i+1, info_pdb, threshold))

    return list_all_neighbors


def access_solvant(info, list_all_neighbors, list_all_atoms_points, pdb_id):
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
    info : list of list,
    list_all_neighbors : list of list, 
    list_all_atoms_points : list of list, 
    pdb_id : string, 
    
    Returns
    -------
    No returns 
    Print the total surface of protein exposed to solvent
    And writh in a "pdb_id_output.txt" file the accessibilit by atom. 
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

                x_n = neighbor[4]
                y_n = neighbor[5]
                z_n = neighbor[6]

                point_n = np.array([x_n, y_n, z_n])

                # calculates Euclidean distances between the neighboring atom and the atom's points
                distances = np.linalg.norm(point_n - points, axis=1)

                occulted = distances <= RADIUS_H2O + atom[7]
                indices_occulted = np.where(occulted)[0]
                occulted_pts[indices_occulted] = -1

            # divise par 3 car 3 coordonnées x,y,z  pour chaque sonde
            accessible_pts = (np.size(points)/3) - np.count_nonzero(occulted_pts == -1)/3
            prct_acc_atom = ((accessible_pts / (np.size(points)/3)) * NB_POINTS)

            acc_area_atom = (prct_acc_atom/NB_POINTS) * 4 * np.pi*((RADIUS_H2O + atom[7])**2)
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


