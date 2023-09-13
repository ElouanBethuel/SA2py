import random
import copy
import pymol
from pymol import cmd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d

import sasa

RADIUS_H2O = 1.4
NB_POINTS = 20

def create_points_graphic(info_pdb, num_atom, pdb_id):
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

    x = info_pdb[num_atom][4]
    y = info_pdb[num_atom][5]
    z = info_pdb[num_atom][6]

    r = info_pdb[num_atom][6] + RADIUS_H2O
    theta = np.random.uniform(0, np.pi, NB_POINTS)
    phi = np.random.uniform(0, 2*np.pi, NB_POINTS)

    points_x = x + r * np.sin(theta) * np.cos(phi)
    points_y = y + r * np.sin(theta) * np.sin(phi)
    points_z = z + r * np.cos(theta)

    plt.figure()
    plt.title("Viewing points on the solvation sphere")
    axes = plt.axes(projection="3d")
    # displays the central atom in red
    axes.scatter(x, y, z, color='red')
    # displays points in red
    axes.scatter(points_x, points_y, points_z)
    plt.savefig(pdb_id + "_points.png")
    plt.close()

    points = np.array([points_x, points_y, points_z]).T

    return points


def plot_pymol_surface(info, pdb_id):
    """function to generate a protein pymol file.
    This function generates a pymol file based of 
    the PDB file of the protein with a specific 
    coloration of of each atom as a function 
    of their solvent surface area. In red, the atoms 
    most acessible and in blue the less accessible.

    Parameters
    ----------
    info : list of list, 
    pdb_id : string, 

    Returns
    -------
    No returns
    Generate a pymol file (.pse) 
    """

    cmd.load("pdb" + pdb_id + ".ent", "proteine")

    list_color_pdb = []

    for atom in info:
        list_color_pdb.append(atom[8])

    array_color_pdb = np.asarray(list_color_pdb)

    for idx, value in enumerate(array_color_pdb):

        if value <= 1.0:
            color = "blue"
        elif value <= 5.0 :
            color = "green"
        elif value <= 15.0 :
            color = "yellow"
        elif value <= 30.0 :
            color = "orange"
        else:
            color = "red"

        cmd.color(color, f'proteine and id {idx + 1}')

    cmd.show('surface', 'proteine')
    cmd.save(pdb_id + "_surface.pse", "proteine")



def plot_pymol_prot_n(info, num_atom, distance, pdb_id):
    """function to generate a pymol file.
    This function generates a pymol file based of 
    the PDB file of the protein with a specific 
    coloration. In red, the atoms neighboring atoms 
    of the selected atom (num_atom), i.e. those 
    whose distance is less than that of the argument. 
    And in blue the atoms who are not neighbors. 

    Parameters
    ----------
    info : list of list,
    num_atom : int, 
    distance : float,
    pdb_id : string, 

    Returns
    -------
    No returns
    Generate a pymol file (.pse) 
    """

    cmd.load("pdb" + pdb_id + ".ent", "proteine")
    list_color_pdb = []
    list_atoms_neighbors = sasa.neighbors(num_atom, info, distance)

    for id_atom, atom in enumerate(info):
        if atom in list_atoms_neighbors:
            cmd.color("red", f'proteine and id {id_atom + 1}')
        else:
            cmd.color("blue", f'proteine and id {id_atom + 1}')
    
    cmd.save(pdb_id + "_neighbors.pse", "proteine")


def stat_by_atom(info_pdb, pdb_id):
    """Function to generate a .png barplot
    of solvent accessibility by atom category. 
    Make th sum of solvent accessibility area
    of each atom category (C, N, O, S). 

    Parameters
    ----------
    info_pdb : list of list,
    pdb_id : string, 

    Returns
    -------
    dict_atoms : a dictionary, 
    """ 

    dict_atoms = {"C": 0, "N": 0, "O": 0, "S": 0}

    for atom in info_pdb:

        a = atom[1][0]

        if a in dict_atoms:

            dict_atoms[a] = dict_atoms[a] + atom[8]

    list_atoms = list(dict_atoms.keys())
    list_values = list(dict_atoms.values())

    fig = plt.figure()
    plt.title("Barplot solvent accessible surface area by atoms")
    width = 0.5
    plt.bar(list_atoms, list_values, width, color='b')
    plt.savefig(pdb_id + "_barplot_atoms.png")
    plt.close()

    return dict_atoms



def stat_residus(info_pdb, pdb_id):
    """Function to generate a .png barplot
    of solvent accessibility by residu category. 
    Make th sum of solvent accessibility area
    of each residu (amino acide) category. 

    Parameters
    ----------
    info_pdb : list of list,
    pdb_id : string, 

    Returns
    -------
    dict_aa : a dictionary, 
    """
    acc_by_amino_acide = {}
    
    for atom in info_pdb:
    
        num_aa = atom[3]
        acc = atom[8]
        
        if num_aa in acc_by_amino_acide:
            acc_by_amino_acide[num_aa] += acc

        else:
            acc_by_amino_acide[num_aa] = acc
    
    num_aa = list(acc_by_amino_acide.keys())
    acc = list(acc_by_amino_acide.values())
            
    fig = plt.figure()
    plt.title("solvent accessible surface area per residue")
    plt.xticks(np.arange(min(num_aa), max(num_aa)+1, 10.0))
    plt.plot(num_aa, acc, marker='o', linestyle='-')
    plt.savefig(pdb_id + "_sasa_amino_acide.png")
    plt.close()


