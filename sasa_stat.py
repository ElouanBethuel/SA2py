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


def stat_by_residus(info_pdb, pdb_id):
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
