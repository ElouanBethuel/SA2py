
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
"""


__authors__ = ("Elouan Bethuel", "M2BI Université Paris Cité")
__contact__ = ("elouan.bethuel@etu-u.paris.fr")
__date__ = "05-09-2023"


import os
import sys
import time
import Bio
from Bio.PDB import PDBList

import sasa
import sasa_stat


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
    info = sasa.parse_pdb_file("pdb" + pdb_id + ".ent")
    list_all_atoms_points = sasa.all_atoms_points(info)
    list_all_neighbors = sasa.all_neighbors(info, 15)
    sasa.access_solvant(info, list_all_neighbors, list_all_atoms_points, pdb_id)

    # Plots
    sasa_stat.create_points_graphic(info, 1, pdb_id)
    sasa_stat.plot_pymol_surface(info, pdb_id)
    sasa_stat.plot_pymol_prot_n(info, 1, 10, pdb_id)

    # Statistics
    sasa_stat.stat_by_atom(info, pdb_id)
    sasa_stat.stat_by_residus(info, pdb_id)

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

