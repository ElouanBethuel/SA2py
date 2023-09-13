
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
    
    argument1: PDB ID of the protein studied
    argument2: Number of points by atom to model the solvation sphere
"""

__authors__ = ("Elouan Bethuel", "M2BI Université Paris Cité")
__contact__ = ("elouan.bethuel@etu-u.paris.fr")
__date__ = "05-09-2023"

import os
import sys
import time
import argparse
import Bio
from Bio.PDB import PDBList

import sasa
import sasa_stat


########## Main ##############

if __name__ == "__main__":

    # Checking user arguments
    parser = argparse.ArgumentParser(description=("Program to calculate the Solvent Accessible Surface Area (SASA) of a protein"))
    parser.usage = "main.py pdb_id number_points -s"
    parser.add_argument("pdb_id", type=str, help="ID of the PDB file")
    parser.add_argument("number_points", type=int, help="Number of points by atom to model the solvation sphere")
    parser.add_argument("-s", "--statistical", action="store_true", help=("produce statisticals graphics"))
    args = parser.parse_args()
    print("\n")
    print(f"Calculates the solvent accessible surface area from the {args.pdb_id} PDB file :\n")
    time.sleep(1)

    # Download PDB file
    print(f"Downloading the PDB file\t")
    time.sleep(1)
    print(f"Loading the PDB file\n")
    time.sleep(1)
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(args.pdb_id, pdir='./', file_format='pdb')

    # Compute the solvent accessible surface area
    print("\n")
    print("Calculation ongoing, please wait a few seconds...\n")
    info = sasa.parse_pdb_file("pdb" + args.pdb_id + ".ent")
    list_all_atoms_points = sasa.all_atoms_points(info)
    list_all_neighbors = sasa.all_neighbors(info, 15)
    sasa.access_solvant(info, list_all_neighbors, list_all_atoms_points, args.pdb_id)

    # Statistics
    if args.statistical: 
        sasa_stat.stat_by_atom(info, args.pdb_id)
        sasa_stat.create_points_graphic(info, 1, args.pdb_id)
        sasa_stat.plot_pymol_surface(info, args.pdb_id)
        sasa_stat.plot_pymol_prot_n(info, 1, 10, args.pdb_id)
        sasa_stat.stat_residus(info, args.pdb_id)

    # Create and delete folders
    if not os.path.exists("pdb_folder"):
        os.mkdir("pdb_folder")
    if not os.path.exists("outputs"):
        os.mkdir("outputs")
    if os.path.exists("obsolete"):
        os.rmdir("obsolete")

    # Move files
    path = os.getcwd()
    os.rename(path + "/pdb" + args.pdb_id + ".ent", path + "/pdb_folder/pdb" + args.pdb_id + ".ent")
    os.rename(path + "/" + args.pdb_id + "_ouput.txt", path + "/outputs/" + args.pdb_id + "_ouput.txt")
    
    if args.statistical: 
        os.rename(path + "/" + args.pdb_id + "_barplot_aa.png", path + "/outputs/" + args.pdb_id + "_barplot_aa.png")
        os.rename(path + "/" + args.pdb_id + "_barplot_atoms.png", path + "/outputs/" + args.pdb_id + "_barplot_atoms.png")
        os.rename(path + "/" + args.pdb_id + "_points.png", path + "/outputs/" + args.pdb_id + "_points.png")
        os.rename(path + "/" + args.pdb_id + "_surface.pse", path + "/outputs/" + args.pdb_id + "_surface.pse")
        os.rename(path + "/" + args.pdb_id + "_neighbors.pse", path + "/outputs/" + args.pdb_id + "_neighbors.pse")

