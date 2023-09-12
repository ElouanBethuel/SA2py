
# SASA Project (Solvent-Accessible Surface Area)

 Program to calculate the Solvent-Accessible Surface Area (SASA) of a protein.   
 Give the the area in Square Angstrom for each atom and for the protein.
 This program is devide in two script python :   
 - with all functions   
 - with the main function   

 

## How tu use the program

<<<<<<< HEAD
Clone the repository   
```sh
git clone https://github.com/ElouanBethuel/SASA.git
```
Move to the project directory   
```sh
cd SASA
```
Create and activate the conda environnement   
```sh
conda env create -f sasa_project.yml
conda activate sasa_project 
```
=======
Clone the repository  
``git clone ``
Move to the new directory   
``cd ``   
Create and activate the conda environnement   
``conda env create -f sasa_project.yml   
  conda acvitate sas_project``   
  
>>>>>>> 57b9d3b7f9899efe90dcb40eff9690510af407fd
Execute the script   
```sh
python project.py 7kh5
```
You can replace 7kh5 by any another PDB ID. This program work with many PDB files. If the PDB file format is not compatible an error message is raise. 

## Ouputs

IIn the shell, the Solvent Accessible Surface Area of the protein is displayed :

```
Calculates the solvent accessible surface area from the 7kh5 PDB file :

Downloading the PDB file	
Loading the PDB file

Downloading PDB structure '7kh5'...


Calculation ongoing, please wait a few seconds...

=================
Final result :	
The solvent accessible surface area of the protein : 6503.51 Å² 
```

### Many files are also generated : 
- The pdb file 
- A text file with the solvent accessibility of each atom
- A png file (graph) to visualize the points created around a atom
- A png file (graph) showing solvent accessibility by atom category
- A png file (graph) showing solvent accessibility by amino acid category
- A pymol file for showing the accessibility of the protein surface
- A pymol file for showing the neighbors selection 

This files are stored in two folders : 
- Outputs : for stored all outputs
- PDB :  for stored all pdb files

The file name format is always the same (exemple with the 7kh5 pdb ID):
- pdb file : pdb7kh5.ent (pdb + pdb_id + .ent)
- png file : pdb_id + name_file + .png 
- texte file : pdb_id + .txt 
- pymol file (neihgbors) : pdb_id + _neighbors.pse 
- pymol file (surface) : pdb_id + _surface.pse

For open pymol files with the shell 
```sh 
pymol outputs/7kh5_surface.pse 
```
![image](.readme_images/pymol_7kh5_surface.png)

## Results 


