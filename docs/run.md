# Run SA2py

## Setup your environnement 

Clone the repository   
```sh
git clone https://github.com/ElouanBethuel/SASA.git
```
Move to the project directory   
```sh
cd SASA
```
Install [conda](https://docs.conda.io/en/latest/)   
Install mamba 
```
conda install mamba -n base -c conda-forge
```
Create and activate the conda environnement   
```sh
conda env create -f sasa_project.yml
conda activate sasa_project 
``` 
<br>

## Run the program
Run with the program :  
```sh
python main.py 7kh5 100
```
Execute the script main.py with the pdf file 7kh5.pdb and with 100 points to model the solvation sphere.   
You can replace 7kh5 by any another PDB ID. This program work with many PDB files.If the PDB file format is not compatible an error message is raise.   

To create also graphics and pymol files for visualisation add the argument -s

```sh
python main.py 7kh5 100 -s 
```

<br>

## Outputs

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

### Files generated

With by default arguments : 
- The pdb file 
- A text file with the solvent accessibility of each atom   

<br> 

With additional arguments (-s) : 
- A png file (graph) to visualize the points created around a atom
- A png file (graph) showing solvent accessibility by atom category
- A png file (graph) showing solvent accessibility by amino acid category
- A pymol file for showing the accessibility of the protein surface
- A pymol file for showing the neighbors selection 

### This files are stored in two folders
- Outputs : for stored all outputs
- PDB :  for stored all pdb files

### The file name format is always the same (exemple with the 7kh5 pdb ID)
- pdb file : pdb7kh5.ent (pdb + pdb_id + .ent)
- png file : pdb_id + name_file + .png 
- texte file : pdb_id + .txt 
- pymol file (neihgbors) : pdb_id + _neighbors.pse 
- pymol file (surface) : pdb_id + _surface.pse

### To visualize pymol files with the shell 
```sh 
pymol outputs/7kh5_surface.pse 
```
![image](.readme_images/pymol_7kh5_surface.png)

This image shows the 7kh5 protein with its surface coloured according to its accessibility to the solvent: the blue parts are the least accessible and the red ones the most accessible. 

<br>
