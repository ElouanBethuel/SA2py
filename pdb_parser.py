
import numpy as np
import random
import copy
import Bio
from Bio.PDB import PDBList
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d

RADIUS_H2O = 1.4
NB_SONDES = 92


# Créez une instance de la classe PDBList
#pdbl = PDBList()
# Spécifiez le code PDB du fichier que vous souhaitez télécharger
#pdb_code = "1ATP" 
# Utilisez la méthode get_pdb_file pour télécharger le fichier PDB
#pdbl.retrieve_pdb_file(pdb_code, file_format="pdb")



# script to parse PDB files and add Radius and area of solvated sphere 

def parse_pdb_file(path_pdb_file):

	info = []

	with open(path_pdb_file , "r") as pdb_file :
		
		for line in pdb_file:
			if line.startswith("ATOM") == True :
				split_line = line.split()
				
				if split_line[1] != "H":
					new_line = [split_line[1], split_line[2], split_line[3], split_line[6], split_line[7], split_line[8]]
					
					# add Van der Wadd radius
					
					if new_line[1].find("N") != -1:
						new_line.append("1.55")
						
					if new_line[1].find("O") != -1:
						new_line.append("1.52")
					
					if new_line[1].find("S") != -1:
						new_line.append("1.8")
					
					if new_line[1].find("C") != -1:
						new_line.append("1.7")
					
					info.append(new_line)

	return info
	



# function for generated 92 sondes

def test_point(x,y,z, radius):

	r = radius + RADIUS_H2O

	theta = np.random.uniform(0, np.pi, NB_SONDES)
	phi = np.random.uniform(0, 2*np.pi, NB_SONDES)
	
	test_x = x + r * np.sin(theta) * np.cos(phi)
	test_y = y + r * np.sin(theta) * np.sin(phi)
	test_z = z + r * np.cos(theta)
	
	point = np.array([test_x, test_y, test_z]).T

	return point
	
	


def test_point_graphic(x,y,z, radius):
	
	r = radius + RADIUS_H2O
	
	theta = np.random.uniform(0, np.pi, NB_SONDES)
	phi = np.random.uniform(0, 2*np.pi, NB_SONDES)

	test_x = x + r * np.sin(theta) * np.cos(phi)
	test_y = y + r * np.sin(theta) * np.sin(phi)
	test_z = z + r * np.cos(theta)
	
	plt.figure()
	axes = plt.axes(projection="3d")
	
	# affiche en rouge l'atome 
	axes.scatter(x, y, z, color='red')
	
	# affiche en bleu les sondes autour de l'atome
	axes.scatter(test_x, test_y, test_z)
	
	plt.show()
	point = np.array([test_x, test_y, test_z]).T
	
	return points




def all_atomes_sondes(info_pdb):

	list_all_points = {}
	
	for atome in info_pdb:
	
		point = test_point(float(atome[3]), float(atome[4]), float(atome[5]), float(atome[6]))
		list_all_points[atome[0]]  = point

	return list_all_points




# retourne la liste des voisins d'un atome à partir de ses coordonnées 
# et d'un  seuil (treeshold : distance en Angstrom) 

def neigbords(num_atome, pdb, threeshold):

	info_pdb = copy.deepcopy(pdb)
	
	x = float(info_pdb[num_atome-1][3])
	y = float(info_pdb[num_atome-1][4])
	z = float(info_pdb[num_atome-1][5])
	
	del(info_pdb[num_atome-1])
	list_neigbords = []
	
	for atome in info_pdb:
			
			x_nb = float(atome[3])
			y_nb = float(atome[4])
			z_nb = float(atome[5])
			
			d = np.sqrt((x - x_nb)**2 + (y - y_nb)**2 + (z - z_nb)**2)
			
			if d < threeshold : 
				list_neigbords.append(atome)
	
	return list_neigbords 
	
	



def all_neigbords(info_pdb , threeshold):

	list_all_neigbords = []
	
	for i in range(len(info_pdb)):
		list_all_neigbords.append(neigbords(i+1, info, threeshold))

	return list_all_neigbords 
	




# CODE POUR PLOTER TT LES POINTS DE LA PROTEINE EN BLEUE
# ET LES VOISINS DE L'ATOME SELECTIONNE EN ROUGE 


def plot_proteine(info, num_atome, distance):

	list_x_pdb = []
	list_y_pdb = []
	list_z_pdb = []

	for atome in info:
	
		list_x_pdb.append(float(atome[3]))
		list_y_pdb.append(float(atome[4]))
		list_z_pdb.append(float(atome[5]))

	array_x_pdb = np.asarray(list_x_pdb)
	array_y_pdb = np.asarray(list_y_pdb)
	array_z_pdb = np.asarray(list_z_pdb)

	list_atomes_neigbords = neigbords(num_atome, info, distance)

	list_x_neigbords = []
	list_y_neigbords = []
	list_z_neigbords = []

	for atome in list_atomes_neigbords:
	
		list_x_neigbords.append(float(atome[3]))
		list_y_neigbords.append(float(atome[4]))
		list_z_neigbords.append(float(atome[5]))
		
	array_x_neigbords = np.asarray(list_x_neigbords)
	array_y_neigbords = np.asarray(list_y_neigbords)
	array_z_neigbords = np.asarray(list_z_neigbords)

	plt.figure()
	axes = plt.axes(projection="3d")
	axes.scatter(array_x_pdb, array_y_pdb, array_z_pdb, color="blue")
	axes.scatter(float(info[num_atome][3]), float(info[num_atome][4]), float(info[num_atome][5]), s= 300, color="green")
	axes.scatter(array_x_neigbords, array_y_neigbords, array_z_neigbords, s=200, color='red')
	plt.show()




# fonction pour déterminer le pourcentage de sondes accessibles au solvant 
# pour chaque atome de la protéine 

def access_solvant(info, list_all_neigbords, list_all_atomes_sondes):

	sum_acc_area = 0
	i = 0 

	for atome in info:
		
		sondes = list_all_atomes_sondes[atome[0]]
		sondes_inacessibles = copy.deepcopy(sondes)
		neigbords = list_all_neigbords[i]

		i += 1 
		nb_sondes = 0 
		
		for neigbord in neigbords:
			
			x_n = float(neigbord[3])
			y_n = float(neigbord[4])
			z_n = float(neigbord[5])
			
			point_n = np.array([x_n, y_n, z_n])
			
			# calcule les distances euclidiennes entre l'atome voisin et les sondes de l'atome 
			distances = np.linalg.norm(point_n -sondes, axis=1)

			points_inacessibles = distances <= 1.4 + float(atome[6])
			indices_points_inacessibles = np.where(points_inacessibles)[0]
			sondes_inacessibles[indices_points_inacessibles] = -1
			
		# divise par 3 car 3 coordonnées x,y,z  pour chaque sonde 
		nb_sondes_acessibles = (np.size(sondes)/3) - np.count_nonzero(sondes_inacessibles == -1)/3
		prct_acc_atome = ((nb_sondes_acessibles / (np.size(sondes)/3)) * 100)
		
		area_acc_atome = (prct_acc_atome/100) * 4*np.pi*((1.4+float(atome[6]))**2)
		sum_acc_area = sum_acc_area + area_acc_atome
		
		print(f"{atome} \t")
		print(f"Le nombre de sondes accessibles au solvant est : {nb_sondes_acessibles} % \t")
		print(f"Le pourcentage d'accessibilité de l'atome est : {round(prct_acc_atome,2)} % \t")
		print(f"L'air d'accessibilité de l'atome est de : {round(area_acc_atome,2)} Angstrom carré \n")
	
	sum_acc_area = round(sum_acc_area, 2)
	print(f"L'accessibilité au solvant de la protéine est de {sum_acc_area} Angstrom carré")
			



info = parse_pdb_file("3i40.pdb")
list_all_atomes_sondes = all_atomes_sondes(info)
list_all_neigbord = all_neigbords(info , 10)
access_solvant(info, list_all_neigbord, list_all_atomes_sondes)

plot_proteine(info, 5, 10)








