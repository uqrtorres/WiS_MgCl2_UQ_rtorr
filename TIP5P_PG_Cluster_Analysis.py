##   Import MDAnalysis, numpy and others   ##
import MDAnalysis as mda
import numpy as np
from MDAnalysis import analysis
from MDAnalysis.analysis import distances
import pandas as pd

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')

from matplotlib import pyplot as plt


##    Load the data/dump file of interest and lists the total atoms in the system   ##
trajectories = ['NPT_TIP5P_2018_6.3M_MgCl2_H2O_dump.lammpstrj', 
                'NVT_TIP5P_2018_6.3M_MgCl2_H2O_dump.lammpstrj',
                '55ns_NVT_TIP5P_2018_6.3M_MgCl2_H2O_dump.lammpstrj']

u = mda.Universe('NPT_TIP5P_2018_6.3M_MgCl2_H2O.data',
                 trajectories,
                 topology_format="data", format="lammpsdump",
                 atom_style = 'id resid type x y z') 

print(u.trajectory, "\n")


##    Groups each atom type and lists total number of each atoms   ## 
#Group the Mg Atoms
Mg = u.select_atoms('type 4')
n_Mg = len(Mg)
print('There are',n_Mg, 'Mg atoms')

#Group the Cl atoms
Cl = u.select_atoms('type 1')
n_Cl = len(Cl)
print('There are',n_Cl, 'Cl atoms')

#Group the O atoms
O = u.select_atoms('type 5')
n_O = len(O)
print('There are',n_O, 'O atoms\n')


##   Group the C atoms based on electrode type (cathode and anode)   ##
#Establish cutoff for cluster coordination
inp_cutoff = float(input("Please enter your desired cutoff for a Mg cluster (in angstroms):"))
cutoff = inp_cutoff
print("The cutoff is", cutoff, "angstroms\n")


##   Creates blank arrays for the collection of the coordination data  ##
Cl_array = [] #Makes an empty list for a Cl-coordination summary array
O_array = [] #Makes an empty list for a O-coordination summary array
Coord_array = [] #Makes an empty list for a coordination number summary array
Type_array = [] #Makes an empty list for a cluster-type array


##   Creates a function that calls a loop to extract data every 10 frames  ##
def cluster_loop(x,y,z):
    yes = [x,y,z]
    for i in range(x,y,z):
        ##  Skips initialisation (i = 0) and duplicate frames (i = 101, 202, etc.) if applicable
        if i == 0:
            i = i + 1
        elif i == 101:
            i = i + 1
        elif i == 302:
            i = i + 1
    
        ##    Lists the current frame   ##
        time = i
        u.trajectory[time]
        print("\n" + "Timestep:", time)
        print("These are the results for", u.trajectory[time], "\n")
        

        ##    Calculate the distances between Mg & Electrodes   ## 

        #Create an array that lists the distances between Mg and Cl
        MgCl_dist_arr = distances.distance_array(Mg.positions, # reference
                                                 Cl.positions, # configuration
                                                 box=u.dimensions)

        #Create an array that lists the distances between Mg and O
        MgO_dist_arr = distances.distance_array(Mg.positions, # reference
                                                O.positions, # configuration
                                                box=u.dimensions)

        
        ##   Establish the variables for each clustering group type   #
        neutral = 0
        cationic = 0
        anionic = 0


        ##    Counts the atoms that satisfy the cutoff for clustering   ##
        Cl_cluster = np.count_nonzero(MgCl_dist_arr <= cutoff, axis = 1)
        O_cluster = np.count_nonzero(MgO_dist_arr <= cutoff, axis = 1)


        ##   Counts how many Cl atoms are coordinated to each Mg atom, and adds it to an array   ##
        print() #Blank space
        for i in range(0,5):
            if i == 2:
                neutral = np.count_nonzero(Cl_cluster == 2)
            elif i <= 1:
                cationic = np.count_nonzero(Cl_cluster <= 1)
            else:
                anionic = np.count_nonzero(Cl_cluster >= 3)
            print("There are", np.count_nonzero(Cl_cluster == i), "Mg atoms that are coordinated to", i, "Cl atoms")
            Cl_array.append(np.count_nonzero(Cl_cluster == i))
        Type_array.extend([neutral,cationic,anionic])
        

        ##   Counts how many O atoms are coordinated to each Mg atom, and adds it to an array   ##
        print() #Blank space
        for i in range(0,7):
            print("There are", np.count_nonzero(O_cluster == i), "Mg atoms that are coordinated to", i, "O atoms")
            O_array.append(np.count_nonzero(O_cluster == i))


        ##    Summarises the cluster groupings and coordinations   ##
        
        ##   Calculates and summarises the total coordination of each Mg cluster, and adds it to an array   ##
        coordination = Cl_cluster + O_cluster
        print() #Blank space
        for i in range(0,7):
            coord_num = np.count_nonzero(coordination == i)
            if coord_num != 0:
                print("There are", coord_num, "Mg clusters coordinated to", i, "atoms")
                Coord_array.append(np.count_nonzero(coordination == i))
            elif i == 5 and coord_num == 0:
                print("There are", coord_num, "Mg clusters coordinated to", i, "atoms")
                Coord_array.append(np.count_nonzero(coordination == i))
            else:
                continue
                
    
        ##   Summarises the types of clusters near the electrodes   #
        print() #Blank space
        print() #Blank space
        print("In total, there are:", neutral, "neutral clusters,", cationic, "cationic clusters, and", anionic, "anionic clusters \n")
        print(neutral, cationic, anionic)

        #Check to see if total clusters = 240.
        if neutral + cationic + anionic != 240:
            print("ERROR: Total clusters does not equal to 240. Check array and code")


#Calls the function cluster_loop and examines the clusters in the specified range, whilst skipping the duplicate frames
cluster_loop(0,101,10)
cluster_loop(100,201,10)
cluster_loop(201,402,10)


##   Creates matrices for all coordination data for data analysis purposes   ##
#Create a matrix for the Cl-coordination numbers for easy transfer to excel
rows = 43 #Create a variable that represents each row of data, replace with how many rows are present in the data
print() #Blank space
print("This is a matrix consisting of all the data for the Cl-coordinations")
Cl_matrix = np.array(Cl_array).reshape(rows,5)
print(Cl_matrix) #Summarises the coordination in arrays within a matrix for easier transfer to excel

#Create a matrix for the O-coordination numbers for easy transfer to excel
print() #Blank space
print("This is a matrix consisting of all the data for the O-coordinations")
O_matrix = np.array(O_array).reshape(rows,7)
print(O_matrix) #Summarises the coordination in arrays within a matrix for easier transfer to excel

#Create a matrix for the overall coordination numbers for easy transfer to excel
print() #Blank space
print("This is a matrix consisting of all the data for the overall coordinations")
Coord_matrix = np.array(Coord_array).reshape(rows,3)
print(Coord_matrix) #Summarises the coordination in arrays within a matrix for easier transfer to excel

#Create a matrix for the cluster types for easy transfer to excel
print() #Blank space
print("This is a matrix consisting of all the data for the cluster types")
Type_matrix = np.array(Type_array).reshape(rows,3)
print(Type_matrix) #Summarises the coordination in arrays within a matrix for easier transfer to excel

# Use ExcelWriter to handle multiple sheets
with pd.ExcelWriter('TIP5P Bulk Cluster Analysis.xlsx', engine='openpyxl') as writer:
    pd.DataFrame(Cl_matrix).to_excel(writer, sheet_name='Cl Coordination', index=False, header=False)
    pd.DataFrame(O_matrix).to_excel(writer, sheet_name='O Coordination', index=False, header=False)
    pd.DataFrame(Coord_matrix).to_excel(writer, sheet_name='Overall Coordination', index=False, header=False)
    pd.DataFrame(Type_matrix).to_excel(writer, sheet_name='Cluster Type', index=False, header=False)
