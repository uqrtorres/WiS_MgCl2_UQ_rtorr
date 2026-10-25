##   Import MDAnalysis, numpy and others   ##
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis import analysis
from MDAnalysis.analysis import distances

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')

from matplotlib import pyplot as plt


##    Load the data/dump file of interest and lists the total atoms in the system   ##

u = mda.Universe('3rd_30ns_Equilibriate_PG_6.3M_MgCl2_nvt.data',
                 'Dump_3rd_2.5V_30ns_PG_MgCl2.lammpstrj',
                 topology_format="data", format="lammpsdump",
                 atom_style = 'id resid type x y z') 

print(u.trajectory, "\n")


##    Groups each atom type and lists total number of each atoms   ## 
#Group the Mg Atoms
Mg = u.select_atoms('type 3')
n_Mg = len(Mg)
print('There are',n_Mg, 'Mg atoms')

#Group the Cl atoms
Cl = u.select_atoms('type 1')
n_Cl = len(Cl)
print('There are',n_Cl, 'Cl atoms')

#Group the O atoms
O = u.select_atoms('type 4')
n_O = len(O)
print('There are',n_O, 'O atoms\n')


##   Group the C atoms based on electrode type (cathode and anode)   ##
#Anode
C_an = u.select_atoms('type 5 and index 5196:5875')
n_C_An = len(C_an)
#print(C_anode) #Uncomment for sanity check - lists the indexed C atoms
print("There are", n_C_An, "C atoms in the anode")

#Cathode
C_cat = u.select_atoms('type 5 and index 5876:6555')
n_C_Cat = len(C_cat)
#print(C_cathode) #Uncomment for sanity check - lists the indexed C atoms
print("There are", n_C_Cat, "C atoms in the cathode")


##   Group the C atoms based on electrode type (cathode and anode)   ##
#Establish cutoff for cluster coordination
inp_cutoff = float(input("Please enter your desired cutoff for a Mg cluster (in angstroms):"))
cutoff = inp_cutoff
print("The cutoff is", cutoff, "angstroms\n")

#Establish cutoff for the clusters within a distance of the electrode walls
inp_ele_cutoff = float(input("Please enter your desired cutoff for a Mg cluster being at an electrode's interface (in angstroms):"))
ele_cutoff = inp_ele_cutoff
print("The cutoff is", ele_cutoff, "angstroms\n")

##   Creates blank arrays for the collection of the coordination data  ##
An_Cl_array = [] #Makes an empty list for a Cl-coordination summary array (Anode)
Cat_Cl_array = [] #Makes an empty list for a Cl-coordination summary array (Cathode)
An_O_array = [] #Makes an empty list for a O-coordination summary array (Anode)
Cat_O_array = [] #Makes an empty list for a O-coordination summary array (Cathode)
An_Coord_array = [] #Makes an empty list for a coordination number summary array (Anode)
Cat_Coord_array = [] #Makes an empty list for a coordination number summary array (Cathode)
An_Type_array = [] #Makes an empty list for a cluster-type array (Anode)
Cat_Type_array = [] #Makes an empty list for a cluster-type array (Cathode)

##   Creates a function that calls a loop to extract data every 10 frames  ##
def cluster_loop(x,y,z):
    yes = [x,y,z]
    for i in range(x,y,z):
        ##  Skips initialisation (i = 0) and duplicate frames (i = 101, 202, etc.) if applicable
        if i == 0:
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
        
        #Create an array that lists the distances between Mg and anode C
        MgAn_dist_arr = distances.distance_array(Mg.positions, # reference
                                                C_an.positions, # configuration
                                                box=u.dimensions)

        #Create an array that lists the distances between Mg and cathode C
        MgCat_dist_arr = distances.distance_array(Mg.positions, # reference
                                                C_cat.positions, # configuration
                                                box=u.dimensions)
        

        ##    Counts the atoms that is within the range of the cathode and anode   ##
        #Anode
        #Counts how many Mg clusters are within the cutoff range of the anode
        C_An_cluster = np.count_nonzero(MgAn_dist_arr <= ele_cutoff, axis = 1)
        Mg_An = np.count_nonzero(C_An_cluster)
        print('There are', Mg_An, "Mg atoms that are within the", ele_cutoff, "range of the pristine graphene anode")

        #Places the Mg within range of the anode into a seperate array
        Mg_An_index =[index for (index, item) in enumerate(C_An_cluster) if item != 0]
        #print(len(Mg_An_index)) #Uncomment for sanity check
        #print(Mg_An_index) #Uncomment to see which Mg indices are close to the cathode

        #Cathode
        #Counts how many Mg clusters are within the cutoff range of the cathode
        C_Cat_cluster = np.count_nonzero(MgCat_dist_arr <= ele_cutoff, axis = 1)
        Mg_Cat = np.count_nonzero(C_Cat_cluster)
        print('There are', Mg_Cat, "Mg atoms that are within the", ele_cutoff, "range of the pristine graphene cathode")

        #Places the Mg within range of the cathode into a seperate array
        Mg_Cat_index =[index for (index, item) in enumerate(C_Cat_cluster) if item != 0]
        #print(len(Mg_Cat_index)) #Uncomment for sanity check
        #print(Mg_Cat_index) #Uncomment to see which Mg indices are close to the cathode

        
        ##   Establish the variables for each clustering group type   #
        Cat_neutral = 0
        An_neutral = 0
        Cat_cationic = 0
        An_cationic = 0
        Cat_anionic = 0
        An_anionic = 0

        
        ##    Counts the atoms that satisfy the cutoff for clustering   ##
        #Anode
        An_Cl_cluster = np.count_nonzero(MgCl_dist_arr[Mg_An_index] <= cutoff, axis = 1)
        #print(len(An_Cl_cluster)) #Uncomment for sanity check
        An_O_cluster = np.count_nonzero(MgO_dist_arr[Mg_An_index] <= cutoff, axis = 1)
        #print(len(An_O_cluster)) #Uncomment for sanity check

        #Cathode
        Cat_Cl_cluster = np.count_nonzero(MgCl_dist_arr[Mg_Cat_index] <= cutoff, axis = 1)
        #print(len(Cat_Cl_cluster)) #Uncomment for sanity check
        Cat_O_cluster = np.count_nonzero(MgO_dist_arr[Mg_Cat_index] <= cutoff, axis = 1)
        #print(len(Cat_O_cluster)) #Uncomment for sanity check


        ##   Counts how many Cl atoms are coordinated to each Mg atom, and adds it to an array   ##
        print() #Blank space
        print() #Double space
        for i in range(0,5):
            if i == 2:
                An_neutral = np.count_nonzero(An_Cl_cluster == 2)
                Cat_neutral = np.count_nonzero(Cat_Cl_cluster == 2)
            elif i <= 1:
                An_cationic = np.count_nonzero(An_Cl_cluster <= 1)
                Cat_cationic = np.count_nonzero(Cat_Cl_cluster <= 1)
            else:
                An_anionic = np.count_nonzero(An_Cl_cluster >= 3)
                Cat_anionic = np.count_nonzero(Cat_Cl_cluster >= 3)
            print("There are", np.count_nonzero(An_Cl_cluster == i), "Mg atoms that are coordinated to", i, "Cl atoms (Anode)")
            print("There are", np.count_nonzero(Cat_Cl_cluster == i), "Mg atoms that are coordinated to", i, "Cl atoms (Cathode)\n")
            An_Cl_array.append(np.count_nonzero(An_Cl_cluster == i))
            Cat_Cl_array.append(np.count_nonzero(Cat_Cl_cluster == i))
        An_Type_array.extend([An_neutral, An_cationic, An_anionic])
        Cat_Type_array.extend([Cat_neutral, Cat_cationic, Cat_anionic])
        

        ##   Counts how many O atoms are coordinated to each Mg atom, and adds it to an array   ##
        print() #Double space
        for i in range(0,7):
            #Anode
            print("There are", np.count_nonzero(An_O_cluster == i), "Mg atoms that are coordinated to", i, "O atoms (Anode)")
            An_O_array.append(np.count_nonzero(An_O_cluster == i))
            
            #Cathode
            print("There are", np.count_nonzero(Cat_O_cluster == i), "Mg atoms that are coordinated to", i, "O atoms (Cathode)\n")
            Cat_O_array.append(np.count_nonzero(Cat_O_cluster == i))


        ##    Summarises the cluster groupings and coordinations   ##
        
    
        ##   Calculates and summarises the total coordination of each Mg cluster, and adds it to an array   ##
        An_coordination = An_Cl_cluster + An_O_cluster
        Cat_coordination = Cat_Cl_cluster + Cat_O_cluster
        print() #Blank space
        for i in range(0,7):
            An_coord_num = np.count_nonzero(An_coordination == i)
            Cat_coord_num = np.count_nonzero(Cat_coordination == i)
            if An_coord_num != 0:
                print("There are", An_coord_num, "Mg clusters coordinated to", i, "atoms (Anode)")
                An_Coord_array.append(np.count_nonzero(An_coordination == i))
            elif i == 5 and An_coord_num == 0:
                print("There are", An_coord_num, "Mg clusters coordinated to", i, "atoms (Anode)")
                An_Coord_array.append(np.count_nonzero(An_coordination == i))
            if Cat_coord_num != 0:
                print("There are", Cat_coord_num, "Mg clusters coordinated to", i, "atoms (Cathode)\n")
                Cat_Coord_array.append(np.count_nonzero(Cat_coordination == i))
            elif i == 5 and Cat_coord_num == 0:
                print("There are", Cat_coord_num, "Mg clusters coordinated to", i, "atoms (Cathode)\n")
                Cat_Coord_array.append(np.count_nonzero(Cat_coordination == i))
            else:
                continue
    
        ##   Summarises the types of clusters near the electrodes   #
        print() #Blank space
        #Anode
        print("In total, there are:", An_neutral, "neutral clusters,", An_cationic, "cationic clusters, and", An_anionic, "anionic clusters (Anode)")
        print(An_neutral, An_cationic, An_anionic)

        #Cathode
        print("In total, there are:", Cat_neutral, "neutral clusters,", Cat_cationic, "cationic clusters, and", Cat_anionic, "anionic clusters (Cathode)")
        print(Cat_neutral, Cat_cationic, Cat_anionic, "\n")

##  Run the cluster analysis code and print the outputs  ##
cluster_loop(0,601,10)

##   Creates matrices for all coordination data for data analysis purposes   ##

#Create a variable that corresponds to the number of rows in a matrix
row = 61 #Edit to how many expected rows are in the matrix

#Create a matrix for the Cl-coordination numbers for easy transfer to excel
print() #Blank space
#Anode
print("This is a matrix consisting of all the data for the Cl-coordinations near the Anode")
An_Cl_matrix = np.array(An_Cl_array).reshape(row,5)
print(An_Cl_matrix, "\n") #Summarises the coordination in arrays within a matrix for easier transfer to excel

#Cathode
print("This is a matrix consisting of all the data for the Cl-coordinations near the Cathode")
Cat_Cl_matrix = np.array(Cat_Cl_array).reshape(row,5)
print(Cat_Cl_matrix, "\n") #Summarises the coordination in arrays within a matrix for easier transfer to excel


#Create a matrix for the O-coordination numbers for easy transfer to excel
#Anode
print("This is a matrix consisting of all the data for the O-coordinations near the Anode")
An_O_matrix = np.array(An_O_array).reshape(row,7)
print(An_O_matrix, "\n") #Summarises the coordination in arrays within a matrix for easier transfer to excel

#Cathode
print("This is a matrix consisting of all the data for the O-coordinations near the Cathode")
Cat_O_matrix = np.array(Cat_O_array).reshape(row,7)
print(Cat_O_matrix, "\n") #Summarises the coordination in arrays within a matrix for easier transfer to excel


#Create a matrix for the overall coordination numbers for easy transfer to excel
#Anode
print("This is a matrix consisting of all the data for the overall coordinations near the Anode")
An_Coord_matrix = np.array(An_Coord_array).reshape(row,3)
print(An_Coord_matrix, "\n") #Summarises the coordination in arrays within a matrix for easier transfer to excel

#Cathode
print("This is a matrix consisting of all the data for the overall coordinations near the Cathode")
Cat_Coord_matrix = np.array(Cat_Coord_array).reshape(row,3)
print(Cat_Coord_matrix, "\n") #Summarises the coordination in arrays within a matrix for easier transfer to excel


#Create a matrix for the cluster types for easy transfer to excel
#Anode
print("This is a matrix consisting of all the data for the cluster types (Anode)")
An_Type_matrix = np.array(An_Type_array).reshape(row,3)
print(An_Type_matrix, "\n") #Summarises the coordination in arrays within a matrix for easier transfer to excel

#Cathode
print("This is a matrix consisting of all the data for the cluster types (Cathode)")
Cat_Type_matrix = np.array(Cat_Type_array).reshape(row,3)
print(Cat_Type_matrix, "\n") #Summarises the coordination in arrays within a matrix for easier transfer to excel

# Use ExcelWriter to handle multiple sheets
with pd.ExcelWriter('PG Interface Cluster Analysis (2.5V).xlsx', engine='openpyxl') as writer:
    pd.DataFrame(An_Cl_matrix).to_excel(writer, sheet_name='Anode Cl', index=False, header=False)
    pd.DataFrame(Cat_Cl_matrix).to_excel(writer, sheet_name='Cathode Cl', index=False, header=False)
    pd.DataFrame(An_O_matrix).to_excel(writer, sheet_name='Anode O', index=False, header=False)
    pd.DataFrame(Cat_O_matrix).to_excel(writer, sheet_name='Cathode O', index=False, header=False)
    pd.DataFrame(An_Coord_matrix).to_excel(writer, sheet_name='Anode Overall', index=False, header=False)
    pd.DataFrame(Cat_Coord_matrix).to_excel(writer, sheet_name='Cathode Overall', index=False, header=False)
    pd.DataFrame(An_Type_matrix).to_excel(writer, sheet_name='Anode Type', index=False, header=False)
    pd.DataFrame(Cat_Type_matrix).to_excel(writer, sheet_name='Cathode Type', index=False, header=False)
