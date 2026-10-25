import MDAnalysis as mda
import numpy as np
from MDAnalysis import analysis
from MDAnalysis.analysis import distances

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')

from matplotlib import pyplot as plt



##    Load the data/dump file of the beginning and end state of the equilibriation, and lists the total atoms in both systems (cross-check)   ##
u = mda.Universe('NPT_TIP5P_2018_6.3M_MgCl2_H2O.data', 
                 "NPT_TIP5P_2018_6.3M_MgCl2_H2O_dump.lammpstrj",
                 topology_format="data", format="lammpsdump",
                 atom_style = 'id resid type x y z') #Beginning 5ns of equilibration

u2 = mda.Universe('NPT_TIP5P_2018_6.3M_MgCl2_H2O.data', 
                 "55ns_NVT_TIP5P_2018_6.3M_MgCl2_H2O_dump.lammpstrj",
                 topology_format="data", format="lammpsdump",
                 atom_style = 'id resid type x y z') #Last 10ns of equilibration



##    Groups each atom type and lists total number of each atoms   ## 

#Group the Mg Atoms for both starting and end points of equilibration
Start_Mg = u.select_atoms('type 4')
Start_n_Mg = len(Start_Mg)
print('There are',Start_n_Mg, 'Mg atoms in the beginning')

End_Mg = u2.select_atoms('type 4')
End_n_Mg = len(End_Mg)
print('There are',End_n_Mg, 'Cl atoms in the end \n')

#Group the Cl atoms for both starting and end points of equilibration
Start_Cl = u.select_atoms('type 1')
Start_n_Cl = len(Start_Cl)
print('There are',Start_n_Cl, 'Cl atoms in the beginning')

End_Cl = u2.select_atoms('type 1')
End_n_Cl = len(End_Cl)
print('There are',End_n_Cl, 'Mg atoms in the end \n')

#Group the O atoms for both starting and end points of equilibration
Start_O = u.select_atoms('type 5')
Start_n_O = len(Start_O)
print('There are',Start_n_O, 'O atoms in the beginning')

End_O = u2.select_atoms('type 5')
End_n_O = len(End_O)
print('There are',End_n_O, 'O atoms in the end \n')

#Establish cutoff for coordination
cutoff = 3.0
print("The cutoff is", cutoff, "angstroms")



##   Creates blank arrays for the collection of the coordination data  ##

Start_Cl_array = [] #Makes an empty list for a Mg-coordination summary array
Start_O_array = [] #Makes an empty list for a H-coordination summary array
Start_Coord_array = [] #Makes an empty list for a coordination number summary array
Start_array = [] #Makes an empty list for a cluster-type array

End_Cl_array = [] #Makes an empty list for a Mg-coordination summary array
End_O_array = [] #Makes an empty list for a H-coordination summary array
End_Coord_array = [] #Makes an empty list for a coordination number summary array
End_array = [] #Makes an empty list for a cluster-type array



##  Extracts data from the start of the equilibration  ##

##    Lists the current frame   ##
time = 1
u.trajectory[time]
print("\n" + "Timestep: 0.05ns")
print("These are the total coordinations of the clusters for the beginning of the equilibration \n")

##    Calculate the distances between Mg & Cl, and Mg & H/O   ## 

#Create an array that lists the distances between Mg and Cl
Start_MgCl_dist_arr = distances.distance_array(Start_Mg.positions, # reference
                                         Start_Cl.positions, # configuration
                                         box=u.dimensions)

#Create an array that lists the distances between Mg and O
Start_MgO_dist_arr = distances.distance_array(Start_Mg.positions, # reference
                                        Start_O.positions, # configuration
                                        box=u.dimensions)
    
    
##    Counts the atoms that satisfy the cutoff for clustering   ##
Start_Cl_cluster = np.count_nonzero(Start_MgCl_dist_arr <= cutoff, axis = 1)
Start_O_cluster = np.count_nonzero(Start_MgO_dist_arr <= cutoff, axis = 1)
    
#Remove the '#' from the below code to check the arrays for each cluster
#print(Start_Cl_cluster, '\n')
#print(Start_O_cluster, '\n')
    
#Calculates and summarises the total coordination of each Mg cluster, and adds it to an array
Start_coordination = Start_Cl_cluster + Start_O_cluster
print(Start_coordination, '\n')



##  Extracts data from the end of the equilibration  ##

##    Lists the current frame   ##
time = 200
u2.trajectory[time]
print("\n" + "Timestep: 35ns")
print("These are the total coordinations of the clusters after 35ns of equilibration \n")

##    Calculate the distances between Mg & Cl, and Mg & H/O for the end of equilibration ## 

#Create an array that lists the distances between Mg and Cl
End_MgCl_dist_arr = distances.distance_array(End_Mg.positions, # reference
                                         End_Cl.positions, # configuration
                                         box=u2.dimensions)

#Create an array that lists the distances between Mg and O
End_MgO_dist_arr = distances.distance_array(End_Mg.positions, # reference
                                        End_O.positions, # configuration
                                        box=u2.dimensions)
    
    
##    Counts the atoms that satisfy the cutoff for clustering for the end of equilibration   ##
End_Cl_cluster = np.count_nonzero(End_MgCl_dist_arr <= cutoff, axis = 1)
End_O_cluster = np.count_nonzero(End_MgO_dist_arr <= cutoff, axis = 1)
    
#Remove the '#' from the below code to check the arrays for each cluster
#print(End_Cl_cluster, '\n')
#print(End_O_cluster, '\n')
    
#Calculates and summarises the total coordination of each Mg cluster at the end of equilibration, and adds it to an array
End_coordination = End_Cl_cluster + End_O_cluster
print(End_coordination, '\n')


##    Compares the coordinations of the clusters between the start and end states of the equilibration + lists each indices' differences   ## 
Diff_array = End_coordination - Start_coordination
print("These are the differences in coordination numbers for the Mg clusters (between the start and end of equilibration)")
print(Diff_array, '\n')

Change = np.count_nonzero(Diff_array)
No_Change = 240 - Change

if Change + No_Change != 240:
    print("There is an error in the count!")

# Lists the indices of the Mg clusters that change/don't change during equilibration #
change_indices = [index for (index, item) in enumerate(Diff_array) if item != 0]
print("If you want to monitor the Mg clusters that change during the equilibration, monitor the following indices:")
print(change_indices, '\n')

nochange_indices = [index for (index, item) in enumerate(Diff_array) if item == 0]
print("If you want to monitor the Mg clusters that don't change during the equilibration, monitor the following indices:")
print(nochange_indices)



