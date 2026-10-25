import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis import analysis
from MDAnalysis.analysis import distances

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')

from matplotlib import pyplot as plt



##    Load the data/dump file of the entire 35ns equilibration   ##

trajectory_files = ['NPT_TIP5P_2018_6.3M_MgCl2_H2O_dump.lammpstrj',
                    'NVT_TIP5P_2018_6.3M_MgCl2_H2O_dump.lammpstrj',
                    '55ns_NVT_TIP5P_2018_6.3M_MgCl2_H2O_dump.lammpstrj']

u = mda.Universe('NPT_TIP5P_2018_6.3M_MgCl2_H2O.data', trajectory_files,
                 topology_format="data", format="lammpsdump",
                 atom_style = 'id resid type x y z') 

print(u.trajectory)



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

#Establish cutoff for coordination
inp_cutoff = float(input("Please enter your desired cutoff (in angstroms):"))
cutoff = inp_cutoff
print("The cutoff is", cutoff, "angstroms\n")



##  Establish which Mg cluster the program will monitor  ##
residue = 7460 + int(input('Enter the Mg index you want to monitor (starting from 0):')) #Identifies which atom is the monitored Mg index 
Mg_index = int(residue - 7460) #Stores the monitored Mg index
print() #Blank space
print('The specific Mg atom in the universe will be atom', residue + 1)
print('Check below output to see if the residue matches the Mg index (where residue = Mg index value + 1)')
print(u.atoms[residue], '\n')



##  Create arrays that store the Mg-Cl, Mg-O and overall Coordinations for Mg[index]
Cl_evolve_time = []
O_evolve_time = []
Coord_evolve_time = []


## Create a loop function to evaluate the Mg-cluster of interest and its coordinations over time  ##
for i in range(1,402, 1):
    #Load the trajectory at timeframe i
    u.trajectory[i]

    #Calculate the distances between Mg and O/Cl
    MgCl_dist_arr = distances.distance_array(Mg.positions, # reference
                                            Cl.positions, # configuration
                                            box=u.dimensions)

    MgO_dist_arr = distances.distance_array(Mg.positions, # reference
                                            O.positions, # configuration
                                            box=u.dimensions)

    #Assesses the Cl-coordinations to the indexed Mg and adds the data to an array
    Cl_cluster = np.count_nonzero(MgCl_dist_arr <= cutoff, axis = 1)
    Cl_evolve_time.append(int(Cl_cluster[Mg_index]))

    #Assesses the O-coordinations to the indexed Mg and adds the data to an array
    O_cluster = np.count_nonzero(MgO_dist_arr <= cutoff, axis = 1)
    O_evolve_time.append(int(O_cluster[Mg_index]))

    #Assesses the overall coordinations to the indexed Mg and adds the data to an array
    Cluster_coord = Cl_cluster[Mg_index] + O_cluster[Mg_index]
    Coord_evolve_time.append(int(Cluster_coord))

    

##   Create a dataframe using Panda to tabulate the data  ##
df = pd.DataFrame()
df['Frame'] = range(1,402,1)
df['Cl-Coordinations'] = Cl_evolve_time
df['O-Coordinations'] = O_evolve_time
df['Overall Coordinations'] = Coord_evolve_time
display(df) #Gives a preview of the table

# Create an excel file with the data table
file = input('Please enter the name for the data file (add .xlsx at the end):')
df.to_excel(file, index=False)

print('File is saved! Check directory for the data.')

# Remove # to see the entire arrays for the coordinations (to cross check with table)
#print(Cl_evolve_time, '\n')
#print(O_evolve_time, '\n')
#print(Coord_evolve_time, '\n')



        

