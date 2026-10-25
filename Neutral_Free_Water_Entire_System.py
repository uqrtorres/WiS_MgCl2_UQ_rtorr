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
                 'Dump_3rd_30ns_Equilibriate_NVT_PG_MgCl2.lammpstrj',
                 topology_format="data", format="lammpsdump",
                 atom_style = 'id resid type x y z') 

print(u.trajectory, "\n")


##   Group the C atoms based on electrode type (cathode and anode)   ##
#Establish cutoff for cluster coordination
inp_cutoff = float(input("Please enter your desired cutoff for a Mg cluster (in angstroms):"))
cutoff = inp_cutoff
print("The cutoff is", cutoff, "angstroms\n")

#Establish cutoff for the clusters within a distance of the electrode walls
inp_ele_cutoff = float(input("Please enter your desired cutoff for a Mg cluster being at an electrode's interface (in angstroms):"))
ele_cutoff = inp_ele_cutoff
print("The cutoff is", ele_cutoff, "angstroms\n")

##   Creates blank arrays for the collection of the free water evaluation ##
Total_Water = []
Free_Water = []


##   Creates a function that calls a loop to extract data every 10 frames  ##
def cluster_loop(x,y,z):
    yes = [x,y,z]
    for i in range(x,y,z):
        ##  Skips initialisation (i = 0) and duplicate frames (i = 101, 202, etc.) if applicable
        if i == 0:
            i = i + 1

        ##  Creates/resets the arrays for free water analysis  ##
        O_array = [] #Makes an empty list for a O-coordination summary array
        Free_O_Array = [] #Makes an empty list for a free O-coordination summary array
        H_array = [] #Makes an empty list for a H-coordination summary array
        Free_H_Array = [] #Makes an empty list for a free H-coordination summary array
    
        ##    Lists the current frame   ##
        time = i
        u.trajectory[time]
        print("\n" + "Timestep:", time)
        print("These are the results for", u.trajectory[time], "\n")


        ##  Identifies the total O atoms (water) that are not coordinated to a Mg ion  ##
        #Group the O atoms
        O = u.select_atoms('type 4')
        n_O = len(O)
        O_array.append(n_O)
        print('There are',n_O, 'O atoms')

        #Group of O atoms that are not within the cutoff of Mg and within the cutoff of the anode
        Free_O = u.select_atoms(f"type 4 and (not around {cutoff} type 3)")
        n_Free_O = len(Free_O)
        Free_O_Array.append(n_Free_O)
        print('There are',n_Free_O, 'free water O atoms\n')

        #Places the indices of each free O-atom into one array
        for o in Free_O:
            #print(o.index) #Uncomment for sanity check
            Free_O_Array.append(int(o.index))
        #print("The length of the O-array is", len(Free_O_Array),'\n')


        ##  Identifies the total H atoms (water) that are not coordinated to a Cl ion, and finds their attached O parent-atom  ##
        #Group of H atoms that are within the cutoff of the anode
        H = u.select_atoms('type 2')
        n_H = len(H)
        print('There are',n_H, 'H water atoms')

        #Group of H atoms that are not within the cutoff of Cl and within the cutoff of the anode
        Free_H = u.select_atoms(f"type 2 and (not around {cutoff} type 1)")
        n_Free_H = len(Free_H)
        print('There are',n_Free_H, 'free H water atoms\n')

        #Finds the O parent-atom for each 'free H' and places the parent O indices into an array (Anode)
        for h in Free_H:
            #print(h.index, h.bonded_atoms.indices)
            Free_H_Array.append(int((h.bonded_atoms.indices)))
            #print(u.atoms[h.bonded_atoms.indices])
        #print("The length of the H-array is", len(Free_H_Array),'\n')
        #print(len(Free_O_Array))
        #print(len(Free_H_Array))

        ##  Finds the shared O-indices and removes any duplicates to create the final free water list (Anode) ##
        #Find shared O-indices between both lists
        Reference_O = set(Free_O_Array)
        Filtered_H_Free_Water = [item for item in Free_H_Array if item in Reference_O]
        #print(Filtered_H_Free_Water) #Uncomment for sanity Check

        #Remove duplicates from the filtered list to get final free water count
        Final_Free_Water = list(dict.fromkeys(Filtered_H_Free_Water))
        #print(Final_Free_Water, '\n') #Uncomment for sanity check
        #print(n_O) #Uncomment for sanity check
        Total_Water.append(n_O)
        #print(len(Final_Free_Water)) #Uncomment for sanity check
        Free_Water.append(len(Final_Free_Water))


##  Run the cluster analysis code and print the outputs  ##
cluster_loop(0,101,2)

## Check the array for the free and total O
print(Total_Water, '\n')
print(Free_Water)


## Collect the total and free water data for excel
Data = {
    "Total O": Total_Water,
    "Total Free O": Free_Water
}

#Use ExcelWriter to write the free water data into an excel sheet
with pd.ExcelWriter('PG Whole System Free Water (Neutral).xlsx', engine='openpyxl') as writer:
    pd.DataFrame(Data).to_excel(writer, sheet_name='Bulk Water', index=True)