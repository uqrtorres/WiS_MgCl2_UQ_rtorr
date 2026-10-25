##   Import MDAnalysis, numpy and others   ##
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis import analysis
from MDAnalysis.analysis import distances
from MDAnalysis import core
from MDAnalysis.core import groups
from MDAnalysis.core.groups import Atom

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
#Anode
C_an = u.select_atoms('type 5 and index 5196:5875')
n_C_An = len(C_an)
#print(C_anode) #Uncomment for sanity check - lists the indexed C atoms
print("There are", n_C_An, "C atoms in the anode")

#Cathode
C_cat = u.select_atoms('type 5 and index 5876:6555')
n_C_Cat = len(C_cat)
#print(C_cathode) #Uncomment for sanity check - lists the indexed C atoms
print("There are", n_C_Cat, "C atoms in the cathode\n")

##   Group the C atoms based on electrode type (cathode and anode)   ##
#Establish cutoff for cluster coordination
inp_cutoff = float(input("Please enter your desired cutoff for a Mg cluster (in angstroms):"))
cutoff = inp_cutoff
print("The cutoff is", cutoff, "angstroms\n")

#Establish cutoff for the clusters within a distance of the electrode walls
inp_ele_cutoff = float(input("Please enter your desired cutoff for a Mg cluster being at an electrode's interface (in angstroms):"))
ele_cutoff = inp_ele_cutoff
print("The cutoff is", ele_cutoff, "angstroms\n")


##  Establish the arrays for the free water and total water near the electrodes
#Anode
An_Free_Water = []
An_Total_Water=[]

#Cathode
Cat_Free_Water = []
Cat_Total_Water=[]


def cluster_loop(x,y,z):
    yes = [x,y,z]
    for i in range(x,y,z):
        if i == 0:
            i = i + 1
        elif i == 101:
            i = i + 1
        elif i == 202:
            i = i + 1
        elif i == 403:
            i = i + 1
        elif i == 604:
            i = i + 1

        ##    Lists the current frame   ##
        time = i
        u.trajectory[time]
        print("\n" + "Timestep:", time)
        print("These are the results for", u.trajectory[time], "\n")

        ##  Create empty arrays for the Free-H and Free-O atoms near the electrodes  ##
        #Anode
        An_Free_H_Array = []
        An_Free_O_Array = []

        #Cathode
        Cat_Free_H_Array = []
        Cat_Free_O_Array = []
        
        
        ##  Identifies the total O atoms (water) that are not coordinated to a Mg ion  ##
        #Group of O atoms that are within the cutoff of the anode
        An_O = u.select_atoms(f"type 4 and (around {ele_cutoff} group C_an)", C_an = C_an)
        n_An_O = len(An_O)
        print('There are',n_An_O, 'O atoms near the Anode')

        #Group of O atoms that are not within the cutoff of Mg and within the cutoff of the anode
        An_Free_O = u.select_atoms(f"type 4 and (around {ele_cutoff} group C_an) and (not around {cutoff} type 3)", C_an = C_an)
        n_An_Free_O = len(An_Free_O)
        print('There are',n_An_Free_O, 'free water O atoms near the Anode')

        #Group of O atoms that are within the cutoff of the anode
        Cat_O = u.select_atoms(f"type 4 and (around {ele_cutoff} group C_cat)", C_cat = C_cat)
        n_Cat_O = len(Cat_O)
        print('There are',n_Cat_O, 'O atoms near the Cathode')

        #Group of O atoms that are not within the cutoff of Mg and within the cutoff of the anode
        Cat_Free_O = u.select_atoms(f"type 4 and (around {ele_cutoff} group C_cat) and (not around {cutoff} type 3)", C_cat = C_cat)
        n_Cat_Free_O = len(Cat_Free_O)
        print('There are',n_Cat_Free_O, 'free water O atoms near the Cathode \n')

        #Places the indices of each free O-atom into one array (anode)
        for oa in An_Free_O:
            #print(o.index) #Uncomment for sanity check
            An_Free_O_Array.append(int(oa.index))
        #print(len(Free_O_Array),'\n')

        #Places the indices of each free O-atom into one array (cathode)
        for oc in Cat_Free_O:
            #print(o.index) #Uncomment for sanity check
            Cat_Free_O_Array.append(int(oc.index))
        #print(len(Free_O_Array),'\n')
        

        ##  Identifies the total H atoms (water) that are not coordinated to a Cl ion, and finds their attached O parent-atom  ##
        #Group of H atoms that are within the cutoff of the anode
        An_H = u.select_atoms(f"type 2 and (around {ele_cutoff} group C_an)", C_an = C_an)
        n_An_H = len(An_H)
        print('There are',n_An_H, 'H atoms near the Anode')

        #Group of H atoms that are not within the cutoff of Cl and within the cutoff of the anode
        An_Free_H = u.select_atoms(f"type 2 and (around {ele_cutoff} group C_an) and (not around {cutoff} type 1)", C_an = C_an)
        n_An_Free_H = len(An_Free_H)
        print('There are',n_An_Free_H, 'free H atoms around the Anode')

        #Group of H atoms that are within the cutoff of the cathode
        Cat_H = u.select_atoms(f"type 2 and (around {ele_cutoff} group C_cat)", C_cat = C_cat)
        n_Cat_H = len(Cat_H)
        print('There are',n_Cat_H, 'H atoms near the Cathode')

        #Group of H atoms that are not within the cutoff of Cl and within the cutoff of the cathode
        Cat_Free_H = u.select_atoms(f"type 2 and (around {ele_cutoff} group C_cat) and (not around {cutoff} type 1)", C_cat = C_cat)
        n_Cat_Free_H = len(Cat_Free_H)
        print('There are',n_Cat_Free_H, 'free H atoms around the Cathode')

        #Finds the O parent-atom for each 'free H' and places the parent O indices into an array (Anode)
        for ha in An_Free_H:
            #print(h.index, h.bonded_atoms.indices)
            An_Free_H_Array.append(int((ha.bonded_atoms.indices)))
            #print(u.atoms[h.bonded_atoms.indices])
        #print(len(Free_H_Array),'\n')

        #Finds the O parent-atom for each 'free H' and places the parent O indices into an array (Cathode)
        for hc in Cat_Free_H:
            #print(h.index, h.bonded_atoms.indices)
            Cat_Free_H_Array.append(int((hc.bonded_atoms.indices)))
            #print(u.atoms[h.bonded_atoms.indices])
        #print(len(Free_H_Array),'\n')


        ##  Uncomment for sanity check - shows the full list of O indices that are considered 'Free'  ##
        #print(Free_O_Array)
        #print(Free_H_Array, '\n')


        ##  Finds the shared O-indices and removes any duplicates to create the final free water list (Anode) ##
        #Find shared O-indices between both lists
        An_Reference_O = set(An_Free_O_Array)
        An_Filtered_H_Free_Water = [item for item in An_Free_H_Array if item in An_Reference_O]
        #print(Filtered_H_Free_Water) #Uncomment for sanity Check

        #Remove duplicates from the filtered list to get final free water count
        An_Final_Free_Water = list(dict.fromkeys(An_Filtered_H_Free_Water))
        #print(Final_Free_Water, '\n') #Uncomment for sanity check
        #print(n_An_O) #Uncomment for sanity check
        An_Total_Water.append(n_An_O)
        #print(len(An_Final_Free_Water)) #Uncomment for sanity check
        An_Free_Water.append(len(An_Final_Free_Water))


        ##  Finds the shared O-indices and removes any duplicates to create the final free water list (Cathode) ##
        #Find shared O-indices between both lists
        Cat_Reference_O = set(Cat_Free_O_Array)
        Cat_Filtered_H_Free_Water = [item for item in Cat_Free_H_Array if item in Cat_Reference_O]
        #print(Filtered_H_Free_Water) #Uncomment for sanity Check

        #Remove duplicates from the filtered list to get final free water count
        Cat_Final_Free_Water = list(dict.fromkeys(Cat_Filtered_H_Free_Water))
        #print(Final_Free_Water, '\n') #Uncomment for sanity check
        #print(n_Cat_O) #Uncomment for sanity check
        Cat_Total_Water.append(n_Cat_O)
        #print(len(Cat_Final_Free_Water)) #Uncomment for sanity check
        Cat_Free_Water.append(len(Cat_Final_Free_Water))
        

##  Run the free water analysis code and print the outputs  ##
cluster_loop(0,101,2)


#Prints the array for the anode interface free/total water
print(An_Total_Water)
print(An_Free_Water,'\n')

#Prints the array for the anode interface free/total water
print(Cat_Total_Water)
print(Cat_Free_Water,'\n')

An_data = {
    "Total O": An_Total_Water,
    "Total Free O": An_Free_Water
}

Cat_data = {
    "Total O": Cat_Total_Water,
    "Total Free O": Cat_Free_Water
}
      
## Use ExcelWriter to write the free water data into an excel sheet
with pd.ExcelWriter('PG Interface Free Water (Neutral).xlsx', engine='openpyxl') as writer:
   pd.DataFrame(An_data).to_excel(writer, sheet_name='Anode Water', index=False)
   pd.DataFrame(Cat_data).to_excel(writer, sheet_name='Cathode Water', index=False)
