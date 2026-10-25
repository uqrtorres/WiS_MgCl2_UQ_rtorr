##   Import MDAnalysis, numpy and others   ##
import MDAnalysis as mda
import numpy as np
from MDAnalysis import analysis
from MDAnalysis.analysis import lineardensity
from MDAnalysis.analysis.lineardensity import LinearDensity
import pandas as pd

import warnings
# suppress some MDAnalysis warnings about PSF files
warnings.filterwarnings('ignore')

from matplotlib import pyplot as plt


##    Load the data/dump file of interest and lists the total atoms in the system   ##
u = mda.Universe('3rd_30ns_Equilibriate_PG_6.3M_MgCl2_nvt.data',
                 'Dump_3rd_2.5V_30ns_PG_MgCl2.lammpstrj',
                 topology_format="data", format="lammpsdump",
                 atom_style = 'id resid type charge x y z')

##   Select the solution atoms for density profiling  ##
solution = u.select_atoms("not type 5")
print(len(solution), "atoms within the solution")

#C_wall = u.select_atoms("type 5")
#print(len(C_wall), "C atoms within the solution")

Mg_density = u.select_atoms("type 3")
print(len(Mg_density), "Mg atoms within the system")

Cl_density = u.select_atoms("type 1")
print(len(Cl_density), "Cl atoms within the system")

H_density = u.select_atoms("type 2")
print(len(H_density), "H atoms within the system")

O_density = u.select_atoms("type 4")
print(len(O_density), "O atoms within the system")

#H2O_density = u.select_atoms("type 2 or type 4 ")
#print(len(H2O_density), "water atoms within the system")

##  Calculate the linear density of the solution in the system  ##
ldens = LinearDensity(solution, grouping='atoms', binsize=0.25)
ldens.run()

#C_ldens = LinearDensity(C_wall, grouping='atoms', binsize=1)
#C_ldens.run()

Mg_ldens = LinearDensity(Mg_density, grouping='atoms', binsize=0.25)
Mg_ldens.run()

Cl_ldens = LinearDensity(Cl_density, grouping='atoms', binsize=0.25)
Cl_ldens.run()

H_ldens = LinearDensity(H_density, grouping='atoms', binsize=0.25)
H_ldens.run()

O_ldens = LinearDensity(O_density, grouping='atoms', binsize=0.25)
O_ldens.run()

#H2O_ldens = LinearDensity(H2O_density, grouping='atoms', binsize=0.25)
#H2O_ldens.run()

# Get box Z length from universe dimensions
z_length = u.dimensions[2]
print(z_length)

# Obtain the bins created for the z-axis
n_bins = len(ldens.results.z.mass_density)
print(n_bins)

##  Plot the results  ##
# ldens.results.z['pos'] contains the mass-weighted density 
z_coordinates = np.linspace(0, z_length, n_bins)
density_profile = ldens.results.z.mass_density
#C_density_profile = C_ldens.results.z['pos']
Mg_density_profile = Mg_ldens.results.z.mass_density
Cl_density_profile = Cl_ldens.results.z.mass_density
H_density_profile = H_ldens.results.z.mass_density
O_density_profile = O_ldens.results.z.mass_density
#H2O_density_profile = H2O_ldens.results.z.mass_density


# Creates the visual density profile graph
#plt.plot(z_coordinates, density_profile, color="black", label="Solution Density")
plt.plot(z_coordinates, Mg_density_profile, color="green", label="Mg Density")
plt.plot(z_coordinates, Cl_density_profile, color="pink", label="Cl Density")
plt.plot(z_coordinates, H_density_profile, color="blue", label="H Density")
plt.plot(z_coordinates, O_density_profile, color="red", label="O Density")
#plt.plot(z_coordinates, H2O_density_profile, color="blue", label="H2O Density")
#plt.plot(z_coordinates, C_density_profile, color="red", label="C Density")
plt.xlabel('Z-axis (Å)')
plt.ylabel('Density ($g / cm^3$)')
plt.title('Z-axis Linear Density Profile for 2.5V PG System')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.xlim(0,46.5)
plt.savefig('2.5V_Density_Profile.png', dpi=600, bbox_inches='tight')
plt.show()

# Creates an excel spreadsheet of the density profiles
df = pd.DataFrame({
    'Z-Axis ($10^-23 * cm^3$)': z_coordinates,
    'Mg': Mg_density_profile,
    'Cl': Cl_density_profile,
    'H': H_density_profile,
    'O': O_density_profile,
    #'H2O Density (g/cm^3)': H2O_density_profile,
    #'Overall Solution Density (g/cm^3)': density_profile 
})
print(df)
df.to_excel("2.5V Density Profiles.xlsx", sheet_name="3rd Run")

