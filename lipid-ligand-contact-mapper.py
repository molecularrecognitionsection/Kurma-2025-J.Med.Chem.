import MDAnalysis as mda
import os

u = mda.Universe('structure.psf', 'swag.dcd')
time_list = []
CHL1_list = []
POPC_list = []
POPE_list = []
POPS_list = []
POPI_list = []
PSM_list = []
conv = 2

for i, ts in enumerate(u.trajectory):
    CHL1 = len(list(u.select_atoms('resname CHL1 and same residue as around 4 resname LIG').residues))/174
    POPC = len(list(u.select_atoms('resname POPC and same residue as around 4 resname LIG').residues))/162
    POPE = len(list(u.select_atoms('resname POPE and same residue as around 4 resname LIG').residues))/102
    POPS = len(list(u.select_atoms('resname POPS and same residue as around 4 resname LIG').residues))/30
    POPI = len(list(u.select_atoms('resname POPI and same residue as around 4 resname LIG').residues))/30
    PSM = len(list(u.select_atoms('resname PSM and same residue as around 4 resname LIG').residues))/90

    time_list.append(i*conv)
    CHL1_list.append(CHL1)
    POPC_list.append(POPC)
    POPE_list.append(POPE)
    POPS_list.append(POPS)
    POPI_list.append(POPI)
    PSM_list.append(PSM)
#    print(water_num)

if os.path.exists('lipid-ligand.csv'):
    os.remove('lipid-ligand.csv')
with open('lipid-ligand.csv','a') as f:
    for t, n1, n2, n3, n4, n5, n6 in zip(time_list, CHL1_list, POPC_list, POPE_list, POPS_list, POPI_list, PSM_list):
        f.write(str(t) + ',' + str(n1) + ',' + str(n2) + ',' + str(n3) + ',' + str(n4) + ',' + str(n5) + ',' + str(n6) + '\n')
    
import numpy as np
from matplotlib import pyplot as plt

time = np.genfromtxt('lipid-ligand.csv', delimiter = ',', usecols=0)
C1 = np.genfromtxt('lipid-ligand.csv', delimiter = ',', usecols=1)
PC = np.genfromtxt('lipid-ligand.csv', delimiter = ',', usecols=2)
PE = np.genfromtxt('lipid-ligand.csv', delimiter = ',', usecols=3)
PS = np.genfromtxt('lipid-ligand.csv', delimiter = ',', usecols=4)
PI = np.genfromtxt('lipid-ligand.csv', delimiter = ',', usecols=5)
PM = np.genfromtxt('lipid-ligand.csv', delimiter = ',', usecols=6)
plt.figure()
plt.plot(time, C1, linewidth=0.5, color='green', label='CHL1')
plt.plot(time, PC, linewidth=0.5, color='tan', label='POPC')
plt.plot(time, PE, linewidth=0.5, color='magenta', label='POPE')
plt.plot(time, PS, linewidth=0.5, color='gold', label='POPS')
plt.plot(time, PI, linewidth=0.5, color='orange', label='POPI')
plt.plot(time, PM, linewidth=0.5, color='pink', label='PSM')
plt.ylabel('Lipid fraction')
plt.xlabel('Simulation time (ns)')
plt.xlim(0, 500)
plt.title('Ligand-lipid contacts', fontweight='bold')
plt.legend()
plt.savefig('lipid-ligand.png', dpi=300)

