import numpy as np
from numpy.linalg import norm
import MDAnalysis as mda
from matplotlib import pyplot as plt
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

u1 = mda.Universe("structure.psf", "swag.dcd")
#u2 = mda.Universe("MD1/dry.pdb", "MD2/dry.dcd")
#u3 = mda.Universe("MD1/dry.pdb", "MD3/dry.dcd")
#u4 = mda.Universe("MD1/dry.pdb", "MD4/dry.dcd")
#u5 = mda.Universe("MD1/dry.pdb", "MD5/dry.dcd")

def calc_theta(ag):
    """Calculate the TM7-H8 angle for A3R in degrees
    
    Parameters
    ----------
    ag : MDAnalysis.Universe or MDAnalysis.AtomGroup
    """
    TM7_start = ag.select_atoms("segid PROA and resid 276 and backbone").center_of_geometry()
    TM7_end = ag.select_atoms("segid PROA and resid 285 and backbone").center_of_geometry()
    TM7_vector = TM7_end - TM7_start
    H8_start = ag.select_atoms("segid PROA and resid 286 and backbone").center_of_geometry()
    H8_end = ag.select_atoms("segid PROA and resid 303 and backbone").center_of_geometry()
    H8_vector = H8_end - H8_start
    theta = np.arccos(np.dot(TM7_vector, H8_vector)/(norm(TM7_vector)*norm(H8_vector)))
    return np.rad2deg(theta)
    
def calc_phi(ag):
    """Calculate the TM1-H8 angle for A3R in degrees
    
    Parameters
    ----------
    ag : MDAnalysis.Universe or MDAnalysis.AtomGroup
    """
    TM1_start = ag.select_atoms("segid PROA and resid 25 and backbone").center_of_geometry()
    TM1_end = ag.select_atoms("segid PROA and resid 40 and backbone").center_of_geometry()
    TM1_vector = TM1_end - TM1_start
    H8_start = ag.select_atoms("segid PROA and resid 286 and backbone").center_of_geometry()
    H8_end = ag.select_atoms("segid PROA and resid 303 and backbone").center_of_geometry()
    H8_vector = H8_end - H8_start
    phi = np.arccos(np.dot(TM1_vector, H8_vector)/(norm(TM1_vector)*norm(H8_vector)))
    return np.rad2deg(phi)    

times = np.empty(u1.trajectory.n_frames)
theta = np.empty(u1.trajectory.n_frames)
phi = np.empty(u1.trajectory.n_frames)

# iterate through the trajectory and gather data
for i, ts in enumerate(u1.trajectory):
    times[i] = i * 2
    theta[i] = calc_theta(u1)
    phi[i] = calc_phi(u1)
    

# import pandas as pd
import pandas as pd
 
df = pd.DataFrame(list(zip(times, theta, phi)),
               columns =['time', 'theta', 'phi'])
df.to_csv('MD-angles.csv')

    

# plot
ax1 = plt.subplot(121)
ax1.scatter(times, theta, label="TM7-H8 angle", s=2)
ax1.scatter(times, phi, label="TM1-H8 angle", s=2)
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("Angle (degrees)")
ax1.legend(loc="best")
ax2 = plt.subplot(122)
ax2.scatter(theta, phi, s=2)
ax2.set_xlabel("TM7-H8 angle")
ax2.set_ylabel("TM1-H8 angle")
ax2.set_aspect(1)
plt.tight_layout()
plt.savefig('MD-angle.png')
#print(np.mean(theta[0:200]), np.mean(theta[1800:-1]))


