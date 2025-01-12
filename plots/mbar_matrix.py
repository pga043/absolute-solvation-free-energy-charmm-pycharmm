# import both Numpy and FastMBAR package
import numpy as np
import sys
from pymbar import *
import matplotlib.pyplot as plt
import seaborn as sns
"""
E = original potential energy (kcal/mol)
U = E/kB*T  => unitless

F = -ln integral{ exp (-U)}dx  => reduced free energy (unitless) : MBAR will print this

Free energy = kB*T*F

"""

nwins  = 25
frames = 250
anadir = 'post'

kB=0.0019872041 # unit: kcal/(mol*K)
T=298.15           # unit: K


# Initialize an empty list to store rows dynamically
matrix = []

for i in range(1, nwins+1):
    win = []
    for j in range(1, nwins+1):
        ener = np.genfromtxt(f'win{j}/{anadir}/win{i}.dat', skip_header=1, usecols=-1, dtype=np.float128, delimiter=' ')
        win.extend(ener)

    matrix.append(win)

# construct the energy matrix A and the vector v
E = np.array(matrix)
v = np.full(int(nwins), int(frames))

U = E / (float(kB) * float(T))

#mbar = MBAR(u_kn = U, N_k = v, initialize='zeros', solver_protocol='robust', verbose=True) # works
mbar = MBAR(u_kn = U, N_k = v, initialize='BAR', verbose=True) # works

results = mbar.compute_free_energy_differences()
S_H = mbar.compute_entropy_and_enthalpy()

#print('MBAR free energy difference is {:.3f} +- {:.3f} kT'.format(results['Delta_f'][0, -1], results['dDelta_f'][0, -1])) # unit: kT
overlapMatrix = mbar.compute_overlap()['matrix']
applyall = np.vectorize(lambda x: x*100)
overlapMatrix = applyall(overlapMatrix)
applyall = np.vectorize(int)
overlapMatrix = applyall(overlapMatrix)
# plt.imshow(overlapMatrix,cmap='PiYG',interpolation='nearest')
sns.heatmap(overlapMatrix, linewidth=0.5,annot=True)
#plt.savefig('overlap_matrix.pdf', dpi=1000)
#plt.close()
plt.show()
quit()

#print(results['dDelta_f']) # unit: kT
#print(S_H['Delta_s'])   # Entropy
#print(S_H['Delta_u'])   # Enthalpy

# ---------- converting free energy to kcal/mol ----------------- #
dG_all = results['Delta_f'] * (float(kB) * float(T)) # unit: kcal/mol
dG = results['Delta_f'][0, -1] * (float(kB) * float(T)) # unit: kcal/mol
ddG = (results['dDelta_f'][0, -1]) * (float(kB) * float(T)) # unit: kcal/mol
plt.plot(dG_all[0])
#plt.savefig('mbar_pmf.pdf')
#plt.close()
plt.show()
#print('MBAR free energy difference is {:.3f} +- {:.3f} kcal/mol'.format(dG, ddG)) # unit: kcal/mol
#with open('results.txt','a') as f:
#    sys.stdout = f
#    print('MBAR free energy difference is {:.3f} +- {:.3f} kcal/mol'.format(dG, ddG)) # unit: kcal/mol

quit()

#------------------------ save raw data in txt files -------------------------#
#np.savetxt('1test.txt', dG, delimiter=' ')
#np.savetxt('3test.txt', np.transpose(dG), delimiter=' ')

for i in range(38):
    np.savetxt(str(i)+ 'test.txt', dG_all[i], delimiter=' ')


