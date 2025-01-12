# import both Numpy and FastMBAR package
import numpy as np
import sys
from pymbar import *
import matplotlib.pyplot as plt
"""
E = original potential energy (kcal/mol)
U = E/kB*T  => unitless

F = -ln integral{ exp (-U)}dx  => reduced free energy (unitless) : MBAR will print this

Free energy = kB*T*F

"""

kB=0.0019872041 # unit: kcal/(mol*K)
T=298.15           # unit: K


# construct the energy matrix A and the vector v
E = np.genfromtxt('matrix.dat', dtype=None, skip_header=0, delimiter=' ')
v = np.genfromtxt('vector.dat', dtype=None, skip_header=0, delimiter=' ')
#print(E.shape)
#print(v.shape)
#quit()
U = E / (float(kB) * float(T)) 

#print(U)
#print(U.shape)
#quit()

#mbar = MBAR(u_kn = U, N_k = v, initialize='zeros', solver_protocol='robust', verbose=True) # works
mbar = MBAR(u_kn = U, N_k = v, initialize='BAR', verbose=True) # works

results = mbar.compute_free_energy_differences()
S_H = mbar.compute_entropy_and_enthalpy()

print('MBAR free energy difference is {:.3f} +- {:.3f} kT'.format(results['Delta_f'][0, -1], results['dDelta_f'][0, -1])) # unit: kT

#print(results['dDelta_f']) # unit: kT
#print(S_H['Delta_s'])   # Entropy
#print(S_H['Delta_u'])   # Enthalpy

# ---------- converting free energy to kcal/mol ----------------- #
dG_all = results['Delta_f'] * (float(kB) * float(T)) # unit: kcal/mol
dG = results['Delta_f'][0, -1] * (float(kB) * float(T)) # unit: kcal/mol
ddG = (results['dDelta_f'][0, -1]) * (float(kB) * float(T)) # unit: kcal/mol

plt.plot(dG_all[0], marker='o', label=f'MBAR = {round(dG, 2)} +/- {round(ddG ,2)} kcal/mol')
x = np.arange(0, len(v), 1)
plt.errorbar(x, dG_all[0], yerr=abs(ddG))

plt.xlabel('window')
plt.ylabel(r'$\Delta$$\Delta$G (kcal/mol)')
plt.legend(frameon=False)
plt.savefig('mbar_pmf.pdf', dpi=600)
plt.close()
with open('results.txt','a') as f:
    sys.stdout = f
    print('MBAR free energy difference is {:.3f} +- {:.3f} kcal/mol'.format(dG, ddG)) # unit: kcal/mol

#------------------------ save raw data in txt files -------------------------#
#np.savetxt('1test.txt', dG, delimiter=' ')
#np.savetxt('3test.txt', np.transpose(dG), delimiter=' ')

with open('junk_mbar.dat', 'a') as f1:
     for i in range(len(dG_all[0])):
        f1.write(f"{round(results['Delta_f'][0, i] * float(kB) * float(T), 3)} +- {round(results['dDelta_f'][0, i] * float(kB) * float(T), 3)} kcal/mol \n")
quit()


