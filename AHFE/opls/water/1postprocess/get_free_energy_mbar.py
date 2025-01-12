import numpy as np
from FastMBAR import *
from pymbar import *
import matplotlib.pyplot as plt
import sys

nwins  = 15
frames = 250
anadir = 'post'

# constants
kB=0.0019872041 # unit: kcal/(mol*K)
T=298           # unit: K

'''
shape of matrix
Data from -> win1     win2     win3     ........ winN
dat files -> win1.dat win1.dat win1.dat ........ win1.dat
             win2.dat win2.dat win2.dat ........ win2.dat
             win3.dat win3.dat win3.dat ........ win3.dat
             .
             .
             winN.dat winN.dat winN.dat ........ winN.dat

The actual simulations were run only for diagonal elements.
All the other matrix elements were obtained via postprocessing of each window.
'''

# Initialize an empty list to store rows dynamically
matrix = []

for i in range(1, nwins+1):
    win = []
    for j in range(1, nwins+1):
        ener = np.genfromtxt(f'../win{j}/{anadir}/win{i}.dat', skip_header=1, usecols=-1, dtype=np.float128, delimiter=' ')
        win.extend(ener)
    
    matrix.append(win)

# construct the energy matrix A and the vector v
E = np.array(matrix) 
v = np.full(int(nwins), int(frames)) 

U = E / (float(kB) * float(T))

#===============================================================================================
#================================================================================================
def fastmbar_free_energy(U, v):
    # construct and initialize a FastMBAR object with the energy matrix E and the vector v
    # set cuda = True in the following command if you want to run the calcuation on GPUs
    fastmbar = FastMBAR(energy = U, num_conf = v, cuda=False, verbose = True, bootstrap = True)

    #------------  convert dG and uncertanity to kcal/mol -------------- #
    dG_all = (fastmbar.F)*float(kB)*float(T)
    ddG_all = fastmbar.F_std *float(kB)*float(T)
    dG = ((fastmbar.F)[-1])*float(kB)*float(T)
    ddG = (fastmbar.F_std[-1])*float(kB)*float(T)

    plt.plot(dG_all, marker='o', label=f'FastMBAR = {round(dG, 2)} +/- {round(ddG ,2)} kcal/mol')
    x = np.arange(0, len(v), 1)
    plt.errorbar(x, dG_all, yerr=abs(ddG_all))

    plt.xlabel('window')
    plt.ylabel(r'$\Delta$G (kcal/mol)')
    plt.legend(frameon=False)
    #plt.savefig('fastmbar_pmf.png', dpi=600)
    plt.show()
    plt.close()
    with open('results.txt','a') as f1:
        f1.write('FastMBAR free energy difference is {:.3f} +- {:.3f} kcal/mol \n'.format(dG, ddG)) # unit: kcal/mol
        f1.close()


def mbar_free_energy(U, v):
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
    plt.ylabel(r'$\Delta$G (kcal/mol)')
    plt.legend(frameon=False)
    #plt.savefig('mbar_pmf.png', dpi=600)
    plt.show()
    plt.close()
    with open('results.txt','a') as f1:
         f1.write('MBAR free energy difference is {:.3f} +- {:.3f} kcal/mol \n'.format(dG, ddG)) # unit: kcal/mol

#=========================================================================================================
#==========================================================================================================
fastmbar_free_energy(U, v)
mbar_free_energy(U, v)
quit()


 
