# import both Numpy and FastMBAR package
import sys
import numpy as np
from FastMBAR import *
import matplotlib.pyplot as plt

"""
test for pytorch and fastmbar installation
# python -c "import torch;  print(torch.cuda.is_available())" 
# python -m FastMBAR.test_installation

E = original potential energy (kcal/mol)
U = E/kB*T  => unitless

F = -ln integral{ exp (-U)}dx  => reduced free energy (unitless) : MBAR will print this

Free energy = kB*T*F

"""

kB=0.0019872041 # unit: kcal/(mol*K)
T=298           # unit: K


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

# construct and initialize a FastMBAR object with the energy matrix E and the vector v
# set cuda = True in the following command if you want to run the calcuation on GPUs
fastmbar = FastMBAR(energy = U, num_conf = v, cuda=False, verbose = True, bootstrap = True)
#print(fastmbar.__dict__.keys())
#print(fastmbar.num_states)
#print(fastmbar.num_conf)
#print(fastmbar.tot_num_conf)
#quit()

# after initialization, the relative free energies of the M states is stored in fastmbar.F
#print(fastmbar.F)
#print(fastmbar.F_std)

#------------  convert dG and uncertanity to kcal/mol -------------- #
dG_all = (fastmbar.F)*float(kB)*float(T)
ddG_all = fastmbar.F_std *float(kB)*float(T)
dG = ((fastmbar.F)[-1])*float(kB)*float(T)
ddG = (fastmbar.F_std[-1])*float(kB)*float(T)

print(dG_all)

x = np.arange(0, len(dG_all), 1)
plt.plot(dG_all, marker='o', label=f'FastMBAR = {round(dG, 2)} +/- {round(ddG ,2)} kcal/mol')
plt.errorbar(x, dG_all, yerr=ddG_all)
plt.xlabel('window')
plt.ylabel(r'$\Delta$$\Delta$G (kcal/mol)')
plt.legend(frameon=False)
plt.savefig('fastmbar_pmf.pdf', dpi=600)
plt.close()
with open('results.txt','w') as f:
    sys.stdout = f
    print('fastMBAR free energy difference is {:.3f} +- {:.3f} kcal/mol'.format(dG, ddG))


with open('junk_fastmbar.dat', 'a') as f1:
     for i in range(len(dG_all)):
        f1.write(f"{round(fastmbar.F[i]*float(kB)*float(T), 3)} +- {round(fastmbar.F_std[i]*float(kB)*float(T), 3)} kcal/mol \n")


quit()


