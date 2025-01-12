import sys
import numpy as np

"""
check here for the equation: https://en.wikipedia.org/wiki/Free_energy_perturbation

For going in forward direction,
dE =  E_i+1 - E_i ; both energies calculated from trajectory of window i using lambdas i and i+i

For going in backward direction,
dE =  E_i-1 - E_i ; both energies calculated from trajectory of window i using lambdas i and i-i

"""

kB=0.0019872041 # unit: kcal/(mol*K)
T=298           # unit: K

# ------------ forward direction --------------------
frwd = np.genfromtxt(sys.argv[1], dtype=None, skip_header=0, delimiter=' ')

exp = np.exp(-frwd/(float(kB)*float(T)))
#print(exp[:10])

average = np.average(exp)
#print(average)
#quit()

deltaF_frwd = -float(kB)*float(T) * np.log(average)

#print(deltaF_frwd)

# ----------- backward direction ----------------
bcwd = np.genfromtxt(sys.argv[2], dtype=None, skip_header=0, delimiter=' ')

exp = np.exp(-bcwd/(float(kB)*float(T)))

average = np.average(exp)

deltaF_bcwd = -float(kB)*float(T) * np.log(average)

print('Forward free energy difference is {:.3f} kcal/mol'.format(deltaF_frwd))
print('Reverse free energy difference is {:.3f} kcal/mol'.format(deltaF_bcwd))

