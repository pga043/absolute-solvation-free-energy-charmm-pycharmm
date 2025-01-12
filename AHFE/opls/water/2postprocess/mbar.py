import sys
import numpy as np
import pymbar
from pymbar import *

# conda activate openmm

"""
E = original potential energy (kcal/mol)
U = E/kB*T  => unitless

F = -ln integral{ exp (-U)}dx  => reduced free energy (unitless) : MBAR will print this

Free energy = kB*T*F

"""

kB=0.0019872041 # unit: kcal/(mol*K)
T=298           # unit: K


# construct the energy matrix A and the vector v
frwd = np.genfromtxt(sys.argv[1], dtype=None, skip_header=0, delimiter=' ')
bcwd = np.genfromtxt(sys.argv[2], dtype=None, skip_header=0, delimiter=' ')
#print(E.shape)
#print(v.shape)
#quit()
Uf = frwd / (float(kB) * float(T)) 
Ub = bcwd / (float(kB) * float(T))

#print(U)
#print(U.shape)
#quit()

#---------------- Bennett acceptance ratio ------------------------ #
bar = pymbar.other_estimators.bar(Uf, Ub, compute_uncertainty=True, maximum_iterations=1000)
dG = bar['Delta_f'] * (float(kB) * float(T)) # kcal/mol
ddG = bar['dDelta_f'] * (float(kB) * float(T))

#print(bar['Delta_f'])
#print(bar['dDelta_f'])
print('BAR free energy difference is {:.3f} +- {:.3f} kT'.format(bar['Delta_f'], bar['dDelta_f']))
print('BAR free energy difference is {:.3f} +- {:.3f} kcal/mol'.format(dG, ddG))

#quit()
#--------------- one-sided (unidirectional) exponential averaging (EXP) ------------- #
forward = pymbar.other_estimators.exp(Uf, compute_uncertainty=True, is_timeseries=False)
dG_f = forward['Delta_f'] * (float(kB) * float(T)) # kcal/mol
ddG_f = forward['dDelta_f'] * (float(kB) * float(T))

print('Forward free energy difference is {:.3f} +- {:.3f} kT'.format(forward['Delta_f'], forward['dDelta_f']))
print('Forward free energy difference is {:.3f} +- {:.3f} kcal/mol'.format(dG_f, ddG_f))


backward = exp(Ub)
dG_b = backward['Delta_f'] * (float(kB) * float(T))  # kcal/mol
ddG_b = backward['dDelta_f'] * (float(kB) * float(T))
print('Reverse free energy difference is {:.3f} +- {:.3f} kT'.format(backward['Delta_f'], backward['dDelta_f']))
print('Reverse free energy difference is {:.3f} +- {:.3f} kcal/mol'.format(dG_b, ddG_b))

# ---------------- Compute overlap between forward and backward ensembles (using MBAR definition of overlap) ---------- #
# overlap – The overlap: 0 denotes no overlap, 1 denotes complete overlap
print('Definition:: The overlap: 0 denotes no overlap, 1 denotes complete overlap')
print('The overlap between'+' '+str(sys.argv[1])+' and '+str(sys.argv[2])+' '+'is {:.6f}'.format(pymbar.other_estimators.bar_overlap(Uf, Ub)))



#------------ gaussian approximation to one-sided (unidirectional) exponential averaging ----------- #
#onward = exp_gauss(Uf, compute_uncertainty=True, is_timeseries=True)
#print('Forward Gaussian approximated free energy difference is {:.3f} +- {:.3f} kT'.format(onward['Delta_f'], onward['dDelta_f']))
#reverse = exp_gauss(Ub)
#print('Reverse Gaussian approximated free energy difference is {:.3f} +- {:.3f} kT'.format(reverse['Delta_f'], reverse['dDelta_f']))

