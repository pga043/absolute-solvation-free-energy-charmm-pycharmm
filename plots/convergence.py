import sys
import numpy as np
import pymbar
from pymbar import *
import seaborn as sns
import matplotlib.pyplot as plt 

temp = 298.15
beta = 1/(0.001987204259*temp)
def check_convergence(i,j,anadir,throttle=1):
    """
    calculate free energies every `throttle` number of frames.
    """
    # Load energies
    w0_0 = np.loadtxt(f'win{i}/{anadir}/win{i}.dat',dtype=float, skiprows=1, usecols=-1)
    w0_1 = np.loadtxt(f'win{i}/{anadir}/win{j}.dat',dtype=float, skiprows=1, usecols=-1)

    w1_0 = np.loadtxt(f'win{j}/{anadir}/win{i}.dat',dtype=float, skiprows=1, usecols=-1)
    w1_1 = np.loadtxt(f'win{j}/{anadir}/win{j}.dat',dtype=float, skiprows=1, usecols=-1)
 
    nframes = len(w0_0) # Assumes same nframes per sim
    
    frameList = []
    ddGList = []
    uncList = []
    overlapList = []
    for frame in range(nframes):
        if frame % throttle == 0 and frame != 0:
            # Capture frame index
            frameList.append(frame+1)
            
            # Assemble MBAR matrix up to frame `frame`
            w0 = np.concatenate((w0_0[:frame], w1_0[:frame]))
            w1 = np.concatenate((w0_1[:frame], w1_1[:frame]))
            

            u_kn = np.stack((w0,w1),axis=0)*beta
            N_k = np.array([len(w0_0[:frame]),len(w1_1[:frame])])

            # Instantiate MBAR object
            mbar = MBAR(u_kn,N_k, initialize='BAR', verbose=True)

            result = mbar.compute_free_energy_differences()

            # Get free energy/uncertainty in kcal/mol
            fe = result['Delta_f'][0,1]/beta
            unc = result['dDelta_f'][0,1]/beta

            # Get overlap
            overlap = mbar.compute_overlap()
            overlap = overlap['matrix'][0,0]
            
            ddGList.append(fe)
            uncList.append(unc)
            overlapList.append(overlap)
            
    return frameList, ddGList, uncList, overlapList

#==================================================================
throttle = 10
nwins = 25
anadir = 'post'

for i in range(1, nwins):
    j = i + 1
    print(f'win{i} -> win{j}')
    frameList, ddGList, uncList, overlapList = check_convergence(i,j,anadir,throttle=throttle)
    dG = ddGList[-1]
    ddG = uncList[-1]
    print('MBAR free energy difference is {:.3f} +- {:.3f} kcal/mol'.format(dG, ddG)) # unit: kcal/mol
    plt.plot(frameList,ddGList, label='MBAR free energy difference is {:.3f} +- {:.3f} kcal/mol'.format(dG, ddG))
    plt.errorbar(frameList, ddGList, yerr=uncList,fmt='o')
    plt.title(f'Convergence Analysis for win{i} -> win{j}',fontsize=15)
    plt.xlabel('Frame Index',fontsize=12)
    plt.ylabel('Free Energy (kcal/mol)',fontsize=12)
    plt.legend()
    plt.tight_layout()
    #plt.savefig(f's1s{sub1}.s2s{sub2}.pdf', dpi=300)
    plt.show()

quit()
#overlap = overlapList[-1]

#with open(f's1s{sub1}.s2s{sub2}_mbar.dat','a') as f:
#    sys.stdout = f
#    print('MBAR free energy difference from FF -> CRN is {:.3f} +- {:.3f} kcal/mol'.format(dG, ddG)) # unit: kcal/mol
#    print(f'Overlap between forward and reverse distributions is: {overlap*100:.2f}%')


