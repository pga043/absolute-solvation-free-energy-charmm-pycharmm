import sys
import numpy as np
import pymbar
from pymbar import *
import seaborn as sns
import matplotlib.pyplot as plt 

temp = 298.15
beta = 1/(0.001987204259*temp)
def check_convergence(ligPath,throttle=1):
    """
    Given a ligPath where ff_e1.csv, ff_e2.csv, crn_e1.csv,
    crn_e2.csv live, 
    calculate free energies every `throttle` number of frames.
    """
    # Load energies
    ## Forward direction
    ffE1 = np.loadtxt(f'{ligPath}crn2ffprot/ff.ener',dtype=float, skiprows=1, usecols=1)
    crnE1 = np.loadtxt(f'{ligPath}/crn2ffprot/renormal.ener',dtype=float, skiprows=1, usecols=1)
 
    ## Backward direction
    crnE2 = np.loadtxt(f'{ligPath}/prot/renormal.ener',dtype=float, skiprows=1, usecols=1)
    ffE2 = np.loadtxt(f'{ligPath}/prot/ff.ener',dtype=float, skiprows=1, usecols=1)
    
    nframes = len(ffE1) # Assumes same nframes per sim
    
    frameList = []
    ddGList = []
    uncList = []
    overlapList = []
    for frame in range(nframes):
        if frame % throttle == 0 and frame != 0:
            # Capture frame index
            frameList.append(frame+1)
            
            # Assemble MBAR matrix up to frame `frame`
            ffE = np.concatenate((ffE1[:frame], ffE2[:frame]))
            crnE = np.concatenate((crnE1[:frame], crnE2[:frame]))
            

            u_kn = np.stack((ffE,crnE),axis=0)*beta
            N_k = np.array([len(ffE1[:frame]),len(crnE2[:frame])])

            # Instantiate MBAR object
            mbar = MBAR(u_kn,N_k)

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

lig1Path = sys.argv[1]

throttle = 10
sub = sys.argv[2]

frameList, ddGList, uncList, overlapList = check_convergence(lig1Path,throttle=throttle)
plt.plot(frameList,ddGList)
plt.errorbar(frameList, ddGList, yerr=uncList,fmt='o')
plt.title('Convergence Analysis for Ligand ' +str(sub) + ' in Protein',fontsize=15)
plt.xlabel('Frame Index',fontsize=12)
plt.ylabel('Free Energy (kcal/mol)',fontsize=12)
#plt.show()
plt.tight_layout()
plt.savefig('s1s' +str(sub) + '.pdf', dpi=300)

dG = ddGList[-1]
ddG = uncList[-1]

overlap = overlapList[-1]

with open('s1s' +str(sub) +'_mbar.dat','a') as f:
    sys.stdout = f
    print('MBAR free energy difference from FF -> CRN is {:.3f} +- {:.3f} kcal/mol'.format(dG, ddG)) # unit: kcal/mol
    print(f'Overlap between forward and reverse distributions is: {overlap*100:.2f}%')


quit()
