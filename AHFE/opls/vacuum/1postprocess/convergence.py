import sys
import numpy as np
import pymbar
from pymbar import *
import seaborn as sns
import matplotlib.pyplot as plt 

T = 298.15
kB=0.0019872041 # unit: kcal/(mol*K)
temp = 298.15
beta = 1/(0.001987204259*temp)

nwins = 20
nframes = 100

E = np.genfromtxt('matrix.dat', dtype=None, skip_header=0, delimiter=' ')
U = E / (float(kB) * float(T))

'''
Shape of MBAR matrix energies
E{win1_traj}_lam1  E{win1_traj}_lam2  E{win1_traj}_lam3 ... E{win1_traj}_lamN

E{win2_traj}_lam1  E{win2_traj}_lam2  E{win2_traj}_lam3 ... E{win2_traj}_lamN
.
.
.

E{winN_traj}_lam1  E{winN_traj}_lam2  E{winN_traj}_lam3 ... E{winN_traj}_lamN
'''

fig, ax = plt.subplots(2, 10, figsize=(8,8))
for i in range(nwins):
    j = i + 1
    if j < nwins:
       start_i = i * nframes 
       end_i   = j * nframes
       end_j   = (j + 1) * nframes
       U_f1 = U[i][start_i:end_i]
       U_f2 = U[i][end_i:end_j]

       U_r2 = U[j][end_i:end_j]
       U_r1 = U[j][start_i:end_i]

#forward = np.concatenate((U_f1, U_r1))
#reverse = np.concatenate((U_f2, U_r2))

#u_kn = np.stack((forward, reverse),axis=0)
#print(u_kn[1])
#print(u_kn.shape)
#quit()

#------------------------------------------------------------
       frameList = []
       ddGList = []
       uncList = []
       overlapList = []
       throttle = 10

       for frame in range(nframes):
           if frame % throttle == 0 and frame != 0:
               # Capture frame index
               frameList.append(frame+1)

               # Assemble MBAR matrix up to frame `frame`
               forward = np.concatenate((U_f1[:frame], U_f2[:frame])) ## resultant matrix shape: 1, 2*nframes
               reverse = np.concatenate((U_r1[:frame], U_r2[:frame])) ## resultant matrix shape: 1, 2*nframes

               u_kn = np.stack((forward, reverse),axis=0) ## output is a 2-D array. They are stacked row-wise.
               N_k = np.array([len(U_f1[:frame]),len(U_r2[:frame])])

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
       if i < int(nwins / 2):
          ax[0, i].plot(frameList,ddGList, label=f'win_{i+1}_{j+1}')
          ax[0, i].errorbar(frameList, ddGList, yerr=uncList,fmt='o') 
          ax[0, i].legend(frameon=False)
       else:
           p = i - int(nwins / 2)
           ax[1, p].plot(frameList,ddGList, label=f'win_{i+1}_{j+1}')
           ax[1, p].errorbar(frameList, ddGList, yerr=uncList,fmt='o')
           ax[1, p].legend(frameon=False)

#       plt.plot(frameList,ddGList, label=f'win_{i+1}_{j+1}')
#       plt.errorbar(frameList, ddGList, yerr=uncList,fmt='o') 
#       plt.title('Convergence Analysis ',fontsize=15)
#       plt.xlabel('Frame Index',fontsize=12)
#       plt.ylabel('Free Energy (kcal/mol)',fontsize=12)
#       plt.show()
#       plt.tight_layout()
#       plt.savefig('s1s' +str(sub) + '.pdf', dpi=300)

       dG = ddGList[-1]
       ddG = uncList[-1]

       overlap = overlapList[-1]

#with open('s1s' +str(sub) +'_mbar.dat','a') as f:
#    sys.stdout = f
       print('MBAR free energy difference is {:.3f} +- {:.3f} kcal/mol'.format(dG, ddG)) # unit: kcal/mol
       print(f'Overlap between forward and reverse distributions is: {overlap*100:.2f}%')



ax[1,0].set_xlabel('Frame Index',fontsize=12)
ax[1,0].set_ylabel('Free Energy (kcal/mol)',fontsize=12)
ax[0,0].set_ylabel('Free Energy (kcal/mol)',fontsize=12)

# Set overall title and tight layout
fig.suptitle('Convergence Analysis ',fontsize=15)

fig.tight_layout()

# Show the plot
plt.show()

quit()
