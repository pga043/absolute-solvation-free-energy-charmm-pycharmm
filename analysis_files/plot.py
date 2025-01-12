import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats

mol_list = np.genfromtxt('mol_list.txt', dtype=str)
msld =  np.genfromtxt('Result.txt', dtype=str)
cc =  np.genfromtxt('book-ending/8mbar/cc.dat', dtype=str)
expt = [] # list of experimental data

# compounds and expt. deltaG
compounds = {}
compounds['expt'] = expt
compounds['mol'] = mol_list[:,0].tolist()
compounds['cc'] = [float(x) for x in cc[:,0]]
compounds['cc_unc'] = [float(x) for x in cc[:,2]]
compounds['msld'] = []
compounds['msld_unc'] = []

for i in range(len(mol_list)):
    condition = (msld[:, 0] == mol_list[i, 1]) & (msld[:, 1] == mol_list[i, 2])
    msld_mols = msld[condition]
    if msld_mols.size > 0:
        for j in msld_mols:
            compounds['msld'].append(float(j[2]))
            compounds['msld_unc'].append(float(j[-1]))

# compounds and expt. deltaG
#compounds = {'mol': [11, 49, 51, 61, 64, 66, 67, 69, 97, 126],
#            'expt': [-7.87, -8.61, -8.65, -8.86, -7.43, -9.27, -8.39, -5.03, -10.05, -11.82],
#            'msld': [0.627, 0.00, 0.010, 0.315, -0.473, 0.191, -0.488, -0.167, -0.443, -0.662],
#            'cc': [-0.022, 0.000, 0.048, -0.005, 0.015, 0.039, 0.047, 0.011, 0.127, -0.038]}

compounds['msld_cc'] = []
compounds['msld_cc_unc'] = []
for i in range(len(compounds['msld'])):
    #print(compounds['msld'][i])
    compounds['msld_cc'].append(compounds['msld'][i] + compounds['cc'][i])
    compounds['msld_cc_unc'].append(math.sqrt((compounds['msld_unc'][i])**2 + (compounds['cc_unc'][i])**2))

average_expt = sum([x for x in compounds['expt']]) / len(compounds['expt'])
average_msld = sum([x for x in compounds['msld_cc']]) / len(compounds['msld_cc'])

compounds['abs_msld']  = []
for i in compounds['msld_cc']:
    #print(i)
    absolute = i + average_expt - average_msld
    compounds['abs_msld'].append(round(absolute, 1))


RMSE = math.sqrt(np.square(np.subtract(compounds['expt'], compounds['abs_msld'])).mean())
pearson = scipy.stats.pearsonr(compounds['expt'], compounds['abs_msld'])[0]
spearman = scipy.stats.spearmanr(compounds['expt'], compounds['abs_msld'])[0]


print([round(x, 1) for x in compounds['expt']])
print([round(x, 1) for x in compounds['abs_msld']],  [round(x, 1) for x in compounds['msld_cc_unc']])
#-----------------------------
print(f'RMSE = {round(RMSE, 2)}')
print(f'pearson = {round(pearson, 2)}')
print(f'spearman = {round(spearman, 2)}')

#------------------------------------------------
fig = plt.figure()

plt.plot(compounds['expt'], compounds['abs_msld'], linestyle='', marker='p', color='b', label=(r"charmm_opls (RMSE = {:.2f}, $\rho$ = {:.2f})".format(RMSE, pearson)))
plt.errorbar(compounds['expt'], compounds['abs_msld'], yerr=compounds['msld_cc_unc'], fmt='', linestyle='', color='red', alpha=0.4)

plt.plot(compounds['expt'], compounds['expt'], linestyle='-', color='k')
plt.plot(compounds['expt'], [x+1 for x in compounds['expt']], linestyle=':', color='r')
plt.plot(compounds['expt'], [x-1 for x in compounds['expt']], linestyle=':', color='r')

fig.suptitle('Human Neutrophil Elastase', fontsize=14)
plt.xlabel('$\Delta$$G_{expt}$ (kcal/mol)', fontsize=10)
plt.ylabel(r'$\Delta$$G_{MS\lambda D}^{FF}$ (kcal/mol)', fontsize=10)
#plt.xlim([-14,-9])
#plt.ylim([-16,-8])
plt.legend(fontsize="12", frameon=True, loc='upper left', ncol=2)
#--------------------------------------
# changing the fontsize of ticks
plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16)

# Customize the linewidth of the x and y axes
ax = plt.gca()
ax.spines['bottom'].set_linewidth(2)  # x-axis
ax.spines['left'].set_linewidth(2)    # y-axis
ax.spines['top'].set_linewidth(2)  # x-axis
ax.spines['right'].set_linewidth(2)    # y-axis

# Customize the tick sizes
ax.tick_params(axis='both', which='major', labelsize=16, width=2, length=8)  # Major ticks
#ax.tick_params(axis='both', which='minor', width=1, length=4)  # Minor ticks
#======================

plt.tight_layout()

plt.savefig('hne.png', dpi=600)
#plt.savefig('hne.pdf', dpi=600)
#plt.savefig('hne.svg', dpi=600)
#plt.savefig('Figure3.tif', dpi=300)

plt.show()



