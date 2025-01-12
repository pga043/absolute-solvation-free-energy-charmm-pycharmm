import sys
import numpy as np
import matplotlib.pyplot as plt

nwins = 25
anadir = 'post'

for i in range(1, nwins):
    j = i + 1
    print(f'win{i} -> win{j}')
    # Load energies
    ## Forward direction
    f_win1 = np.loadtxt(f'win{i}/{anadir}/win{j}.dat',dtype=float, skiprows=1, usecols=-1)
    f_win0 = np.loadtxt(f'win{i}/{anadir}/win{i}.dat',dtype=float, skiprows=1, usecols=-1)

    ## Backward direction
    b_win1 = np.loadtxt(f'win{j}/{anadir}/win{j}.dat',dtype=float, skiprows=1, usecols=-1)
    b_win0 = np.loadtxt(f'win{j}/{anadir}/win{i}.dat',dtype=float, skiprows=1, usecols=-1)

    # difference always, final - initial
    frwd = np.subtract(f_win1, f_win0)
    bcwd = np.subtract(b_win0, b_win1)

    plt.hist(frwd, bins=50, density=True, label='forward')
    plt.hist(-bcwd, bins=50, density=True, label='backward')
    plt.title(f'Energy Difference Distributions')
    plt.xlabel(r'$\Delta$U (kcal/mol)',fontsize=12)
    plt.ylabel('Count',fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.show()

quit()
