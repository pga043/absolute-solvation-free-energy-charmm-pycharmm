import numpy as np
import matplotlib.pyplot as plt

fwd = np.genfromtxt('forward.txt',dtype=float)
rev = np.genfromtxt('reverse.txt',dtype=float)

plt.plot(np.cumsum(fwd))
plt.savefig('forward_pmf.pdf')
plt.close()

plt.plot(np.cumsum(rev))
plt.savefig('reverse_pmf.pdf')
plt.close()

with open('free_energy.txt','w') as f:
    f.write(f'forward free energy: {sum(fwd)} kcal/mol\n')
    f.write(f'reverse free energy: {sum(rev)} kcal/mol')
