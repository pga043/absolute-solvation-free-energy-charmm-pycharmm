import sys
import numpy as np
import matplotlib.pyplot as plt

#m = np.genfromtxt('matrix.dat', dtype=None, skip_header=0, delimiter=' ')

#win1 = m[0, 51:102] - m[0,0:51]
#win2 = m[1, 51:102] - m[1, 0:51]
#print(win1)
#print(win2)
#quit()
#print(m1.shape)
#quit()

win1 = np.genfromtxt(sys.argv[1], dtype=None, skip_header=0, delimiter=' ')

win2 = np.genfromtxt(sys.argv[2], dtype=None, skip_header=0, delimiter=' ')

plt.hist(win1, bins=50, density=True)
plt.hist(-win2, bins=50, density=True)

plt.show()
