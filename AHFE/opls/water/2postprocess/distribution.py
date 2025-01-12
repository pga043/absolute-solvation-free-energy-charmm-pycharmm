import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm

"""
Probability distribution: Gaussian distribution
P_0(dU) = -1/(sqrt(2*pi*sigma^2)) * exp{ -[dU - <dU>_0]^2/2*sigma^2}

sigma^2 = <dU^2>_0 - (<dU>_0)^2

"""

dU = np.genfromtxt(sys.argv[1], dtype=None, skip_header=0, delimiter=' ')

average = np.average(dU)

sigma = np.sqrt(np.average(dU**2) - average**2 )

# plotting -ve of probability distribution (check the original formula above
probability = 1/np.sqrt((2*np.pi*sigma**2)) * np.exp(-(dU - average)**2/(2*sigma**2))

#print(probability)
plt.hist(dU, bins=50, density=True)
plt.plot(dU, probability, '.', color='red' )
plt.xlabel('dU')
plt.ylabel('Probability Density')
plt.suptitle(str(sys.argv[1]))
plt.show()

#---------------- confirm using seaborn ----------------
#sns.distplot(dU, probability, hist=False, kde=True)
#sns.histplot(probability, kde=True)
#plt.show()



