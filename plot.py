#!/usr/local/bin

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("solution.dat")

print(np.shape(data))

x = data[:, 0]
rho = data[:, 1]
mom = data[:, 2]

plt.plot(x, rho, ".")
plt.plot(x, rho)
plt.show()
