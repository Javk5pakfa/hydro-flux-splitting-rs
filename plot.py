#!/usr/local/bin

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("riemann_exact_solver/solution.dat")

print(np.shape(data))

x = data[:, 0]
rho = data[:, 1]
v   = data[:, 2]
eps = data[:, 3]

plt.subplot(131)
plt.plot(x, rho, "-")

plt.subplot(132)
plt.plot(x, v, "-")

plt.subplot(133)
plt.plot(x, eps, "-")

plt.show()
