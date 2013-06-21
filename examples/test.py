#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

from fourier_DVR_1D import Domain_Fourier_DVR_1D


# settings
m = 1.0
x_min = -10.0
x_max = 10.0
n_DVR = 200
n_g = 201
V = lambda x: 0.5 * x * x
try:
    n_plot = int(sys.argv[1])
except IndexError:
    n_plot = 4

# solve
domain = Domain_Fourier_DVR_1D(x_min, x_max, n_DVR)
E, E_four = domain.solve(m, V)

# evaluate eigenstates on grid
x = np.linspace(x_min, x_max, n_g)
psi_x = domain.grid(x, E_four[:,:n_plot])

# print energies
for e in E:
    print '%12.6f' % e

# plot eigenstates
plt.figure(figsize=(16, 12))
plt.subplots_adjust(left=0.03, right=0.97,
                    bottom=0.03, top=0.92-0.016*((1+n_plot)/2))

for i in range(n_plot):
    plt.plot(x, psi_x[i], label="E = %f" % E[i])
plt.legend(bbox_to_anchor=(0, 1.0+0.016*((1+n_plot)/2),
                           1.0, 0.016*((n_plot+1)/2)),
           ncol=2, mode="expand", borderaxespad=0.0)

plt.show()
