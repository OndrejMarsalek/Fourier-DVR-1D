#!/usr/bin/env python

import sys
import logging
import numpy as np
import matplotlib.pyplot as plt

from fourier_DVR_1D import Domain_Fourier_DVR_1D


# settings
m = 1.0
x_min = -5.0
x_max = 5.0
n_DVR = 500
n_g = 1001
D0 = 5.0
a = 1.0
V = lambda x: (D0 / a**4) * (x**2 - a**2)**2
n_states = 15
n_plot = 15
scale = 2.0

# enable debug info printing
FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG, format=FORMAT)

# solve
domain = Domain_Fourier_DVR_1D(x_min, x_max, n_DVR)
E, E_four = domain.solve(m, V, n_states=n_states)

# evaluate eigenstates on grid
x = np.linspace(x_min, x_max, n_g)
psi_x = domain.grid(x, E_four[:,:n_plot])

# print energies
print 'Eigenenergies'
print '-------------'
print
for i, e in enumerate(E[:n_plot]):
    print '%3i %12.6f' % (i, e)

# plot eigenstates
plt.figure()
plt.subplots_adjust(left=0.05, right=0.95,
                    bottom=0.05, top=0.95)
plt.plot(x, V(x), 'k-', lw=2)
for i in range(n_plot):
    plt.plot([x[0], x[-1]], [E[i], E[i]], '--', color='gray')
for i in range(n_plot):
    plt.plot(x, scale * psi_x[i] + E[i])
plt.xlim(x_min, x_max)
plt.ylim(-E[0], E[n_plot])
plt.show()
