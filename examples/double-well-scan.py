#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from fourier_DVR_1D import Domain_Fourier_DVR_1D


# settings
m0 = 1.0
m1 = 5.0
n_m = 100
x_min = -5.0
x_max = 5.0
n_DVR = 400
n_g = 1000
n_plot = 50
D0 = 50.0
a = 1.5
scale = 1.5
V = lambda x: (D0 / a**4) * (x**2 - a**2)**2
E_max = 2 * D0


m = np.linspace(m0, m1, n_m)
x = np.linspace(x_min, x_max, n_g)
V_x = V(x)

# prepare domain
domain = Domain_Fourier_DVR_1D(x_min, x_max, n_DVR)

# calculate eigenenergies and eigenstates for all masses
print 'calculate'
E = []
for i in range(len(m)):

    # solve
    e = domain.solve(m[i], V, n_states=n_plot, calc_eigenstates=False)
    E.append(e)

    # print progress
    print '%6.2f %%' % (100.0 * (i+1) / n_m)
E = np.array(E)
print

print 'plot'
for i in range(n_m):

    plt.figure(figsize=(10, 6))

    plt.subplot(121)
    plt.plot([m[i], m[i]], [0, E_max], 'k--', color='gray')
    for j in range(n_plot):
        plt.plot(m, E[:,j], '-')
    plt.xlim(m0, m1)
    plt.ylim(0, E_max)
    plt.xlabel('mass')
    plt.ylabel('energy')

    plt.subplot(122)
    plt.xlabel('x')
    plt.yticks([])
    plt.twinx()
    plt.ylabel('energy')
    for j in range(n_plot):
        E_x = E[i, j] * np.ones_like(x)
        E_x[E_x < V_x] = None
        plt.plot(x, E_x)
    plt.plot(x, V_x, 'k-', lw=2)
    plt.xlim(-3, 3)
    plt.ylim(0, E_max)

    plt.tight_layout()

    plt.savefig('%03d.png' % i, dpi=150)
    plt.close()
    print '%6.2f %%' % (100.0 * (i+1) / n_m)