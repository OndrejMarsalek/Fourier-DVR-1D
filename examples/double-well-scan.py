#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from fourier_DVR_1D import Domain_Fourier_DVR_1D


# settings
m0 = 1.0
m1 = 5.0
n_lam = 100
x_min = -5.0
x_max = 5.0
n_DVR = 400
n_g = 1000
n_plot = 50
D0 = 50.0
a = 1.5
scale = 1.5
V = lambda x: (D0 / a**4) * (x**2 - a**2)**2


def m_l(lam, m0, m1):
    """Mass switching function."""

    return m0 * m1 / (lam * np.sqrt(m0) + (1-lam) * np.sqrt(m1))**2


m = [m_l(lam, m0, m1) for lam in np.linspace(0, 1, n_lam)]
x = np.linspace(x_min, x_max, n_g)
V_x = V(x)

# prepare domain
domain = Domain_Fourier_DVR_1D(x_min, x_max, n_DVR)

# calculate eigenenergies and eigenstates for all masses
E = []
E_four = []
for i in range(len(m)):

    # solve
    e, e_four = domain.solve(m[i], V)
    E.append(e)
    E_four.append(e_four)

    # print progress
    print '%6.2f' % (100.0 * i / n_lam)

E = np.array(E)

for i in range(len(m)):

    plt.figure(figsize=(10, 6))

    plt.subplot(121)
    plt.plot([m[i], m[i]], [0, 2*D0], 'k--', color='gray')
    for j in range(n_plot):
        plt.plot(m, E[:,j], '-')
    plt.xlim(m0, m1)
    plt.ylim(0, 2*D0)
    plt.xlabel('mass')
    plt.ylabel('energy')

    plt.subplot(122)
    plt.xlabel('x')
    plt.yticks([])
    plt.twinx()
    plt.ylabel('energy')
    for j in range(n_plot):
        E_x = E[i, j] * np.ones_like(x)
        E_x /= E_x > V_x
        plt.plot(x, E_x)
    plt.plot(x, V_x, 'k-', lw=2)
    plt.xlim(-3, 3)
    plt.ylim(0, 2*D0)

    plt.tight_layout()

    plt.savefig('%03d.png' % i, dpi=75)
    print i
