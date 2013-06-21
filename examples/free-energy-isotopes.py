#!/usr/bin/env python

import itertools
import numpy as np
import matplotlib.pyplot as plt

from fourier_DVR_1D import Domain_Fourier_DVR_1D


# settings
T0 = 0.1
T1 = 20
n_T = 30
m0 = 1.0
m1 = 2.0
n_l = 51
x_min = -5.0
x_max = 5.0
n_DVR = 400

# # potential
# k = 1.0
# V = lambda x: 0.5 * k * x * x
# V_id = 'harmonic, k=%f' % k

# # potential
# D0 = 10.0
# a = 1.0
# V = lambda x: (D0 / a**4) * (x**2 - a**2)**2
# V_id = 'double well, D0=%f, a=%f' % (D0, a)

# potential
k = 1.0
V = lambda x: k * x**4
V_id = 'quartic, k=%f' % k


# mass switching function
def m_l(lam, m0, m1):
    return m0 * m1 / (lam * np.sqrt(m0) + (1-lam) * np.sqrt(m1))**2


print V_id
print

# masses
lam = np.linspace(0.0, 1.0, n_l)
dl = 1.0 / (n_l - 1.0)
m = m_l(lam, m0, m1)

# colors to cycle through
colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

# prepare domain
domain = Domain_Fourier_DVR_1D(x_min, x_max, n_DVR)

# solve for each m
print 'solving for %i lambda values' % n_l
E = []
for i in range(n_l):
    EE, psi_four = domain.solve(m[i], V)
    E.append(EE)
    print '%.3f' % lam[i]
E = np.array(E)
print

# match the ground state energy of the harmonic oscillator
k_harm = 4 * E[0, 0]**2 * m[0]
print 'harmonic k =', k_harm
print

# some thermodynamics
dA_dlam = []
Ts = np.linspace(T0, T1, n_T)
for T in Ts:

    Q = np.exp(-E / T).sum(axis=1)
    A = - T * np.log(Q)
    dA_dlam.append(np.gradient(A, dl))

    print '               T =', T
    print '    full Delta A =', A[-1] - A[0]
    print 'midpoint Delta A =', dA_dlam[-1][n_l/2]
    print

# plot double well and HO eigenenergies as function of mass
plt.figure()
for i in range(n_DVR):
    c = colors.next()
    E_harm = (0.5 + i) * np.sqrt(k_harm / m)
    plt.plot(lam, E_harm, '--', color=c)
    plt.plot(lam, E[:, i], '-', color=c)
plt.xlim(0, 1)
plt.ylim(0, 20)
plt.xlabel('lambda')
plt.ylabel('energy')
plt.title('eigenenergies | %s, dashed - HO (k=%f)' % (V_id, k_harm))

plt.figure()
for i in range(n_T):
    plt.plot(lam, dA_dlam[i], 'k-')
plt.xlim(0, 1)
plt.ylim(ymax=0)
plt.xlabel('lambda')
plt.ylabel('d A / d lambda')
plt.title(V_id)

plt.show()
