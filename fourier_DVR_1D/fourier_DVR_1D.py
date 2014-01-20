import logging

import numpy as np
import scipy.linalg


class Domain_Fourier_DVR_1D(object):
    """
    Solves the Schroedinger equation on a finite one-dimensional interval using
    the discrete variable representation method with the Fourier sine basis:

    phi_k = sqrt(2 / (b-a)) * sin(k * pi * (x-a) / (b-a)), k=1..n_DVR

    For details, see for example Appendix 2 of:
    J. Chem. Phys. 104, 8433 (1996)
    http://dx.doi.org/10.1063/1.471593
    """

    def __init__(self, a, b, n_DVR):
        """Constructs the domain with given end points and basis size.

        :param a: lower bound of domain
        :param b: upper bound of domain
        :param n_DVR: number of basis functions

        """

        # store domain parameters
        self._a = a
        self._b = b
        self._n_DVR = n_DVR

        # no previous calculation was performed
        self._m = None
        self._V = None

        # update spectral decomposition of the position operator
        self._update_X_in_Fourier_basis()

    def _update_X_in_Fourier_basis(self):
        """Decompose the position operator in the Fourier sine basis.

        Eigenvalues and eigenvectors of the position operator
        in the Fourier sine basis are stored in `self`.

        """

        logging.debug('X decomposition | start')

        a = self._a
        b = self._b
        n_DVR = self._n_DVR

        # construct position operator in Fourier sine basis
        i = np.arange(1, n_DVR+1, dtype=float)[:, np.newaxis].repeat(n_DVR, axis=1)
        j = i.transpose()
        ipj = i + j
        ipjmod1 = (ipj) % 2
        div = np.where(ipjmod1,  ((i-j) * ipj)**2, 1)
        fact = -8.0 * (b-a) / np.pi**2
        X_four = fact * np.where(ipjmod1, (i*j) / div, 0.0)

        x, phi_four = scipy.linalg.eigh(X_four)
        x = x + (a+b) / 2.0

        self._x = x
        self._phi_four = phi_four

        logging.debug('X decomposition | done')

    def _update_T_four(self, m):
        """Build the kinetic energy operator and store it in `self`.

        :param m: mass

        """

        if self._m is None or (m != self._m):

            self._m = m

            l = self._b - self._a
            t = (0.5 / m) * (np.pi / l)**2 * np.arange(1, self._n_DVR + 1)**2
            self._T_four = np.diagflat(t)

            logging.debug('build T | done')

    def _update_V_four(self, V):
        """Build the potential energy operator and store it in `self`.

        :param V: potential energy - real-space vectorized function

        """

        if self._V is None or (V != self._V):

            self._V = V

            phi_four = self._phi_four
            V_x = np.diagflat(V(self._x))
            self._V_four = np.dot(np.dot(phi_four, V_x), phi_four.transpose())

            logging.debug('build V | done')

    def solve(self, m, V, n_states=None, calc_eigenstates=True):
        """Solve the Schroedinger equation on this domain.

        :param m: mass
        :param V: potential energy - real-space vectorized function
        :param n_states: number of states to calculate
        "param calc_eigenstates: whether to return eigenstates as well

        Returns eigenenergies and (optionally) eigenstates of the Hamiltonian
        sorted by eigenenergy magnitude.

        """

        logging.debug('solve | start')

        if n_states is None:
            eigvals = None
        else:
            eigvals = (0, n_states)

        # update kinetic energy and potential operators, if needed
        self._update_T_four(m)
        self._update_V_four(V)

        # construct the Hamiltonian
        H_four = self._T_four + self._V_four
        logging.debug('solve | Hamiltonian built')

        # solve
        eigvals_only = not calc_eigenstates
        result = scipy.linalg.eigh(H_four,
                                   eigvals=eigvals,
                                   eigvals_only=eigvals_only,
                                   overwrite_a=True)
        logging.debug('solve | done')

        return result


    def grid(self, x, psi_four):
        """Evaluate states on a real-space grid.

        :param x: real-space grid - 1D array
        :param psi_four: states in Fourier basis

        Returns states on grid x.

        """

        logging.debug('grid | start')

        n_DVR, n_out = psi_four.shape

        if n_DVR != self._n_DVR:
            data = (self._n_DVR, n_DVR)
            err = 'Wrong dimension of states. Expected %d, got %d.' % data
            raise ValueError(err)

        # convenience
        a = self._a
        b = self._b

        # flag points outside the [a, b] interval
        outside = np.logical_and((x > a), (x < b))

        # construct Fourier sine basis functions on real-space grid
        norm = np.sqrt(2 / (b-a))
        four_1_grid = np.exp(1.0j * np.pi * (x-a) / (b-a))
        four_1_grid *= outside
        k = np.arange(1, n_DVR + 1, dtype=complex)
        four_grid = norm * np.imag(four_1_grid[:, np.newaxis]**k)

        # project states to real-space grid
        psi_grid = np.dot(four_grid, psi_four).transpose()

        logging.debug('grid | done')

        return psi_grid
