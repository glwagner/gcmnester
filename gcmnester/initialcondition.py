from __future__ import division

import numpy as np
from scipy.interpolate import griddata





class InitialCondition:
    def __init__(self, hgrid, zgrid):
        """A class for an MITgcm initial condition."""

        self.nx = hgrid.nx
        self.ny = hgrid.ny
        self.nz = zgrid.nz

        self.xC = hgrid.xC
        self.yC = hgrid.yC

        self.xU = hgrid.xU
        self.yU = hgrid.yU

        self.xV = hgrid.xV
        self.yV = hgrid.yV

        self.hgrid = hgrid
        self.zgrid = zgrid

        self.nx = hgrid.nx
        self.ny = hgrid.ny
        self.nz = zgrid.nz


    def set_constant_stratification(self, T, S):

        # Uniform velocity initial condition
        self.U = np.zeros((self.nz, self.ny, self.nx))
        self.V = np.zeros((self.nz, self.ny, self.nx))

        # Constant temperature and salinity
        self.T = T.reshape(self.nz, 1, 1)*np.ones((self.nz, self.ny, self.nx))
        self.S = S.reshape(self.nz, 1, 1)*np.ones((self.nz, self.ny, self.nx))
