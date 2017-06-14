from __future__ import division

import numpy as np
from scipy.interpolate import griddata

# Schematic:
#
#           north
#            ____
#    west   |    |  east
#           |____|
#
#           south
#
# Note the following MITgcm Truths According to JMC:
#
#   1. On a southern boundary, tracers and U (the tangential velocity)
#       are specified at the first wet point, while V (normal vel)
#       is specified at the second wet point.
#
#   2. On a northern boundary, tracers, U, and V are all specified at
#       the first wet point.
#
#   3. On a western boundary, tracers and V (the tangential velocity)
#       are specified at the first wet point, while U (normal vel)
#       is specified at the second wet point.
#
#   4. On an eastern boundary, tracers, U, and V are all specified at
#       the first wet point.


class obc:
    def __init__(self, edge, hgrid, zgrid):
        """A class that describes open boundaries in MITgcm, 
        which is on an Arakawa C-grid.

        Args:
            edge (str): A string specifying on which edge of the hgrid to put
                the open boundary. The options are 'south', 'north', 'west', and
                'east'.

            hgrid: A gcmgrid hgrid object with attributes xC, yC, 
                xG, and yG.
            
            zgrid: A gcmgrid zgrid object with attributes zC and zF.
        """

        if edge is 'south':
            self.xC = hgrid.xC[0, :] 
            self.yC = hgrid.yC[0, :] 

            self.xU = hgrid.xU[0, :]
            self.yU = hgrid.yU[0, :]

            # Second wet point
            self.xV = hgrid.xV[1, :]
            self.yV = hgrid.yV[1, :]

        elif edge is 'north':
            self.xC = hgrid.xC[-1, :] 
            self.yC = hgrid.yC[-1, :] 

            self.xU = hgrid.xU[-1, :]
            self.yU = hgrid.yU[-1, :]

            self.xV = hgrid.xV[-1, :]
            self.yV = hgrid.yV[-1, :]

        elif edge is 'west':
            self.xC = hgrid.xC[:, 0]
            self.yC = hgrid.yC[:, 0]
        
            # Second wet point
            self.xU = hgrid.xU[:, 1]
            self.yU = hgrid.yU[:, 1]

            self.xV = hgrid.xV[:, 0]
            self.yV = hgrid.yV[:, 0]

        elif edge is 'east':
            self.xC = hgrid.xC[:, -1]
            self.yC = hgrid.yC[:, -1]
        
            self.xU = hgrid.xU[:, -1]
            self.yU = hgrid.yU[:, -1]

            self.xV = hgrid.xV[:, -1]
            self.yV = hgrid.yV[:, -1]
         

        self.zC = zgrid.zC
        self.zF = zgrid.zF

        self.nb = self.xC.size
        self.nz = self.zC.size


 
    def get_tidal_data(self, tidaldriver, constits=['m2', 's2']):
        """Load tidal amplitudes and phases into obc object."""

        (self.tidalconstits, latout, lonout, self.amps, self.phases
            ) = tidaldriver.extract_amp_phase(
            self.yC, self.xC, constits=constits, var='z', ocegeo='oce')



    def interp_tracer_field(self, fieldname, sdata, sxC, syC, szC, 
        method='linear', fill_value=np.nan):
        
        """Interpolate data onto the open boundary.
    
        Args:
            fieldname (str): The name of the field to be interpolated.
                typically 'T' for temperature or 'S' for salinity.

            sdata (array): The three-dimensional source data to interpolate from.
    
            sxC (array): The tracer longitudes on which the source
                data lives. If sxC is one dimensional, the data is assumed
                to be on a uniform grid, which permits much faster interpolation.

            syC (array): Like sxC but for latitude.
            
            szC (array): The vertical cell centers where the source data is 
                defined. szC must be a 1d array.

            method (str): The interpolation method to use.
        """
    
        snx, sny, snz = sdata.shape

        # Construct 3D array of points
        X = np.zeros((snx, sny, snz), dtype=np.float32)
        Y = np.zeros((snx, sny, snz), dtype=np.float32)
        Z = np.zeros((snx, sny, snz), dtype=np.float32)

        # Fill
        for k in range(snz):
            X[:, :, k] = sxC
            Y[:, :, k] = syC
            Z[:, :, k] = szC(k)

        for k in range(self.nz):
            Xi[:, k] = self.xC[:]
            Yi[:, k] = self.yC[:]
            Zi[:, k] = self.zgrid.zC[k]

        # Use scipy's griddata function. Could be costly...
        setattr(self, fieldname, 
            griddata((X, Y, Z), data, 
            (Xi.ravel(), Yi.ravel(), Zi.ravel()), 
            method=method, fill_value=fill_value))

