from __future__ import division

import numpy as np



class HorizontalGrid:
    def __init__(self, xC, yC, xG, yG):
        """A class of horizontal grids that correspond to grids commonly
        used in MITgcm.
    
        Args:
            xC, yC: Two-dimensional arrays that give the latitude
                and longitudes of the center of every horizontal cell.

            xG, yG: Two-dimensional arrays that give the latitudes and longitudes
                of the vertices of every horizontal cell.

        Note: xG, yG have an extra point in either dimension, because there is
            one more vertex than cell center along the two coordinates..
        """

        self.nx, self.ny = xC.shape 

        if xC.shape != yC.shape:
            raise ValueError("The input xC and yC must have the same shape.")
        elif xG.shape != yG.shape:
            raise ValueError("The input xG and yG must have the same shape.")
        elif (self.nx+1, self.ny+1) != xG.shape:
            raise ValueError("The input xG and yG must have an extra grid point"
                "in both dimensions.")

        # Bind basic grid properties.
        self.xC = xC
        self.yC = yC
        self.xG = xG
        self.yG = yG

        # Store bindings for velocity points to obviate confusion.
        self.xV = self.xC
        self.yV = self.yG[:-1]

        self.xU = self.xG[:-1]
        self.yU = self.yC

        self.dxG = self.xG[:, 1:] - self.xG[:, 0:-1]
        self.dyG = self.yG[1:, :] - self.yG[0:-1, :]

        self.dxC = self.xC[:, 1:] - self.xC[:, 0:-1]
        self.dyC = self.yC[1:, :] - self.yC[0:-1, :]





class LatLonGrid(HorizontalGrid):
    def __init__(self, domain, dx, dy=None, tilesize=1):
        """Define properties and a few functions associated with a 
        constant Latitude/Longitude MITgcm grid on the domain
        domain = [south, north, west, east] and with grid spacing (dx, dy).       
        If dy is not given as a keyword parameter it is assumed equal
        to dx.
        
            Args:
                domain (array-like): An array like structure that defines
                    the four boundaries of a latitude-longitude box with
                    the convention domain = [south, north, west, east].

                dx (float): The spacing of the grid. If the keyword parameter
                    'dy' is not defined, this is both the latitudinal and 
                    longitudinal spacing. If dy is defined, dx is just the 
                    longitudinal spacing.
                
                dy (float): The optional latitudinal spacing. If dy is not 
                    specified, it is set equal to dx.

                tilesize (int): The size of a square tile in integer grid points.
                    The default is 1, which produces no change in the size of the
                    grid.
        """

        if dy is None: dy = dx

        self.domain = domain
        self.dx = dx
        self.dy = dy

        [south, north, west, east] = self.domain

        if not -180.0 < west < 180.0 or not -180.0 < east < 180.0:
            raise ValueError("Longitudes must lie between -180 and +180.")
        elif not -90.0 < south < 90.0 or not -90.0 < north < 90.0:
            raise ValueError("Latitudes must lie between -90 and +90.")

        if west > east:
            # Domain crosses dateline
            acrossdateline = True
            east += 360 
        else:
            acrossdateline = False
            
        Lx = east-west
        Ly = north-south

        # Ensure actual computational domain contains the specification
        nx0 = np.ceil(Lx/self.dx)
        ny0 = np.ceil(Ly/self.dy)

        # Round up to nearest multiple of the tilesize
        self.nx = int(tilesize*np.ceil(nx0/tilesize))
        self.ny = int(tilesize*np.ceil(ny0/tilesize))

        # Now set actual values. Anchor in the southwest.
        self.west = west
        self.east = west + self.nx*self.dx

        self.south = south
        self.north = south + self.ny*self.dy

        # Adjust domain so that center is (almost) invariant to the choice of 
        # 'tilesize'
        self.west  = self.west  - self.dx*(self.nx-nx0)//2
        self.east  = self.east  - self.dx*(self.nx-nx0)//2
        self.north = self.north - self.dy*(self.ny-ny0)//2
        self.south = self.south - self.dy*(self.ny-ny0)//2

        self.Lx = self.west-self.east
        self.Ly = self.north-self.south

        # Build the grid: xG, yG, xC, yC 
        xG = np.arange(west,  east,  dx, dtype=np.float64) 
        yG = np.arange(south, north, dy, dtype=np.float64) 

        xC = 0.5*(xG[1:] + xG[0:-1])
        yC = 0.5*(yG[1:] + yG[0:-1])

        # Re-correct coordinates if domain crosses dateline
        if acrossdateline:
            xC[xC > 180] -= 360
            xG[xG > 180] -= 360
        
        # Form and store arrays
        xG, yG = np.meshgrid(xG, yG)
        xC, yC = np.meshgrid(xC, yC)

        HorizontalGrid.__init__(self, xC, yC, xG, yG)
        





class VerticalGrid:
    def __init__(self, zF):
        """Initialize an MITgcm vertical grid.
        
        Arg:
            zF (array): An array-like object that defines the faces of the 
                vertical grid. The first element is the bottom of the vertical
                domain and the last element is the surface (so normally, 
                zgrid[-1]=0). The convention used is positive upwards so that
                the bottom at zgrid[0] should be negative. If zgrid is not given, 
                zgrid defaults to the same zgrid as the parent model, if it can
                be found.
        """

        # Make simple calculations of the vertical grid.
        self.zF  = zF
        self.zC  = 0.5*(zF[:-1] + zF[1:])
        self.dzF = zF[1:] - zF[:-1] 

        # Calculate dz(dzF), or the rate of change of the
        # vertical grid spacing
        self.dzdzF = self.dzF[1:] - self.dzF[:-1] 

        self.nz = self.zC.size

### From "GRID.h":
#
#   | (3) Views showing nomenclature and indexing
#   |     for grid descriptor variables.
#   |
#   |      Fig 3a. shows the orientation, indexing and
#   |      notation for the grid spacing terms used internally
#   |      for the evaluation of gradient and averaging terms.
#   |      These varaibles are set based on the model input
#   |      parameters which define the model grid in terms of
#   |      spacing in X, Y and Z.
#   |
#   |      Fig 3b. shows the orientation, indexing and
#   |      notation for the variables that are used to define
#   |      the model grid. These varaibles are set directly
#   |      from the model input.
#   |
#   | Figure 3a
#   | =========
#   |       |------------------------------------
#   |       |                       |
#   |"PWY"********************************* etc...
#   |       |                       |
#   |       |                       |
#   |       |                       |
#   |       |                       |
#   |       |                       |
#   |       |                       |
#   |       |                       |
#   |
#   |       .                       .
#   |       .                       .
#   |       .                       .
#   |       e                       e
#   |       t                       t
#   |       c                       c
#   |       |-----------v-----------|-----------v----------|-
#   |       |                       |                      |
#   |       |                       |                      |
#   |       |                       |                      |
#   |       |                       |                      |
#   |       |                       |                      |
#   |       u<--dxF(i=1,j=2,k=1)--->u           t          |
#   |       |/|\       /|\          |                      |
#   |       | |         |           |                      |
#   |       | |         |           |                      |
#   |       | |         |           |                      |
#   |       |dyU(i=1,  dyC(i=1,     |                      |
#   | ---  ---|--j=2,---|--j=2,-----------------v----------|-
#   | /|\   | |  k=1)   |  k=1)     |          /|\         |
#   |  |    | |         |           |          dyF(i=2,    |
#   |  |    | |         |           |           |  j=1,    |
#   |dyG(   |\|/       \|/          |           |  k=1)    |
#   |   i=1,u---        t<---dxC(i=2,j=1,k=1)-->t          |
#   |   j=1,|                       |           |          |
#   |   k=1)|                       |           |          |
#   |  |    |                       |           |          |
#   |  |    |                       |           |          |
#   | \|/   |           |<---dxV(i=2,j=1,k=1)--\|/         |
#   |"SB"++>|___________v___________|___________v__________|_
#   |       <--dxG(i=1,j=1,k=1)----->
#   |      /+\                                              .
#   |       +
#   |       +
#   |     "WB"
#   |
#   |   N, y increasing northwards
#   |  /|\ j increasing northwards
#   |   |
#   |   |
#   |   ======>E, x increasing eastwards
#   |             i increasing eastwards
#   |
#   |    i: East-west index
#   |    j: North-south index
#   |    k: up-down index
#   |    u: x-velocity point
#   |    V: y-velocity point
#   |    t: tracer point
#   | "SB": Southern boundary
#   | "WB": Western boundary
#   |"PWX": Periodic wrap around in X.
#   |"PWY": Periodic wrap around in Y.
# -----------------------------------------------------------------------------  





