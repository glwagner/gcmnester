import numpy as np
import pyotps
import globotopo




# ----------------------------------------------------------------------------- 

class nestedmodel:
    def __init__(self, domain, dx, dates, bcs,
        name       = 'test',        # Just the name.
        parent     = 'ECCO',        # Name of the parent model
        grid       = 'LatLon',      
        dy         = None
    ):
        """Initialize the nested model.

        Args:
            domain (list): Latitutde/Longitude limits of the nested model
                in the form domain = [south, north, west, east]. The specified
                domain will be adjusted if necessary.

            dx: longitudinal resolution of the model in degrees. If dy is not
                provided, it is set to dx.

            dates: Start and end date of the model 
                in the form dates = [start, end]

            bcs: Boundary conditions for the model in the form
                bcs = [bc_south, bc_north, bc_west, bc_east].
                Each bc is a string indicating 'open' or 'land'.

        """

        if dy is None: dy = dx
        
        self.domain = domain
        self.dx     = dx
        self.dy     = dy
        self.dates  = dates
        self.bcs    = bcs
        self.name   = name
        self.parent = parent
        self.grid   = grid

        if grid is 'LatLon':
            self.init_latlon_grid()
        else:
            raise NotImplementedError("The only type of grid we "
                "can handle now is a LatLon grid.")



    def init_latlon_grid(self):
        """Initialize a latitude-longitude grid with constant spacing."""

        # Pre-calculate some basic properties of the domain
        [south, north, east, west] = self.domain

        Lx = east-west
        Ly = north-south

        # Ensure actual computational domain contains the specification
        self.nx = np.ceil(Lx/self.dx)
        self.ny = np.ceil(Ly/self.dy)

        # Now set actual values. Anchor in the southwest.
        self.west = west
        self.east = west + nx*self.dx

        self.south = south
        self.north = south + ny*self.dy

        self.Lx = self.west-self.east
        self.Ly = self.north-self.south

        # Build the grid: XG, YG, XC, YC 
        
    


    def build_grid(self):
        """Build the horizontal grid of the child model."""


        #self.YG = 

        # Constant Lat/Lon grid. 


    def extract_tide_bcs(self):
        """Extract the tidal amplitudes and phases along the open boundary."""


    def load_bathymetry(self):
        """Load the Smith/Sandwell bathymetry for the nested model region."""

        globalbathy = globotopo.smithsandwell()

        # Pad the domain for interpolation-safety
        (lonpad, latpad) = (1.0, 1.0)

        bathydom = self.domain 
        bathydom[0] -= latpad
        bathydom[1] += latpad
        bathydom[2] -= lonpad
        bathydom[3] += lonpad

        lat, lon, bathy = globalbathy.get_region(bathydom)

        return lat, lon, bathy


    def interpolate_bathymetry(self, bathyLat, bathyLon, bathy):
        """Interpolate the inputted bathymetry onto the child grid."""

        # Linearly interpolate for now.
        #self.bathy = scipy
        


    
