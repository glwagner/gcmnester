import numpy as np
import pyotps
import globotopo




# ----------------------------------------------------------------------------- 

class nestedmodel(object):
    def __init__(self, domain, res, dates, bcs,
        name       = 'test',        # Just the name.
        parent     = 'ECCO',        # Name of the parent model
    ):
    """Initialize the nested model.

        Args:
            domain (list):  Latitutde/Longitude limits of the nested model
                            in the form domain = [south, north, west, east].

            res:            Resolution of the model in degrees.

            dates:          Start and end date of the model 
                            in the form dates = [start, end]

            bcs:            Boundary conditions for the model in the form
                            bcs = [bc_south, bc_north, bc_west, bc_east].
                            Each bc is a string indicating 'open' or 'land'.

    """

        self.domain = domain
        self.res = res
        self.dates = dates
        self.parent = parent
        self.name = name

    def build_grid(self):
        """Build the horizontal grid of the child model."""

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
        


    
