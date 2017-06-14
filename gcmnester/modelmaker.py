from __future__ import division

import os
import numpy as np
import pyotps, globotopo

from gcmnester import climatology
from gcmnester import gcmgrids
from gcmnester import openboundaries



class nestedmodel:
    def __init__(self, domain, dx, zF, dates,
        name       = 'test',        # Just the name.
        hgridtype  = 'LatLon',      
        dy         = None,
        oblist     = None,
    ):
        """A nested model class.

        Args:
            domain (list): Latitutde/Longitude limits of the nested model
                in the form domain = [south, north, west, east]. The specified
                domain will be adjusted if necessary.

            dx: longitudinal resolution of the model in degrees. If dy is not
                provided, it is set to dx.

            zF (array): An array-like object that defines the faces of the 
                vertical grid. The first element is the bottom of the vertical
                domain and the last element is the surface (so normally, 
                zgrid[-1]=0). The convention used is positive upwards so that
                the bottom at zgrid[0] should be negative. If zgrid is not given, 
                zgrid defaults to the same zgrid as the parent model, if it can
                be found.

            dates: Start and end date of the model in the form 
                dates = [start, end]
        
            name (str): The name of the nested child model.
        
            hgridtype (str): A string identifying the type of grid on which the 
                nested model will be integrated by MITgcm. The only option right
                now is 'LatLon'.

            dy (float): An optional latitudinal resolution. If left unset it 
                defaults to dx.

            #oblist: List of open boundaries. For example, if the southern boundary
            #    is open then oblist = ['south']. If all four boundaries are open
            #    then oblist = ['south', 'north', 'east', 'west']
        """

        if dy is None: dy = dx
        
        self.domain  = domain
        self.dx      = dx
        self.dy      = dy
        self.dates   = dates
        self.oblist  = oblist
        self.name    = name

        self.hgrid   = self.init_hgrid(hgridtype)
        self.zgrid   = gcmgrids.zgrid(zF)

        if oblist is not None:
            self.obcs = {}
            for ob in self.oblist: 
                self.obcs[ob] = openboundaries.obc(ob, self.hgrid, self.zgrid)


    def init_hgrid(self, hgridtype):
        """Initialize the nested model grid. This function takes 
        a single string as argument which specifies the type of the grid
        to initialize.
        """

        if hgridtype is 'LatLon':
            return gcmgrids.latlongrid(self.domain, self.dx, dy=self.dy)
        else:
            raise NotImplementedError("The only type of grid we "
                "can handle now is a LatLon grid.")


    def load_bathymetry(self, datapath=None):
        """Load bathymetry for the model domain using globotopo."""

        if datapath is not None:
            topodata = globotopo.smithsandwell(datapath=datapath)
        else:
            topodata = globotopo.smithsandwell()

        if (not topodata.minlat < self.hgrid.south < topodata.maxlat or
            not topodata.minlat < self.hgrid.north < topodata.maxlat):
            raise ValueError("The grid extends beyond reach of the "
                "Smith-Sandwell bathymetric map.")

        # Interpolate bathymetry
        self.bathy = topodata.interp_to_grid(self.hgrid.xC, self.hgrid.yC,
            method='linear')
    
                
    def load_tidal_bcs(self, otpspath=None, otpstype=None, 
            constits=['m2', 's2']):
        """Load tidal boundary conditions for the specified open boundaries
        using pyotps.

        Args:
            otpspath (str): Path to the OTPS data. See the pyotps docs for 
                more information.

            otpstype (str): Type of the OTPS data; either 'v1' or 'nc'. See
                the pyotps docs for more information.

            constits (list): A list of strings identifying the tidal constituents
                to load. If specified as 'all', all available tidal constituents
                are loaded.
        """

        # Initialize the tidal driver
        if otpspath is not None and otpstype is not None:
            self.tidaldriver = pyotps.tidaldriver(
                otpspath=otpspath, otpstype=otpstype)
        elif otpspath is not None:
            self.tidaldriver = pyotps.tidaldriver(otpspath=otpspath)
        elif otpstype is not None:
            self.tidaldriver = pyotps.tidaldriver(otpstype=otpstype)
        else:
            self.tidaldriver = pyotps.tidaldriver()
        

        # Load tidal data for each boundary condition
        for ob in self.obcs.keys():
            self.obcs[ob].get_tidal_data(self.tidaldriver, constits=constits)







class simpletides(nestedmodel):
    def __init__(self, domain, dx, zF,
        name        = 'test',        
        hgridtype   = 'LatLon',      
        dates       = None,
        dy          = None,
        oblist      = None,
    ):
        """Initialize a 'simple tides' nested model.

        Args:
            domain (list): Latitutde/Longitude limits of the nested model
                in the form domain = [south, north, west, east]. The specified
                domain will be adjusted if necessary.

            dx: longitudinal resolution of the model in degrees. If dy is not
                provided, it is set to dx.

            zF (array): An array-like object that defines the faces of the 
                vertical grid. The first element is the bottom of the vertical
                domain and the last element is the surface (so normally, 
                zgrid[-1]=0). The convention used is positive upwards so that
                the bottom at zgrid[0] should be negative. If zgrid is not given, 
                zgrid defaults to the same zgrid as the parent model, if it can
                be found.

            name (str): The name of the nested child model.
        
            hgridtype (str): A string identifying the type of grid on which the 
                nested model will be integrated by MITgcm. The only option right
                now is 'LatLon'.

            dates: Start and end date of the model in the form 
                dates = [start, end]. Defaults to [2005, 2015]
            
            dy (float): An optional latitudinal resolution. If left unset it 
                defaults to dx.

            oblist: List of open boundaries. For example, if the southern boundary
                is open then oblist = ['south']. If all four boundaries are open
                then oblist = ['south', 'north', 'east', 'west']

            otpspath (str): A string specifying the path that contains both the 
                Oregon Tidal Prediction System compiled executables and TPXO data.

            otpstype (str): The string 'v1' or 'nc' specifying whether the OTPS is
                version 'v1' or version 'nc'. The two versions are incompatible.

        """

        # Defaults
        if dates is None:
            dates = [2005, 2015]
        if dy is None:
            dy = dx
        if oblist is None:
            oblist = ['south', 'north', 'west', 'east']

        # Initialize the nested model.
        nestedmodel.__init__(self, domain, dx, zF, dates, name=name, 
            hgridtype=hgridtype, dy=dy, oblist=oblist) 





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

    # Cap latitudinal limits
    if bathydom[0] < -90: bathydom[0] = -90
    if bathydom[1] < 90:  bathydom[1] = 90

    lat, lon, bathy = globaltopo.get_region(bathydom)

    return lat, lon, bathy



def get_ecco_z():
    """Return the zF values corresponding to the ECCO vertical grid."""

    zF =  [    0.0,    -10.0,    -20.0,   -30.0,     -40.0,
             -50.0,    -60.0,    -70.0,   -80.0,     -90.0,
            -100.0,   -110.6,   -121.5,  -133.4,    -147.1,
            -163.4,   -183.6,   -208.7,   -240.1,   -278.7,
            -325.3,   -380.3,   -443.7,   -515.1,   -593.7,
            -678.6,   -768.4,   -862.1,   -958.4,  -1056.5,
           -1155.7,  -1255.9,  -1357.7,  -1463.1,  -1575.6,
           -1699.7,  -1839.6,  -1999.1,  -2180.1,  -2383.7,
           -2610.2,  -2859.8,  -3132.3,  -3427.8,  -3746.3,
           -4087.8,  -4452.3,  -4839.8  ,-5250.3,  -5683.8, 
           -6117.3
    ]

    return np.array(zF)



