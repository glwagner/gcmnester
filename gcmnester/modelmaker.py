from __future__ import division

import os
import numpy as np
import pyotps, globotopo

from scipy.interpolate import interp1d

from gcmnester import climatology
from gcmnester import gcmgrids
from gcmnester import openboundaries
from gcmnester import initialcondition



class NestedModel:
    def __init__(self, domain, dx, zF, dates,
        name       = 'test',        # Just the name.
        hgridtype  = 'LatLon',      
        dy         = None,
        oblist     = None,
        tilesize   = 1,
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

            oblist: List of open boundaries. For example, if the southern boundary
                is open then oblist = ['south']. If unset, oblist defaults to   
                oblist = ['south', 'north', 'east', 'west'].

            tilesize (int): The size of a square tile in integer grid points. The
                default is 60x60, which means the smallest default grid
                that can be specified is nz x 60 x 60.
        """

        if dy is None: dy = dx
        
        self.domain   = domain
        self.dx       = dx
        self.dy       = dy
        self.dates    = dates
        self.oblist   = oblist
        self.name     = name
        self.tilesize = tilesize

        self.hgrid   = self.init_hgrid(hgridtype, tilesize=tilesize)
        self.zgrid   = gcmgrids.VerticalGrid(zF)

        self.init_obcs()
        self.init_ic()
        


    def init_ic(self):
        """Initialize the initial condition."""
        self.ic = initialcondition.InitialCondition(self.hgrid, self.zgrid)



    def init_obcs(self):
        """Initialize the open boundary conditions."""
        self.obcs = {}
        for ob in self.oblist: 
            self.obcs[ob] = openboundaries.OpenBoundaryCondition(
                ob, self.hgrid, self.zgrid)



    def init_hgrid(self, hgridtype, tilesize=1):
        """Initialize the nested model grid. 

        Args:
            hgridtype (str): A string specifying the type of grid to initialize.
                The option are 'LatLonGrid' only at the moment.

            tilesize (int): An optional specification that forces the size of the
                generated grid to be a multiple of 'tilesize'.
        """

        if hgridtype is 'LatLon':
            return gcmgrids.LatLonGrid(self.domain, self.dx, dy=self.dy,
                tilesize=tilesize)
        else:
            raise NotImplementedError("The only type of grid we "
                "can handle now is a LatLon grid.")



    def load_bathymetry(self, datapath=None):
        """Load bathymetry for the model domain using globotopo."""

        if datapath is not None:
            topodata = globotopo.SmithSandwell(datapath=datapath)
        else:
            topodata = globotopo.SmithSandwell()

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
            self.tidaldriver = pyotps.TidalDriver(
                otpspath=otpspath, otpstype=otpstype)
        elif otpspath is not None:
            self.tidaldriver = pyotps.TidalDriver(otpspath=otpspath)
        elif otpstype is not None:
            self.tidaldriver = pyotps.TidalDriver(otpstype=otpstype)
        else:
            self.tidaldriver = pyotps.TidalDriver()
        

        # Load tidal data for each boundary condition
        for ob in self.obcs.keys():
            self.obcs[ob].get_tidal_data(self.tidaldriver, constits=constits)



    def interp_ecco_stratification(self, eccopath=None, month=0, 
        method='linear'):
        """Load the climatological, domain-averaged temperature and
        salinity for the indicated month from ECCO and linearly interpolate
        to the model grid.

        Args:
            eccopath (str): Path to ECCO data. See climatology.ECCO for more
                information.

            month (int): A integer indicating which month among the 
                monthly-averaged climatological fields to extract. For example,
                January=0 and December=11.

            method (str): The type of interpolation to use. The string 
                correpsonds to the 'kind' keyword argument in scipy's
                interp1d.

        Returns:
            A tuple of arrays of temperature and salinity on the model's
            vertical grid.
        """

        if type(month) is not int:
            raise ValueError("The parameter 'month' must be an integer\n"
                "corresponding to the month of the year, where January=0\n"
                "and December=11, for example")
        elif eccopath is not None and type(eccopath) is not str:
            raise ValueError("The specified eccopath must be a string that\n"
                "points to downloaded ECCO data.")

        fields = ['THETA', 'SALT']

        if eccopath is not None:
            ecco = climatology.ECCO(eccopath=eccopath)
        else:
            ecco = climatology.ECCO()

        ecco.init_data(fields)
        eccolat, eccolon, eccoz, yday, data = ecco.extract_region(
            self.domain, fields, [month])

        # Horiontally-averaged temperature and salinity
        # ECCO fields are outputted with format THETA[z, lat, lon]
        eccoT = np.nanmean(data['THETA'], axis=(1, 2))
        eccoS = np.nanmean(data['SALT'],  axis=(1, 2))

        # Extend temperature and depth uniformly to 20,000 meters
        nz = eccoz.size
        eccoT.resize(nz+1)
        eccoS.resize(nz+1)
        eccoz.resize(nz+1)

        eccoz[-1] = -20000
        eccoT[-1] = eccoT[-2]
        eccoS[-1] = eccoS[-2]

        # Construct objects to perform interpolation to vertical grid
        T = interp1d(eccoz, eccoT, kind='linear') 
        S = interp1d(eccoz, eccoS, kind='linear') 
        
        # Interpolate to model vertical grid and return
        return T(self.zgrid.zC), S(self.zgrid.zC)


        def init_gcm_setup(self, setuppath=None):
            """Initialize the directory structure of an MITgcm run.
    
            Args:
                setuppath (str): Path to create for the the MITgcm setup.
            """

            # Default setuppath is the current one
            if setuppath is None:
                setuppath = os.getcwd()

            codepath  = setuppath + "/code"
            buildpath = setuppath + "/build"
            inputpath = setuppath + "/input"
            namepath  = setuppath + "/namelists"

            # Create the necessary directories, checking to see if they exist.
            if not os.path.exists(setuppath):
                os.makedirs(setuppath)

            if not os.path.exists(codepath):
                os.makedirs(codepath)

            if not os.path.exists(buildpath):
                os.makedirs(buildpath)

            if not os.path.exists(inputpath):
                os.makedirs(inputpath)

            if not os.path.exists(namepath):
                os.makedirs(namepath)

            self.setuppath = setuppath
            self.codepath  = codepath
            self.buildpath = buildpath
            self.inputpath = inputpath
            self.namepath  = namepath




class SimpleTideModel(NestedModel):
    def __init__(self, domain, dx, zF,
        name        = 'test',        
        hgridtype   = 'LatLon',      
        dates       = None,
        dy          = None,
        oblist      = None,
        tilesize    = 1,
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

            tilesize (int): The size of a square tile in integer grid points. The
                default is 60x60, which means the smallest default grid
                that can be specified is nz x 60 x 60.
        """

        # Defaults
        if dates is None:
            dates = [2005, 2015]
        if dy is None:
            dy = dx
        if oblist is None:
            oblist = ['south', 'north', 'west', 'east']

        # Initialize the nested model.
        NestedModel.__init__(self, domain, dx, zF, dates, name=name, 
            hgridtype=hgridtype, dy=dy, oblist=oblist, tilesize=tilesize) 





    def set_obcs(self, T, S):
        """Set the constant stratification open boundary conditions for
        the simple tide model. 

        Args:
            T (array): Vertical temperature profile for the stratification.
                Must be 1d and have the same size as z.
            S (array): Vertical salinity profile for the stratification. 
                Must be 1d and have the same size as z.
        """
        for ob in self.oblist: 
            self.obcs[ob].set_constant_stratification(T, S)



    def set_ic(self, T, S):
        """Initialize the initial condition for the simple tide model.

        Args:
            T (array): Vertical temperature profile for the stratification.
                Must be 1d and have the same size as z.
            S (array): Vertical salinity profile for the stratification. 
                Must be 1d and have the same size as z.
        """

        self.ic.set_constant_stratification(T, S)



    def set_stratification(self, T, S):
        """Initialize the open boundaries and initial condition for the simple
        tide model with with a constant stratification.

        Args:
            T (array): Vertical temperature profile for the stratification.
                Must be 1d and have the same size as z.
            S (array): Vertical salinity profile for the stratification. 
                Must be 1d and have the same size as z.
        """

        self.set_obcs(T, S)
        self.set_ic(T, S)










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



