import os
import numpy as np
import netCDF4 as nc


class ecco:
    def __init__(self, 
        eccopath='../../data/ECCOv4r2'
    ):
        """Initialize an ecco object giving access to the
            ecco climatological data."""


        if not os.path.isdir(eccopath):
            raise ValueError("Check your eccopath. "
                "The directory %s does not exist" % eccopath)

        self.eccopath = eccopath
        self.get_field_descriptions()


    def init_data(self, fields=['THETA', 'SALT'], grid='LatLon'):
        """Initialize a bunch of NetCDF objects corresponding to climatological 
            ECCO fields.

            Args:
                grid (str): The grid on which the data is stored. The options
                    are 'LatLon' for ECCO data interpolated onto a 0.5 deg
                    Latitude-Longitude grid, or 'LLC' for ECCO data on its
                    native nominal 1 deg LLC0090 grid.

                fields (list): A list of strings corresponding to data fields to 
                    initialize in self.data.

            Note:
                ECCO data is stored in NetCDF3 files. Each variable has its
                own file. Using the python netCDF4 package, the data is loaded
                by referencing the grid with 

                    fileObj['field']['field'][time, dep, lat, lon], 
                
                where time, dep, lat, and lon are indices to be extracted. The grid
                has shape (12, 50, 360, 720), corresponding to a data from 12 years 
                on a 0.5 deg with 50 vertical levels. The convention is that depth is 
                *positive*, which means it gets larger and more positive towards
                negative-z. Each file contains the fields dep, lat, lon, tim, 
                FIELD (the actual data), and index vectors i1, i2, i3, i4.
        """

        self.get_field_descriptions()
        validFields = self.fieldDescriptions.keys()

        # Check inputs and assign path to data
        if type(fields) is not list:
            raise ValueError("The parameter 'fields' must be a list of "
                "strings corresponding to ECCO climatological fields")
        elif not set(fields).issubset(validFields):
            raise ValueError("One or more of the elements of 'fields' is not "
                "a valid ECCO climatological field.")
        elif grid is 'LatLon':
            climatopath = self.eccopath + "/interp_climatology"
        elif grid is 'LLC':
            climatopath = self.eccopath + "/nctiles_climatology"
            raise ValueError("The LLC grid is not yet supported.")
        else:
            raise ValueError("Specified grid must be 'LatLon' or LLC "
                "(Latitude-Longitude-Cap).")

        if not os.path.isdir(climatopath):
            raise ValueError("Your eccopath exists, but it does "
                "not contain climatology data in %s." % climatopath)

        if hasattr(self, 'datagrid') and self.datagrid is not None:
            if self.datagrid is not grid:
                raise ValueError("Cannot initialize new data on "
                    "the {} grid alongside prexisting data ".format(grid) +
                    "on the {} grid.\n".format(self.datagrid) +
                    "Either delete the existing data with 'close_data()' "
                    "or select grid = {}".format(self.datagrid))
        else:
            self.datagrid = grid

        # Initialize dictionaries if they do not exist.
        if not hasattr(self, 'data'):
            self.data = {}
            self.datafiles = {}

        # Load the data!
        for fld in fields:

            fieldpath = "%s/%s.0001.nc" % (climatopath, fld)

            self.datafiles[fld] = nc.Dataset(fieldpath, 'r',
                format='NETCDF4')

            self.data[fld] = self.datafiles[fld][fld]

        # Load Lat, Lon, Dep, Tim from the first field
        self.data['lat'] = self.datafiles[fields[0]]['lat']
        self.data['lon'] = self.datafiles[fields[0]]['lon']
        self.data['dep'] = self.datafiles[fields[0]]['dep']
        self.data['tim'] = self.datafiles[fields[0]]['tim']

        self.data['ntim'] = len(self.datafiles[fields[0]]['i1'])
        self.data['ndep'] = len(self.datafiles[fields[0]]['i2'])
        self.data['nlat'] = len(self.datafiles[fields[0]]['i3'])
        self.data['nlon'] = len(self.datafiles[fields[0]]['i4'])




    def extract_global_data(self, fields=['THETA', 'SALT'], 
        months=[0], grid=None):

        # If initialized data exists, match grid of initialized data
        if hasattr(self, 'datagrid') and self.datagrid is not None:
            if grid is None or grid is self.datagrid:
                grid = self.datagrid
            else:
                raise ValueError("Initialized data exists on a "
                    "different grid than the one specified")

        # If no input is given, set to default LatLon grid.
        elif grid is None:
            grid = 'LatLon'

        # Initialize data if it hasn't been already
        if not hasattr(self, 'data') or self.data is None:
            self.init_data(fields=fields, grid=grid)
            dataInitializedHere = True
        else:
            dataInitializedHere = False

        # Extract data from NetCDF files by referencing indices.
        data = {}
        for fld in fields:
            data[fld] = self.data[fld][months, :, :, :].squeeze()

        # Use convention that positive z-values are upwards.
        z   = -self.data['dep'][:]
        lat =  self.data['lat'][:, :]
        lon =  self.data['lon'][:, :]

        if len(months) == 1:
            yday = self.data['tim'][months[0]]
        else:
            yday = self.data['tim'][0]

        if dataInitializedHere is True:
            self.close_data()

        return lat, lon, z, yday, data




    def close_data(self):
        """Close all data files and set associated values to None."""
        for fld in datafiles.keys():
            self.datafiles.close()
        self.datafiles = None
        self.datagrid = None
        self.data = None


    def get_field_descriptions(self):
        """Load descriptions of the fields given in the README in
            ECCO's interp_climatology directory."""

        self.fieldDescriptions = {
            'ADVeHEFF'  :  "Eastward Advective Flux of eff ice thickn",
            'ADVe_SLT'  :  "Eastward Advective Flux of Salinity",
            'ADVeSNOW'  :  "Eastward Advective Flux of eff snow thickn",
            'ADVe_TH'   :  "Eastward Advective Flux of Pot.Temperature",
            'ADVnHEFF'  :  "Northward Advective Flux of eff ice thickn",
            'ADVn_SLT'  :  "Northward Advective Flux of Salinity",
            'ADVnSNOW'  :  "Northward Advective Flux of eff snow thickn",
            'ADVn_TH'   :  "Northward Advective Flux of Pot.Temperature",
            'ADVr_SLT'  :  "Vertical   Advective Flux of Salinity",
            'ADVr_TH'   :  "Vertical   Advective Flux of Pot.Temperature",
            'ADVxHEFF'  :  "U Comp. Advective Flux of eff ice thickn",
            'ADVx_SLT'  :  "U Comp. Advective Flux of Salinity",
            'ADVxSNOW'  :  "U Comp. Advective Flux of eff snow thickn",
            'ADVx_TH'   :  "U Comp. Advective Flux of Pot.Temperature",
            'ADVyHEFF'  :  "V Comp. Advective Flux of eff ice thickn",
            'ADVy_SLT'  :  "V Comp. Advective Flux of Salinity",
            'ADVySNOW'  :  "V Comp. Advective Flux of eff snow thickn",
            'ADVy_TH'   :  "V Comp. Advective Flux of Pot.Temperature",
            'DFeEHEFF'  :  "Eastward Diffusive Flux of eff ice thickn",
            'DFeE_SLT'  :  "Eastward Diffusive Flux of Salinity",
            'DFeESNOW'  :  "Eastward Diffusive Flux of eff snow thickn",
            'DFeE_TH'   :  "Eastward Diffusive Flux of Pot.Temperature",
            'DFnEHEFF'  :  "Northward Diffusive Flux of eff ice thickn",
            'DFnE_SLT'  :  "Northward Diffusive Flux of Salinity",
            'DFnESNOW'  :  "Northward Diffusive Flux of eff snow thickn",
            'DFnE_TH'   :  "Northward Diffusive Flux of Pot.Temperature",
            'DFrE_SLT'  :  "Vertical Diffusive Flux of Salinity    (Explicit part)",
            'DFrE_TH'   :  "Vertical Diffusive Flux of Pot.Temperature (Explicit part)",
            'DFrI_SLT'  :  "Vertical Diffusive Flux of Salinity    (Implicit part)",
            'DFrI_TH'   :  "Vertical Diffusive Flux of Pot.Temperature (Implicit part)",
            'DFxEHEFF'  :  "U Comp. Diffusive Flux of eff ice thickn",
            'DFxE_SLT'  :  "U Comp. Diffusive Flux of Salinity",
            'DFxESNOW'  :  "U Comp. Diffusive Flux of eff snow thickn",
            'DFxE_TH'   :  "U Comp. Diffusive Flux of Pot.Temperature",
            'DFyEHEFF'  :  "V Comp. Diffusive Flux of eff ice thickn",
            'DFyE_SLT'  :  "V Comp. Diffusive Flux of Salinity",
            'DFyESNOW'  :  "V Comp. Diffusive Flux of eff snow thickn",
            'DFyE_TH'   :  "V Comp. Diffusive Flux of Pot.Temperature",
            'DRHODR'    :  "Stratification: d.Sigma/dr (kg/m3/r_unit)",
            'ETAN'      :  "Free Surface Height Anomaly (Ocean-Ice Interface)",
            'EVELMASS'  :  "Eastward Mass-Weighted Comp of Velocity (m/s)",
            'EVELSTAR'  :  "Eastward Comp of Bolus Velocity (m/s)",
            'GM_PsiX'   :  "GM Bolus transport stream-function' : U component",
            'GM_PsiY'   :  "GM Bolus transport stream-function' : V component",
            'MXLDEPTH'  :  "Mixed-Layer Depth (>0)",
            'NVELMASS'  :  "Northward Mass-Weighted Comp of Velocity (m/s)",
            'NVELSTAR'  :  "Northward Comp of Bolus Velocity (m/s)",
            'oceFWflx'  :  "net surface Fresh-Water flux into the ocean (+=down), >0 decreases salinity",
            'oceQnet'   :  "net surface heat flux into the ocean (+=down), >0 increases theta",
            'oceQsw'    :  "net Short-Wave radiation (+=down), >0 increases theta",
            'oceSPflx'  :  "net surface Salt flux rejected into the ocean during freezing, (+=down),",
            'oceSPtnd'  :  "salt tendency due to salt plume flux >0 increases salinity",
            'oceTAUE'   :  "Eastward surface wind stress, >0 increases eVel",
            'oceTAUN'   :  "Northward surface wind stress, >0 increases nVel",
            'oceTAUX'   :  "U Comp. surface wind stress, >0 increases uVel",
            'oceTAUY'   :  "V Comp. surface wind stress, >0 increases vVel",
            'PHIBOT'    :  "Bottom Pressure Pot.(p/rho) Anomaly",
            'PHIHYD'    :  "Hydrostatic Pressure Pot.(p/rho) Anomaly",
            'RHOAnoma'  :  "Density Anomaly (=Rho-rhoConst)",
            'SALT'      :  "Salinity",
            'SFLUX'     :  "total salt flux (match salt-content variations), >0 increases salt",
            'SIarea'    :  "SEAICE fractional ice-covered area [0 to 1]",
            'SIatmFW'   :  "Net freshwater flux from atmosphere & land (+=down)",
            'SIatmQnt'  :  "Net atmospheric heat flux, >0 decreases theta",
            'sIceLoad'  :  "sea-ice loading (in Mass of ice+snow / area unit)",
            'SIheff'    :  "SEAICE effective ice thickness",
            'SIhsnow'   :  "SEAICE effective snow thickness",
            'TFLUX'     :  "total heat flux (match heat-content variations), >0 increases theta, W/m2",
            'THETA'     :  "Potential Temperature",
            'UVELMASS'  :  "U Mass-Weighted Comp of Velocity (m/s)",
            'UVELSTAR'  :  "U Comp of Bolus Velocity (m/s)",
            'VVELMASS'  :  "V Mass-Weighted Comp of Velocity (m/s)",
            'VVELSTAR'  :  "V Comp of Bolus Velocity (m/s)",
            'WVELMASS'  :  "Vertical Mass-Weighted Comp of Velocity",
            'WVELSTAR'  :  "Vertical Comp of Bolus Velocity (m/s)",
        }


    def describe_fields(self, fields=['THETA']):
        """Print a description of the stored ECCO fields."""

        self.get_field_descriptions()

        if fields is 'all':
            fields = fieldDescriptions.keys()

        for fld in fields:
            print("{:>12s} : {}".format(fld, self.fieldDescriptions[fld]))




# ----------------------------------------------------------------------------- 


class worldoceanatlas:
    def __init__(self, atlaspath=''):
        """A class for accessing and viewing data from the 
            world ocean atlas."""

