import os
import numpy as np
import netCDF4 as nc


class ECCO:
    def __init__(self, 
        eccopath='../../data/ECCOv4r2'
    ):
        """Initialize an ecco object giving access to the
            ECCO data."""


        if not os.path.isdir(eccopath):
            raise ValueError("Check your eccopath. "
                "The directory %s does not exist" % eccopath)

        self.eccopath = eccopath
        self.fieldDescriptions = get_ecco_field_descriptions()


    def init_data(self, fields, grid='LatLon'):
        """Initialize a bunch of NetCDF objects corresponding to climatological 
            ECCO fields.

            Args:
                fields (list): A list of strings corresponding to data fields to 
                    initialize in self.data.

                grid (str): The grid on which the data is stored. The options
                    are 'LatLon' for ECCO data interpolated onto a 0.5 deg
                    Latitude-Longitude grid, or 'LLC' for ECCO data on its
                    native nominal 1 deg LLC0090 grid.


            Note:
                ECCO data is stored in NetCDF3 files. Each variable has its
                own file. Using the python netCDF4 package, the data is loaded
                by referencing the grid with 

                    fileObj['field']['field'][time, dep, lat, lon], 
                
                where time, dep, lat, and lon are indices to be extracted. The grid
                has shape (12, 50, 360, 720), corresponding to a data from 12 years 
                on a 0.5 deg with 50 vertical levels. The convention is that depth  
                is *positive*, which means it gets larger and more positive towards
                negative-z. Each file contains the fields dep, lat, lon, tim, 
                FIELD (the actual data), and index vectors i1, i2, i3, i4.
        """

        self.fieldDescriptions = get_ecco_field_descriptions()
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




    def extract_globe(self, fields, months, grid=None):
        """Extract global, full-depth data for the indicated fields and months
        on the indicated grid.


        Args:
            fields (list): A list of strings corresponding to data fields to 
                extract.

            months (list): A list of integers corresponding to months to
                extract. The count starts are 0 and each value must be between
                0 and 11.

            grid (str): The grid on which the data is stored. The options
                are 'LatLon' for ECCO data interpolated onto a 0.5 deg
                Latitude-Longitude grid, or 'LLC' for ECCO data on its
                native nominal 1 deg LLC0090 grid.


        Returns:
            A tuple of numpy arrays containing the extracted latitude, 
            longitude, z-coordinate, yday, and a dictionary whos keys
            are the extracted fields and whos values are arrays of the data.


        """
        # Check the specified months
        if type(months) is not list:
            raise ValueError("Invalid months parameter. Parameter must be a\n"
                "list even if it only has one element.")
        elif not set(months).issubset(range(11)):
            raise ValueError("Invalid months parameter. Each value in then"
                "list of months to extract must lie between 0 (January 15)\n"
                "and 11 (December 15).")
            
        # If initialized data exists, match grid of initialized data
        if hasattr(self, 'datagrid') and self.datagrid is not None:
            if grid is None or grid is self.datagrid:
                grid = self.datagrid
            else:
                raise ValueError("Invalid grid parameter. There is existing\n"
                    "initialized data on another grid. The grid parameter must be\n"
                    "identical or None, or the existing data must be cleared.")

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

        yday = self.data['tim'][months]
        if len(yday) == 1: yday = yday[0]

        if dataInitializedHere:
            self.close_data()

        return lat, lon, z, yday, data

    

    def extract_regional_stratification(self, region, months, grid=None):
        """Extract profiles of temperature and salinity for a specified 
        region and time-frame.

        Args:
            region: An array-like input of the form 
                region = [south, north, east, west] that defines the 
                latitude-longitude limits of the region to be extracted.
                The input latitudes must lie between -90 and 90 and the
                input longitudes must lie between -180 and 180.
                For example, to extract a region between 20S and 40N, and
                30W and 5E, set box = [-20, 40, -30, 5].

            months (list): A list of integers corresponding to months to
                extract. The count starts are 0 and each value must be between
                0 and 11.

            grid (str): The grid on which the data is stored. The options
                are 'LatLon' for ECCO data interpolated onto a 0.5 deg
                Latitude-Longitude grid, or 'LLC' for ECCO data on its
                native nominal 1 deg LLC0090 grid.

        Returns:
            Two 1 or 2d arrays T, S with dimension T[time, dep]. If only one 
            month is specified, the arrays are one-dimensional.
        """
        lat, lon, z, yday, data = self.extract_region(
            region, ['THETA', 'SALT'], months, grid=grid)

        T = np.nanmean(data['THETA'], axis=(-1, -2))
        S = np.nanmean(data['SALT'],  axis=(-1, -2))

        return T, S



    def extract_region(self, region, fields, months, grid=None):
        """Extract a rectangular region from the LatLon data.

        Args:
            region: An array-like input of the form 
                region = [south, north, east, west] that defines the 
                latitude-longitude limits of the region to be extracted.
                The input latitudes must lie between -90 and 90 and the
                input longitudes must lie between -180 and 180.
                For example, to extract a region between 20S and 40N, and
                30W and 5E, set box = [-20, 40, -30, 5].

            fields (list): A list of strings corresponding to data fields to 
                extract.

            months (list): A list of integers corresponding to months to
                extract. The count starts are 0 and each value must be between
                0 and 11.

            grid (str): The grid on which the data is stored. The options
                are 'LatLon' for ECCO data interpolated onto a 0.5 deg
                Latitude-Longitude grid, or 'LLC' for ECCO data on its
                native nominal 1 deg LLC0090 grid.


           Returns:
            A tuple of numpy arrays containing the extracted latitude, 
            longitude, z-coordinate, yday, and a dictionary whos keys
            are the extracted fields and whos values are arrays of the data.
            
            NOTE: If the region spans the International Date Line (lon=+/-180),
            then values east of the International Date Line will be assigned 
            values greater than +180 to preserve continuity of the grid.
            """

        if not hasattr(region, 'index') or len(region) < 4:
            raise ValueError("The region parameter must be an indexed "
                "collection (a list, numpy array, tuple, etc.) of length "
                "four (2 latitude and 2 longitude bounds each).")
        else:
            try: 
                iswest = region[0] < 0
            except:
                raise ValueError("The region parameter must be an indexed "
                    "collection of latitudes and longitudes.")

        # Check the specified months
        if not set(months).issubset(range(11)):
            raise ValueError("Invalid months parameter. Each value in the\n"
                "list of months to extract must lie between 0 (January 15)\n"
                "and 11 (December 15).")
            
        # If initialized data exists, match grid of initialized data
        if hasattr(self, 'datagrid') and self.datagrid is not None:
            if grid is None or grid is self.datagrid:
                grid = self.datagrid
            else:
                raise ValueError("Invalid grid parameter. There is existing\n"
                    "initialized data on another grid. The grid parameter must be\n"
                    "identical or None, or the existing data must be cleared.")

        # If no input is given, set to default LatLon grid.
        elif grid is None:
            grid = 'LatLon'

        # Initialize data if it hasn't been already
        if not hasattr(self, 'data') or self.data is None:
            self.init_data(fields=fields, grid=grid)
            dataInitializedHere = True
        else:
            dataInitializedHere = False

        # Extract sides of the box, flipping south/north coordinates if need be
        south, north = np.sort(np.array(region)[[0, 1]])
        west, east = np.array(region)[[2, 3]]

        # Raise hell if something is amiss
        if not -180 <= east <= 180 or not -180 <= west <= 180:
            raise ValueError("Longitudes must lie between +/-180 degrees")
        elif south < -90 or north > 90:
            raise ValueError("Latitudes must lie between +/-90 degrees.")
        elif east == west or south == north:
            raise ValueError("Latitudes and longitudes must not be unique!")

        # Wrapping is needed if coordinates cross the prime meridian
        if west > east: 
            acrossdateline = True
        else:
            acrossdateline = False


        # Find indices for cutting, taking care not to produce bad indices.
        if grid is 'LatLon':

            # Square grid
            jsouth = searchsorted_left (self.data['lat'][:, 0], south)
            jnorth = searchsorted_right(self.data['lat'][:, 0], north)

            iwest = searchsorted_left (self.data['lon'][0, :], west)
            ieast = searchsorted_right(self.data['lon'][0, :], east)

        elif grid is 'LLC':
            raise NotImplementedError("Extracting regional domains from an "
                "LLC grid is not implemented.")

        # Cut out depth and ydays.
        z    = -self.data['dep'][:]
        yday = self.data['tim'][months]
        if len(yday) == 1: yday = yday[0]

        # Cut out latitude, longitude, and data, taking into account 
        # whether or not the indicated region crosses the dateline.
        if not acrossdateline: 
            # Easy case
            rlon = self.data['lon'][jsouth:jnorth, iwest:ieast]
            rlat = self.data['lat'][jsouth:jnorth, iwest:ieast]

            # Extract data from NetCDF files by referencing indices.
            data = {}
            for fld in fields:
                data[fld] = self.data[fld][months, :, 
                    jsouth:jnorth, iwest:ieast].squeeze()

        elif acrossdateline:

            # Wrapping case
            nrlat = jnorth - jsouth
            (nreast, nrwest) = (ieast, self.data['nlon']-iwest)

            rlon  = np.zeros((nrlat, nreast+nrwest), dtype=np.float64)
            rlat  = np.zeros((nrlat, nreast+nrwest), dtype=np.float64)

            # Shift longitude coordinate to preserve monotonicity of data
            rlon[:, nrwest:] = self.data['lon'][jsouth:jnorth, :ieast] + 360
            rlon[:, :nrwest] = self.data['lon'][jsouth:jnorth, iwest:]

            rlat[:, nrwest:] = self.data['lat'][jsouth:jnorth, :ieast]
            rlat[:, :nrwest] = self.data['lat'][jsouth:jnorth, iwest:]

            data = {}
            for fld in fields:
                data[fld] = np.zeros(
                    (len(months), self.data['ndep'], nrlat, nreast+nrwest),
                     dtype=np.float64)

                data[fld][:, :, :, nrwest:] = self.data[fld][months, :, 
                    jsouth:jnorth, :ieast].squeeze()

                data[fld][:, :, :, :nrwest] = self.data[fld][months, :, 
                    jsouth:jnorth, iwest:].squeeze()


        if dataInitializedHere:
            self.close_data()

        return rlat, rlon, z, yday, data


    def close_data(self):
        """Close all data files and set associated values to None."""
        for fld in self.datafiles.keys():
            self.datafiles[fld].close()

        self.datafiles = None
        self.datagrid = None
        self.data = None


    def describe_fields(self, fields=None):
        """Print a description of the stored ECCO fields."""

        self.fieldDescriptions = get_ecco_field_descriptions()

        if fields is None:
            fields = fieldDescriptions.keys()

        for fld in fields:
            print("{:>12s} : {}".format(fld, self.fieldDescriptions[fld]))















def searchsorted_left(data, leftside):
    """Return the index of data so that data[ileft] is either the first 
    index in data or lying just left of leftside."""

    if len(data.shape) > 1:
        raise(ValueError, "Input data must be one-dimensional.")

    ileft = np.max([
        np.searchsorted(data, leftside,  side='left')-1, 0])

    return ileft




def searchsorted_right(data, rightside):
    """Return the index of data so that data[iright] is either the last
    index in data or lying just right of rightside"""

    if len(data.shape) > 1:
        raise(ValueError, "Input data must be one-dimensional.")

    iright = np.min([
        np.searchsorted(data, rightside,  side='right')+1, np.size(data)-1])

    return iright




def get_ecco_field_descriptions():
    """Return a dictionary of descriptions of the fields given in the README 
        in ECCO's interp_climatology directory."""

    fieldDescriptions = {
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


    return fieldDescriptions


