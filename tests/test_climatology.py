import sys; sys.path.append('../gcmnester')
import climatology
import matplotlib.pyplot as plt
import numpy as np

# Path to ECCO data
eccopath = '/data5/glwagner/Numerics/patches/data/ECCOv4r2'

# Fields to extract
fields = ['THETA', 'SALT', 'NVELMASS', 'EVELMASS']

# Initialize the ecco climatology object and describe fields
ecco = climatology.ecco(eccopath=eccopath)
ecco.describe_fields(fields=fields)

# Initialize and extract data
ecco.init_data(fields=fields, grid='LatLon')
lat, lon, z, yday, data = ecco.extract_global_data(fields, months=[0], grid='LatLon')

# Make a plot!
zlev = 0
titles = [['Temperature', 'Salinity'], ['E velocity', 'N velocity']]

# Generate masked arrays.
T = np.ma.masked_array(data['THETA'],    np.isnan(data['THETA']))
S = np.ma.masked_array(data['SALT'],     np.isnan(data['SALT']))
U = np.ma.masked_array(data['EVELMASS'], np.isnan(data['EVELMASS']))
V = np.ma.masked_array(data['NVELMASS'], np.isnan(data['NVELMASS']))

fig, axs = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True)
axs[0, 0].pcolormesh(lon, lat, T[zlev, :, :])
axs[0, 1].pcolormesh(lon, lat, S[zlev, :, :])
axs[1, 0].pcolormesh(lon, lat, U[zlev, :, :])
axs[1, 1].pcolormesh(lon, lat, V[zlev, :, :])

fig.suptitle('Global fields')
axs[i, j].set_xlabel('Longitude')
axs[i, j].set_ylabel('Latitude')

for i in range(2): 
    for j in range(2):
        axs[i, j].set_title(titles[i][j])



plt.show()
