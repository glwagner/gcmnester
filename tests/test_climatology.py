import sys; sys.path.append('../gcmnester')
import climatology
import matplotlib.pyplot as plt
import numpy as np

# Path to ECCO data
eccopath = '/data5/glwagner/Numerics/patches/data/ECCOv4r2'

# Fields to extract
fields = ['THETA', 'SALT', 'NVELMASS', 'EVELMASS']

# Initialize the ecco climatology object and describe fields
ecco = climatology.ECCO(eccopath=eccopath)
ecco.describe_fields(fields=fields)

# Initialize data
ecco.init_data(fields, grid='LatLon')

# Extract global climatological data for the month of January
lat, lon, z, yday, data = ecco.extract_globe(
    fields, [0])

# Extract regional North Pacific climatological data for the month of January
region = [5, 55, 100, -100]
paclat, paclon, z, pacyday, pacdata = ecco.extract_region(
    region, fields, [0, 6])

# Close dataset
ecco.close_data()

# Make a plot!
zlev = 0
titles = [['Temperature', 'Salinity'], ['E velocity', 'N velocity']]

# Global plot
T = np.ma.masked_array(data['THETA'],    np.isnan(data['THETA']))
S = np.ma.masked_array(data['SALT'],     np.isnan(data['SALT']))
U = np.ma.masked_array(data['EVELMASS'], np.isnan(data['EVELMASS']))
V = np.ma.masked_array(data['NVELMASS'], np.isnan(data['NVELMASS']))

fig, axs = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True)
axs[0, 0].pcolormesh(lon, lat, T[zlev, :, :])
axs[0, 1].pcolormesh(lon, lat, S[zlev, :, :])
axs[1, 0].pcolormesh(lon, lat, U[zlev, :, :])
axs[1, 1].pcolormesh(lon, lat, V[zlev, :, :])

fig.suptitle("Climatological global fields for one-month average around "
    "yearday {}".format(yday))
axs[0, 0].set_ylabel('Latitude')
axs[1, 0].set_ylabel('Latitude')
axs[1, 0].set_xlabel('Longitude')
axs[1, 1].set_xlabel('Longitude')

for i in range(2): 
    for j in range(2):
        axs[i, j].set_title(titles[i][j])


# Pacific plot
pactitles = ["yearday {}".format(pacyday[0]), "yearday {}".format(pacyday[1])]
pacT = np.ma.masked_array(pacdata['THETA'],    np.isnan(pacdata['THETA']))
pacS = np.ma.masked_array(pacdata['SALT'],     np.isnan(pacdata['SALT']))
pacU = np.ma.masked_array(pacdata['EVELMASS'], np.isnan(pacdata['EVELMASS']))
pacV = np.ma.masked_array(pacdata['NVELMASS'], np.isnan(pacdata['NVELMASS']))
pacSp = pacU**2.0 + pacV**2.0

fig2, axs2 = plt.subplots(ncols=2, sharex=True, sharey=True)
axs2[0].pcolormesh(paclon, paclat, pacSp[0, zlev, :, :])
axs2[1].pcolormesh(paclon, paclat, pacSp[1, zlev, :, :])

fig2.suptitle("Climatological monthly-averaged velocity "
    "around yeardays {} and {}".format(pacyday[0], pacyday[1]))
axs2[0].set_ylabel('Latitude')
axs2[0].set_xlabel('Longitude')
axs2[1].set_xlabel('Longitude')

for i in range(2): 
    axs2[i].set_title(pactitles[i])

plt.show()
