import sys; sys.path.append('..')

import matplotlib.pyplot as plt
import numpy as np
import time

from gcmnester import modelmaker

# Some paths
otpspath = '/data5/glwagner/Numerics/patches/pyotps/OTPS2'
topopath = '/data5/glwagner/Numerics/patches/globotopo/data/topo_18.1.img'


# Define the domain
south = -25
north = -15
west = -27
east = 3

domain = [south, north, west, east]
dx = 1.0/12.0 # resolution in degrees

zF = modelmaker.get_ecco_z()

# Make the model
tidemodel = modelmaker.simpletides(domain, dx, zF, 
    name = 'test', 
)

tidemodel.load_bathymetry(datapath=topopath)
tidemodel.load_tidal_bcs(otpspath=otpspath, otpstype='v1', 
    constits=['m2', 's2'])


# Plot
fig, axes = plt.subplots(ncols=1, nrows=2, figsize=(16, 8), 
    sharex=True, sharey=True)

axes[0].imshow(tidemodel.hgrid.xC)
axes[1].imshow(tidemodel.bathy)

plt.show()
