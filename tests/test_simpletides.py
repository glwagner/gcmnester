import sys; sys.path.append('..')

import matplotlib.pyplot as plt
import numpy as np
import time

from gcmnester import modelmaker
from gcmnester import climatology

# Some paths
otpspath = '/data5/glwagner/Numerics/patches/pyotps/OTPS2'
topopath = '/data5/glwagner/Numerics/patches/globotopo/data/topo_18.1.img'
eccopath = '/data5/glwagner/Numerics/patches/data/ECCOv4r2'


# Define the domain
south = -25
north = -15
west  = -27
east  = 3

domain = [south, north, west, east]
dx = 1.0/2.0 # resolution in degrees
z  = modelmaker.get_ecco_z()

# Instantiate the model
tidemodel = modelmaker.SimpleTideModel(domain, dx, z, name='test')

# Tasks:    
#   1. Load bathymetry
#   2. Load tidal bcs
#   3. Load parent model and set open boundary conditions
#       and initial condition to domain-averaged stratification.
# ----------------------------------------------------------------------------- 
tidemodel.load_bathymetry(datapath=topopath)
tidemodel.load_tidal_bcs(otpspath=otpspath, constits=['m2'])

T, S = tidemodel.interp_ecco_stratification(eccopath=eccopath, month=0)
tidemodel.set_stratification(T, S)


# Plot
fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(16, 8)) 
    #sharex=True, sharey=True)

#axes[0].imshow(tidemodel.bathy)
axes[1].plot(T, tidemodel.zgrid.zC)
axes[2].plot(S, tidemodel.zgrid.zC)

plt.show()
