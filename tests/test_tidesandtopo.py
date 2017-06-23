import sys; sys.path.append('..')

import globotopo
import pyotps

import matplotlib.pyplot as plt
import numpy as np
import time

# Subsample topography for speed
subsample = 64

# Define an Atlantic domain
south = -45
north = 55
west = -80
east = 20

domain = [south, north, west, east]
times = {}   

t0 = time.time()

# Extract topography
topodata = globotopo.SmithSandwell(
    datapath='/data5/glwagner/Numerics/patches/globotopo/data/topo_18.1.img')

topolat, topolon, topo = topodata.get_region(domain, subsample=subsample)

times['topo'] = time.time() - t0



t0 = time.time()

# Extract tides
driver = pyotps.TidalDriver(
    otpspath = '/data5/glwagner/Numerics/patches/pyotps/OTPS2',
    otpstype = 'v1')

constits, tidelat, tidelon, amps, phases = driver.extract_amp_phase(
    topolat, topolon)     

times['tides'] = time.time() - t0


# Print timing results
print("(subsample {:d}) {:>24s}: {:7.4f} s".format(
    subsample, "time for topo", times['topo']))
print("(subsample {:d}) {:>24s}: {:7.4f} s".format(
    subsample, "time for tides", times['tides']))


# Plot
fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(24, 10), 
    sharex=True, sharey=True)

im = {}
im['topo']   = axes[0].imshow(topo,         origin='lower', aspect='auto', 
    vmax=0)
im['amps']   = axes[1].imshow(amps['m2'],   origin='lower', aspect='auto', 
    vmax=5)
im['phases'] = axes[2].imshow(phases['m2'], origin='lower', aspect='auto')

i = 0
for k in im.keys():
    plt.colorbar(im[k], ax=axes[i])
    i += 1

plt.show()
