import os
import numpy as np
import matplotlib.pyplot as plt

def look_at_input(inputpath, icNames=None, bathyName='bathy.bin'):

    # Default initial condition names
    if icNames is None:
        ics = {
            'T': 'initial_T.bin', 
            'S': 'initial_S.bin', 
            'U': 'initial_U.bin', 
            'V': 'initial_V.bin' 
        }

