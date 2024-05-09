# -*- coding: utf-8 -*-
"""
Doing some basic Ray Tracing through the Ionosphere

Given a single ray, trace its path through the Ionosphere until it hits the ground
We choose to ignore strictly reflected waves, considering only the refracted wave

We're going to start with the following implementation of basic ionosphere raytracing:
    https://github.com/MIST-Experiment/iricore/commit/f6c726582cd7c5877f0cca02bdbfa7fd8377e55a
"""

import numpy as np

import iri2016.profile as iri  # Ionosphere Sim
from datetime import datetime

import matplotlib.pyplot as plt
import cartopy.crs as ccrs     # Maps
import cartopy

############################# USER INPUTS #####################################
# RF Wave - Viewed As Free and Globally Tracked
SIGNAL = {
            'freq': 10e6,         # 10MHz
            'mag': 1,             # 1dB, by definition
            'phase': 0,           # 0 Degrees of phase, by definition
            
            'pos': np.array([[],[],[]]), # Lat, Lon, Alt - All Histories
            'dir': [0, 0, 0],            # dLat, dLon, dAlt - normalized direction vector relative to EARTH
            };
DIRECTION_THETA = np.deg2rad(90); # Degrees off of boresight. 0 points to Sky, 90 points across land.
DIRECTION_PHI = np.deg2rad(90-45);   # Polar degrees. 90-0 points along equator, 90-90 points up.
TX_DIRECTION = np.array([np.sin(DIRECTION_THETA)*np.cos(DIRECTION_PHI),
                         np.sin(DIRECTION_THETA)*np.sin(DIRECTION_PHI),
                         np.cos(DIRECTION_THETA)]);
# RF Transmitter
TRANSMITTER = {
            'pos': np.array([[40.7128],[-74.0060],[0.0]]), # Lat, Lon, Alt
            'dir': TX_DIRECTION, # Relative to the Globe. deltaLat, dLon, dAlt
            'txWave': SIGNAL,
            };

ALTITUDES = np.linspace(0, 1000, 1); # Layers of the Ionosphere allowed to be simulated
DATETIME = datetime(2020, 6, 1, 0, 0, 0); # An input to the IRI2016 model

############################### CONSTANTS #####################################
EARTH_R = 6.3781e6; # (m)

############################### FUNCTIONS #####################################
def raytrace(TXprop, signal, altitudes, dt):
    """
        For a given Transmitter with properties `TXprop`, transmitting `signal` at 
        a time `dt`, we trace the signal until it either reaches space or reaches 
        the ground
    """
    
    # Pull Transmitting Direction, then change local direction to global direction
    sigDir = np.add(signal['dir'], TXprop['dir']); # Add preexisting signal direction to transmitting direction
    sigPos = np.append(signal['pos'], TXprop['pos'], axis=1);
    sigPos = np.tile(sigPos, 30);
    
    # TODO: Transfer local direction to global direction.

    for i in range(1, sigPos.shape[1]):
        sigPos[:, i] = np.add(sigPos[:, i-1], 1*sigDir);
    
    signal['pos'] = sigPos;
    signal['dir'] = sigDir;
    return signal; #With a frequency, magnitude, and phase


############################### DRIVER ########################################
# Kick off the Raytracing Alg:
signalPath = raytrace(TRANSMITTER, TRANSMITTER['txWave'], ALTITUDES, DATETIME);

# Plot the Path:
plt.figure();
plt.title("Ray Path for Signal with Direction Along Theta: " + str(DIRECTION_THETA*180/np.pi)
          + " and Longitude: " + str(DIRECTION_PHI*180/np.pi));
ax = plt.axes(projection=ccrs.PlateCarree());
ax.coastlines();
ax.add_feature(cartopy.feature.BORDERS);

lats = signalPath['pos'][0];
lons = signalPath['pos'][1];

ax.scatter(lons, lats, label='Equirectangular Straight Line');
ax.set_global();