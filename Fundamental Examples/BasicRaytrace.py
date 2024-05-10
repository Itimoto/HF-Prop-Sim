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

########################## ADDN'L SERVICE FUNCS ###############################
def globalDirCosines(theta, phi, deg=1):
    """
        Return the direction cosines for a given Theta, Phi. Default: Degrees
        
        Adding the Direction Cosines to [Latitude, Longitude, Altitude] will
        yield a change in Theta off of boresight, Phi along the polar plane.
        For Theta: 0 Points to Sky, 90 points across land.
        For Phi: 0 Points along equator, 90 points due North
    """
    if deg:
        theta = np.deg2rad(theta);
        phi = np.deg2rad(phi);
        
    # We subtract from 90 degrees to convert Phi to Longitude.
    return np.array([np.sin(theta)*np.cos((np.pi/2)-phi),
                     np.sin(theta)*np.sin((np.pi/2)-phi),
                     np.cos(theta)]);

############################# USER INPUTS #####################################
# RF Wave - Viewed As Free and Globally Tracked
SIGNAL = {
            'freq': 10e6,         # 10MHz
            'mag': 1,             # 1dB, by definition
            'phase': 0,           # 0 Degrees of phase, by definition
            
            'pos': np.array([[],[],[]]), # Lat, Lon, Alt - All Histories
            'dir': [0, 0, 0],            # dLat, dLon, dAlt - normalized direction vector relative to EARTH
            };

DIRECTION_THETA = 45; # Degrees off of boresight. 0 points to Sky, 90 points across land.
DIRECTION_PHI = 45;   # Polar degrees. 0 points along equator, 90 points due North.
TX_DIRECTION = globalDirCosines(DIRECTION_THETA, DIRECTION_PHI);

# RF Transmitter
TRANSMITTER = {
            'pos': np.array([[40.7128],[-74.0060],[0.0]]), # Lat, Lon, Alt
            'dir': TX_DIRECTION, # Relative to the Globe. deltaLat, dLon, dAlt
            'txWave': SIGNAL,
            };

ALTITUDES = np.linspace(60, 1000, 20); # Layers of the Ionosphere allowed to be simulated
DATETIME = datetime(2020, 6, 1, 0, 0, 0); # An input to the IRI2016 model

############################### CONSTANTS #####################################
EARTH_R = 6.3781e3;     # (km)
MAX_IRI_ALT = 2000;     # (km) Max Altitude for IRI2016.
MIN_IRI_ALT = 10;       # (km) Min Altitude for IRI2016. 
GND_HEIGHT = 0;         # (km) Sea Level. Collision with Ground considered here.

############################### FUNCTIONS #####################################
def collisionPoint(signal, targetAlt):
    """
    Return Collision Point

    Parameters
    ----------
    signal : SIGNAL
        Dict with 'pos' and 'dir' arrays, each of which denote movement along
        the geodetic coordinate system: [Latitude, Longitude, Altitude]
        Assume 'pos' in km. 'dir' must be normalized.
    targetAlt : FLOAT
        Altitude (km) of layer of interest.

    Returns
    -------
    Tuple (New Pos[Lat, Lon, Alt], Distance Traveled)
    'pos' at Intersection Point, as well as the Distance Traveled
        by the ray.
    If no collision, Distance Traveled = -1
    """
    
    sigPos = signal['pos']; sigDir = signal['dir'];
    collisionAlt = targetAlt - sigPos[2, -1]; # Look at most recent signal data
    # + Difference: target above. - Difference: target below. 
    
    if collisionAlt*sigDir[2] > 0: # + Diff + Dir = Collision with Time
        # Find the intersection point:    
        newPos = np.add(sigPos[:, -1], (collisionAlt/sigDir[2])*sigDir); # Intersection Point Occurs here
        #disttravel = collisionAlt / sigPos[-1, 2]; # (Rise)/(Rise/Run)
        
        # Use the Great-Circle Distance Formula to find Distance Traveled
        # https://www.geeksforgeeks.org/great-circle-distance-formula/
        a, b, x, y = np.deg2rad((sigPos[0, -1], newPos[0], 
                                 sigPos[1, -1], newPos[1]));
        distTraveled = EARTH_R * np.arccos(np.cos(a)*np.cos(b)*np.cos(x-y) + np.sin(a)*np.sin(b));
        return (newPos, distTraveled); # Collision occurs at newPos after travelling distTraveled
    else:
        return (signal['pos'], -1); # No collision occurs.
    
    

def raytrace(TXprop, signal, altitudes, dt):
    """
        For a given Transmitter with properties `TXprop`, transmitting `signal` at 
        a time `dt`, we trace the signal until it either reaches space or reaches 
        the ground
    """
    
    # Pull Transmitting Position & Direction, apply to Initial Signal Info
    sigDir = np.add(signal['dir'], TXprop['dir']);
    sigPos = np.append(signal['pos'], TXprop['pos'], axis=1);
    
    signal['pos'] = sigPos; signal['dir'] = sigDir;
    
    # Simulate all given Ionospheric Layers, inidicated by `altitudes`
        # We can do this later in more complicated examples.
    print(collisionPoint(signal, 2000));
    # Begin Transmission
    # while ray hasn't hit the ground or space yet:
        # Determine where next Intersection occurs
        # Find the closest Ionospheric layer
        # Find the intersection position
            # Find path distance to intersection (for transmission path loss)
        # Pull the Refraction Angle at the intersection point
        # Apply the outgoing angle at the interface using Snell's law
    
    sigPos = np.tile(sigPos, 30);

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
plt.title("Sample Path for Signal with Direction Along Theta: " + str(DIRECTION_THETA)
          + " and Longitude: " + str(DIRECTION_PHI));
ax = plt.axes(projection=ccrs.PlateCarree());
ax.coastlines();
ax.add_feature(cartopy.feature.BORDERS);

lats = signalPath['pos'][0];
lons = signalPath['pos'][1];

ax.scatter(lons, lats, label='Equirectangular Straight Line');
ax.set_global();