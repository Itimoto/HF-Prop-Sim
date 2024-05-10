# -*- coding: utf-8 -*-
"""
Doing some basic Ray Tracing through the Ionosphere

Given a single ray, trace its path through the Ionosphere until it hits the ground
We choose to ignore strictly reflected waves, considering only the refracted wave

We're going to start with the following implementation of basic ionosphere raytracing:
    https://github.com/MIST-Experiment/iricore/commit/f6c726582cd7c5877f0cca02bdbfa7fd8377e55a
"""
import pdb # For debugging

import numpy as np
import copy

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
            'freq': 3e6,         # 10MHz
            'mag': 1,             # 1dB, by definition
            'phase': None,        # (Not implemented yet)
            
            'pos': np.array([[],[],[]]), # Lat, Lon, Alt - All Histories
            'dir': [0, 0, 0],            # dLat, dLon, dAlt - normalized direction vector relative to EARTH
            
            'epsilon_r': 1.0006,  # Relative Permittivity in last recorded medium. Approx at Sea Level.
            };

DIRECTION_THETA = 80; # Degrees off of boresight. 0 points to Sky, 90 points across land.
DIRECTION_PHI = -15;   # Polar degrees. 0 points along equator, 90 points due North.
TX_DIRECTION = globalDirCosines(DIRECTION_THETA, DIRECTION_PHI);

# RF Transmitter
TRANSMITTER = {
            'pos': np.array([[40.7128],[-74.0060],[0.01]]), # Lat, Lon, Alt
            'dir': TX_DIRECTION, # Relative to the Globe. deltaLat, dLon, dAlt
            'txWave': SIGNAL,
            };

ALTITUDES = np.linspace(80, 1000, 200); # Layers of the Ionosphere allowed to be simulated
DATETIME = datetime(2020, 6, 1, 0, 0, 0); # An input to the IRI2016 model

############################### CONSTANTS #####################################
EARTH_R = 6.3781e3;     # (km)
MAX_IRI_ALT = 2000;     # (km) Max Altitude for IRI2016.
MIN_IRI_ALT = 80;       # (km) Min Altitude for IRI2016. 
GND_HEIGHT = 0;         # (km) Sea Level. Collision with Ground considered here.

E = 1.60217662e-19;     # (C) Electron Charge
M_E = 9.10938356e-31;   # (kg) Electron Mass
EPSILON_0 = 8.85418782e-12; # (F/m) Vacuum Permittivity
EPSILON_R_ATM = 1.0006; # (F/m) Relative Permittivity of air 
MU_0 = np.pi*4e-7;      # (N/A^2) Vacuum Permeability

############################### FUNCTIONS #####################################
def wrapLatLon(pos, deg=1):
    max_lat = 90.0; min_lat = -90.0;
    max_lon = 180.0; min_lon = -180.0;
    
    # Lat, Lon, Alt
    _pos = np.array(pos);
    if not deg:
        _pos[0], _pos[1] = np.rad2deg([_pos[0], _pos[1]]);
    
    # Wrap latitude within valid range
    _pos[0] = ((_pos[0] - min_lat) % (max_lat - min_lat)) + min_lat
    
    # Wrap longitude within valid range
    _pos[1] = ((_pos[1] - min_lon) % (max_lon - min_lon)) + min_lon

    return _pos;

def polarToCartesian(pos, deg=1):
    # Lat, Lon, Alt
    _pos = np.array(pos);
    if deg:
        _pos[0], _pos[1] = np.deg2rad([_pos[0], _pos[1]]);
        
    x = pos[2] * np.cos(_pos[0])*np.cos(_pos[1]);
    y = pos[2] * np.cos(_pos[0])*np.sin(_pos[1]);
    z = pos[2] * np.sin(_pos[0]);
    
    return np.array([x, y, z]);

def closestVals(arr, curr):
    """
    Search array, return closest values to curr
    We can use this to find the next height above and below a ray
    Return (below, above)
    """
    idx = np.searchsorted(arr, curr);
    
    # If we've got a direct hit:
    if (arr[idx] == curr):
        if idx > 0: # Not at the very beginning:
            below = arr[idx - 1];
            if idx < len(arr):
                # In between layers
                above = arr[idx + 1];
            else:
                # We're at the last layer
                above = -1;
        else: # We're right at the beginning
            below = 0;
            above = 0 + 1;
        
        return below, above;
            
    # Else, find closest matches.
    below = arr[idx - 1] if idx > 0 else 0
    above = arr[idx] if idx < len(arr) else -1
    return below, above

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
    If no collision, Distance Traveled = None
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
        #a, b, x, y = np.deg2rad((sigPos[0, -1], newPos[0], 
        #                         sigPos[1, -1], newPos[1]));
        #dist = EARTH_R * np.arccos(np.cos(a)*np.cos(b)*np.cos(x-y) + np.sin(a)*np.sin(b));
        # Find distance between two points in Spherical Coords
        # Considering latitude as Theta, longitude as Phi, radius as R:
        theta1, phi1, theta2, phi2 = np.deg2rad((sigPos[0, -1], sigPos[1, -1],
                                                 newPos[0]    , newPos[1]));
        r1, r2 = np.add(EARTH_R, (sigPos[2, -1], newPos[2]));
        dist = np.sqrt(r1**2 + r2**2 
                       - 2*r1*r2*
                           (np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)
                           + np.cos(theta1)*np.cos(theta2)));
        # ^^^ Not 100% sure about this. Investigate in future.
        # Clamp Latitude and Logitude
        newPos = wrapLatLon(newPos);
        return (newPos, dist); # Collision occurs at newPos after travelling distTraveled
    else:
        return (signal['pos'], None); # No collision occurs.
    
def getRelativePermittivity(pos, freq, dt):
    """
        Return relative permittivity at a `pos1 [lat, lon, alt] at DateTime dt
        Uses IRI2016 to pull electron density at said point.
    """
    sim = iri.IRI(glat=pos[0], glon=pos[1], 
                  altkmrange=[pos[2], pos[2], 1], time=dt); # Wrapper requires min, max, step for altkmrange argument.
    e_density = sim['ne'].item(); # Pull the value itself.
    assert e_density > 0, ("Electron density cannot be < 0. Altitude of " +
                           str(pos[2]) + " likely out of range.");
    
    # Next, determine the index of refraction
    # https://www.waves.utoronto.ca/prof/svhum/ece422/notes/20c-ionosphere.pdf
    w_p = np.sqrt((e_density * (E**2)) / (M_E * EPSILON_0));    # Angular Plasma Freq
    w = 2*np.pi*freq;                                           # Signal Frequency
    epsilon_r = 1 - (w_p/w)**2;                                 # Relative Permittivity
    
    return epsilon_r;

def getRefractionCoeffs(eps_1, eps_2, theta_i, theta_t, deg=0):
    """
    Pull Reflection and Transmission Coefficients for given Relative Permittivities
    Assume Relative Permability is the same.
    https://core.ac.uk/download/pdf/130141591.pdf
    https://ocw.mit.edu/courses/6-007-electromagnetic-energy-from-motors-to-lasers-spring-2011/6e5c4b74150da32c3b98964855bbdccb_MIT6_007S11_lec29.pdf
    
    Assumes lossless medium

    Parameters
    ----------
    eps_1 : FLOAT
        Relative permittivity of medium from Ray Origin
    eps_2 : FLOAT
        Relative permittivity of medium from Ray Destination
    theta_i : np.array
        Incident Angles
    theta_t : np.array
        Transmitted Angles (found via Snell's law)
    deg: BOOLEAN
        1 if input angles are degrees. Assume radians by default.

    Returns
    -------
    refleCoeff, transCoeff - FLOATS. Complex magnitude.
        Multiply by Original Ray to pull Transmitted and Reflected Magnitudes
    """
    eta_1 = np.sqrt(MU_0/eps_1);
    eta_2 = np.sqrt(MU_0/eps_2);
    
    if(deg):
        theta_i, theta_t = np.deg2rad([theta_i, theta_t]);
    
    """
    # Parallel Component
    R_parallel = (eta_2 * np.cos(theta_i) - eta_1 * np.cos(theta_t)) / \
                 (eta_2 * np.cos(theta_i) + eta_1 * np.cos(theta_t))
    T_parallel = np.add(1, -1*R_parallel)  # 1 - R
    
    # Perpendicular Component, just in case.
    R_perp = (eta_2 * np.cos(theta_t) - eta_1 * np.cos(theta_i)) / \
             (eta_2 * np.cos(theta_i) + eta_1 * np.cos(theta_t))
    T_perp = np.add(1, -1*R_perp)  # 1 - R
    """
    R = (eta_2 - eta_1) / (eta_2 + eta_1);
    T = 2*eta_2 / (eta_2 + eta_1);

    return R, T
    
def refractSignal(signal, intersectPos, dt):
    """
        For a given signal with a magnitude, position [lat, lon, alt], and 
        direction along unit cosines with [lat, lon, alt], return the reflected
        and refracted signal when treating the Ionosphere at the `intersectPos`
        as a homogeneous medium.

    Parameters
    ----------
    signal : SIGNAL
        Dict with 'pos' and 'dir' arrays, each of which denote movement along
        the geodetic coordinate system: [Latitude, Longitude, Altitude]
        Assume 'pos' in km. 'dir' must be normalized.
        Possesses a 'mag' attribute for signal magnitude.
    intersectPos : Geodetic 'pos' array
        [Latitude, Longitude, Altitude]
    dt : DateTime object
        To feed into IRI2016 module for Ionospheric Measurements

    Returns
    -------
    signalReflect, signalRefract -- two signals, reflected and refracted.
    """
    # Obtain relative permittivity in previous medium:       
    if (signal['epsilon_r'] is None):
        epsilon_pre = getRelativePermittivity(signal['pos'], signal['freq'], dt);
    else:
        epsilon_pre = signal['epsilon_r'];
    n1 = np.sqrt(epsilon_pre); # Refraction Index. Assume mu constant btwn both media
    
    # Find Relative Permittivity of next medium:
    # Check bounds:
    if (intersectPos[2] < MIN_IRI_ALT): # Consider it the troposphere
        print("Last interface before re-entering troposphere");
        epsilon_new = EPSILON_R_ATM; # Approximate with air
    elif (intersectPos[2] > MAX_IRI_ALT): # Exiting the ionosphere, escaping to space
        raise Exception("Ray is escaping to space. This should not be resolved here.");
    else:                               # Propagating thru ionosphere
        epsilon_new = getRelativePermittivity(intersectPos, signal['freq'], dt);
    n2 = np.sqrt(epsilon_new);
        
    #pdb.set_trace() if epsilon_new < 0 else None
    """
    # Apply Snell's Law: https://eng.libretexts.org/Bookshelves/Materials_Science/Supplemental_Modules_(Materials_Science)/Optical_Properties/Snell%27s_Law
    # Define Incident Angles along Lat & Long Slices
    incidentLat, incidentLon    = np.deg2rad([signal['dir'][0], signal['dir'][1]]);
    # Apply Snell's Law to compute refracted angles
    refractLat = np.arcsin((n1 / n2) * np.sin(incidentLat));
    refractLon = np.arcsin((n1 / n2) * np.sin(incidentLon));
    refractAlt = np.cos(refractLat) + np.cos(refractLon); # Add normal components.    
    
    # Pull Reflection and Transmission Coefficients
    refleCoeff, transCoeff = getRefractionCoeffs(epsilon_pre, epsilon_new,   \
                                                 [incidentLat, incidentLon], \
                                                 [refractLat, refractLon]);
    
    # Package Refracted Direction
    refractLat, refractLon      = np.rad2deg([refractLat, refractLon]);
    refractDir = np.array([refractLat, refractLon, refractAlt]);
    refractDir = refractDir / np.linalg.norm(refractDir); # Normalize vector
    
    # Reflect along normal. Magnitude should be preserved.
    reflectDir = np.array([signal['dir'][0], signal['dir'][1], -1*signal['dir'][2]]); 
        
        
    # Update Position, Direction, and Magnitude
    signal['pos'] = np.append(signal['pos'], intersectPos.reshape(-1,1), axis=1);
    signalReflect = copy.deepcopy(signal);
    signalRefract = copy.deepcopy(signal);
    
    signalReflect['dir'] = reflectDir;
    reflectSignalMag     = signal['mag']*np.array([refleCoeff[0]*signal['dir'][0], \
                            refleCoeff[1]*signal['dir'][1], \
                            2*(refleCoeff[0] + refleCoeff[1])*signal['dir'][2]]);
    signalReflect['mag'] = np.linalg.norm(reflectSignalMag);
        
    signalRefract['dir'] = refractDir;
    refractSignalMag     = signal['mag']*np.array([transCoeff[0]*signal['dir'][0], \
                            transCoeff[1]*signal['dir'][1], \
                            2*(transCoeff[0] + transCoeff[1])*signal['dir'][2]]);
    signalRefract['mag'] = np.linalg.norm(refractSignalMag);
    """
    
    # Prepare Output Vars
    signal['pos'] = np.append(signal['pos'], intersectPos.reshape(-1,1), axis=1);
    signalReflect = copy.deepcopy(signal);
    signalRefract = copy.deepcopy(signal);
    
    # Apply Snell's Law: https://eng.libretexts.org/Bookshelves/Materials_Science/Supplemental_Modules_(Materials_Science)/Optical_Properties/Snell%27s_Law        
    # Reasoning behind using incident permittivity as Great Divider for Total Reflection:
    # https://www.waves.utoronto.ca/prof/svhum/ece422/notes/20c-ionosphere.pdf on pg3    
    if epsilon_new > 0: # Propagation constant is real, refract.
        # Define Incident Angles along Lat & Long Slices
        incidentLat, incidentLon    = np.deg2rad([signal['dir'][0], signal['dir'][1]]);
        # Apply Snell's Law to compute refracted angles
        refractLat = np.arcsin((n1 / n2) * np.sin(incidentLat));
        refractLon = np.arcsin((n1 / n2) * np.sin(incidentLon));
        refractAlt = np.sin(refractLat) + np.sin(refractLon);
    
        # Pull Reflection and Transmission Coefficients
        refleCoeff, transCoeff = getRefractionCoeffs(epsilon_pre, epsilon_new,   \
                                                     [incidentLat, incidentLon], \
                                                     [refractLat, refractLon]);
        """    
        refractSigMag = signal['mag']*np.linalg.norm(np.add(
                (transCoeff[0]**2)*polarToCartesian([refractLat, 0, 1]), \
                (transCoeff[1]**2)*polarToCartesian([0, refractLon, 1])
                ));
        reflectSigMag = signal['mag']*np.linalg.norm(np.add(
                (refleCoeff[0]**2)*polarToCartesian([incidentLat, 0, 1]), \
                (refleCoeff[1]**2)*polarToCartesian([0, incidentLon, 1])
                ));   
        """
        reflectSigMag = (transCoeff - 1)**2 * signal['mag'];
        refractSigMag = (1 - (transCoeff-1)**2) * signal['mag'];
            
        # Package Refracted Signal
        refractLat, refractLon      = np.rad2deg([refractLat, refractLon]);
        refractDir = np.array([refractLat, refractLon, refractAlt]);
        refractDir = refractDir / np.linalg.norm(refractDir); # Normalize vector
        
        signalRefract['dir'] = refractDir;
        signalRefract['mag'] = refractSigMag;
        signalRefract['epsilon_r'] = epsilon_new;
    else: # Relative Permittivity is negative, complex propagation -- totally reflected
        print("RefractSignal: Negative Permittivity");
        # Reflect completely
        reflectSigMag = signal['mag'];
        
        # Refracted Wave Does Not Exist    
        signalRefract['dir'] = [None, None, None];    
        signalRefract['mag'] = 0;
        signalRefract['epsilon_r'] = None;
    
    # Package Reflected Signal
    # Reflect along normal. Magnitude should be preserved.
    reflectDir = np.array([signal['dir'][0], signal['dir'][1], -1*signal['dir'][2]]); 

    # Update Position, Direction, and Magnitude
    signalReflect['dir'] = reflectDir;
    signalReflect['mag'] = reflectSigMag;
    signalReflect['epsilon_r'] = epsilon_pre;    
    
    print("\n POSITION: ", signalRefract['pos'][:, -1], "\nReflected Power ", signalReflect['mag'], " vs. Transmitted Power ", signalRefract['mag'],
          "\nRX DIR ", signalReflect['dir'], " vs. TX DIR ", signalRefract['dir']);
    
    return signalReflect, signalRefract;
    

def raytrace(TXprop, signal, altitudes, dt):
    """
        For a given Transmitter with properties `TXprop`, transmitting `signal` at 
        a time `dt`, we trace the signal until it either reaches space or reaches 
        the ground
    """
    
    # Pull Transmitting Position & Direction, apply to Initial Signal Info
    signal['dir'] = np.add(signal['dir'], TXprop['dir']);
    signal['pos'] = np.append(signal['pos'], TXprop['pos'], axis=1);
    
    # Simulate all given Ionospheric Layers, inidicated by `altitudes`
        # We can do this later in more complicated examples.

    # Begin Transmission
    # While Ray hasn't hit Ground or Space:
    while((signal['pos'][2, -1] > GND_HEIGHT) & (signal['pos'][2, -1] < MAX_IRI_ALT)):  
        # Determine where next Intersection occurs
        # Find the closest Intersection Point (Ionosphere Layer or Ground)
        alt_idx = np.searchsorted(altitudes, signal['pos'][2, -1]);
        if(signal['dir'][2] > 0):   # Gaining Altitude
            if(signal['pos'][2, -1] >= altitudes[-1]):
                # If we're escaping into space, quit
                print("Ray has escaped to space.");
                break;
            else:
                # Check to see if the closest layer is above or below
                if altitudes[alt_idx] > signal['pos'][2, -1]:
                    targetAlt = altitudes[alt_idx];
                elif alt_idx < len(altitudes) - 1: # We have space to move up
                    targetAlt = altitudes[alt_idx + 1];
                elif altitudes[alt_idx] == signal['pos'][2, -1]:
                    # We're directly at the altitude already
                    if alt_idx < len(altitudes) - 1:    # We have space to move up
                        targetAlt = altitudes[alt_idx + 1];
                    else:
                        print("Ray has escaped to space");
                        break;                          # Escaping to space
                    
        else:                       # Losing Altitude
            if(signal['pos'][2, -1] <= altitudes[0]):
                # If we're hitting the ground, the target is the ground.
                targetAlt = GND_HEIGHT;
                print("Ray to hit ground.");
            else:
                # Check to see if the closest layer is above or below
                if altitudes[alt_idx] < signal['pos'][2, -1]:
                    # There's wiggle room between the signal and the next lowest layer
                    targetAlt = altitudes[alt_idx];
                elif altitudes[alt_idx] == signal['pos'][2, -1]:
                    # We're directly at the altitude already.
                    if alt_idx > 1: # We have space to move down
                        targetAlt = altitudes[alt_idx - 1];
                    else:           # We're on the lowest layer
                        targetAlt = GND_HEIGHT; # Hitting the ground
                        print("Ray to hit ground");
            
        # Find intersection point and distance traveled
        newPos, dist = collisionPoint(signal, targetAlt);
        
        # If we're not hitting the ground:
        if targetAlt == GND_HEIGHT:
            # Now propagating back through the Atmosphere
            # TODO Calculate Path Loss...
            # Return Endpoint, Break out.
            signal['pos'] = np.append(signal['pos'], newPos.reshape(-1,1), axis=1);
            break;
        else:
            # Propagating through the ionosphere
            signalReflect, signalRefract = refractSignal(signal, newPos, dt);
            # Choose signal with the greatest magnitude
            signal = max([signalReflect, signalRefract], key=lambda x: x['mag']);
      
    
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