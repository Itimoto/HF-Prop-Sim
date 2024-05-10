"""
Mapping Data from IRI 2016 onto Map
"""

import numpy as np
import matplotlib.pyplot as plt
from cmocean import cm# For nice color maps

from datetime import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iri2016.profile as iri

# USER INPUTS:
ALTITUDE = 410; # Example altitude
TIME = datetime(2020, 6, 1, 0, 0, 0);
STEP = 5.0; # FLOAT64 STEP
LATRANGE = (-90, 90, STEP); # FLOAT. (Start, Stop, Step)

# Grid of all latitudes and longitudes
lats = np.linspace(LATRANGE[0], LATRANGE[1], int((LATRANGE[1]-LATRANGE[0])/STEP));
lons = np.linspace(-180, 180, int((360)/STEP));
lat2D, lon2D = np.meshgrid(lats, lons);

# Set up Ionosphere Sim
# Initialize arrays to store ionospheric parameters
e_density = np.zeros_like(lon2D);

# Initialize Plot to PLot as soon as New Data is Available:
plt.figure();
ax = plt.axes(projection=ccrs.PlateCarree());
ax.coastlines();
ax.add_feature(cfeature.BORDERS);

plt.title("Ionosphere E Density at {}km".format(ALTITUDE));
heatmap = ax.contourf(lon2D, lat2D, e_density, cmap=cm.tempo);  # Initialize empty heatmap
cbar = plt.colorbar(heatmap, ax=ax, label='Electron Density (e-/cm^3)');
plt.pause(0.001);

for lon_ind in range(np.size(lons)):
    starttime = datetime.now();
    # Specify the longitude
    lon = lons[lon_ind];
    # Retrieve ionospheric profile for current longitude slice
    sim = iri.geoprofile(latrange=LATRANGE, glon=int(lon), altkm=ALTITUDE, time=TIME);
    # We pull the following vars from the IRI program
    # SIMOUT = ["ne", "Tn", "Ti", "Te", "nO+", 
    #           "nH+", "nHe+", "nO2+", "nNO+","nCI", "nN+"]
    e_density[lon_ind, :] = np.array(sim["ne"]).flatten(); # Removes the extra column
    
    duration = datetime.now() - starttime;
    print("Longitude ", lon, " complete. Average: ", np.mean(e_density[lon_ind, :]),
          " Time Elapsed: ", duration);
    
    # Update plot with new data
    heatmap.remove(); cbar.remove(); # Remove previous contour plot & colorbar
    heatmap = ax.contourf(lon2D, lat2D, e_density, cmap=cm.tempo);
    cbar = plt.colorbar(heatmap, ax=ax, label='Electron Density (e-/cm^3)');
    plt.pause(0.001); # Pause for a short time to update plot
    
plt.show(); # Show final plot after loop completes
