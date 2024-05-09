"""
Mapping Data from IRI 2016 onto Map
"""

import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iri2016.profile as iri

# USER INPUTS:
ALTITUDE = 410; # Example altitude
TIME = datetime(2020, 6, 1, 0, 0, 0);
STEP = 10.0; # FLOAT64 STEP
LATRANGE = (-90, 90, STEP); # FLOAT. (Start, Stop, Step)

# Grid of all latitudes and longitudes
lats = np.linspace(LATRANGE[0], LATRANGE[1], int((LATRANGE[1]-LATRANGE[0])/STEP));
lons = np.linspace(-180, 180, int((360)/STEP));
lat2D, lon2D = np.meshgrid(lats, lons);

# Initialize arrays to store ionospheric parameters
e_density = np.zeros_like(lon2D);

for lon_ind in range(np.size(lons)):
    # Specify the longitude
    lon = lons[lon_ind];
    # Retrieve ionospheric profile for current longitude slice
    sim = iri.geoprofile(latrange=LATRANGE, glon=int(lon), altkm=ALTITUDE, time=TIME);
    # We pull the following vars from the IRI program
    # SIMOUT = ["ne", "Tn", "Ti", "Te", "nO+", 
    #           "nH+", "nHe+", "nO2+", "nNO+","nCI", "nN+"]
    e_density[lon_ind, :] = np.array(sim["ne"]).flatten(); # Removes the extra column
    print("Longitude ", lon, " complete. Average: ", np.mean(e_density[lon_ind, :]));
    
# Plot
plt.figure();
ax = plt.axes(projection=ccrs.PlateCarree());
ax.coastlines();
ax.add_feature(cfeature.BORDERS);

heatmap = ax.contourf(lon2D, lat2D, e_density, cmap='jet');
plt.colorbar(heatmap, ax=ax, label='Electron Density (e-/cm^3)');
plt.title("Ionosphere E Density at {}km".format(ALTITUDE));
plt.show();