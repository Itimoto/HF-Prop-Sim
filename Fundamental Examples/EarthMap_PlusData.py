# -*- coding: utf-8 -*-
"""
Quick Earth Map Example

Python 3.8.10

https://earth-env-data-science.github.io/lectures/mapping_cartopy.html
"""

# Plotting Libraries
import matplotlib.pyplot as plt
# Map Libararies
import cartopy.crs as ccrs
import cartopy

import numpy as np # Always helpful for linspace & other MATLAB-adjacent funcs

## Plotting 2D (Raster) Data ##------------------------------------------------
lon = np.linspace(-80, 80, 25);
lat = np.linspace(30, 70, 25);
lon2d, lat2d = np.meshgrid(lon, lat);
data = np.cos(np.deg2rad(lat2d) * 4) + np.sin(np.deg2rad(lon2d) * 4);

plt.figure();
ax = plt.axes();
ax.contourf(lon2d, lat2d, data); # <-- Map Directly on Map
# Note: No transformations are necessary for PlateCarree projection

## Setting up a map of Earth: ##-----------------------------------------------
plt.figure();
ax = plt.axes(projection=ccrs.Orthographic());

ax.coastlines()
ax.add_feature(cartopy.feature.OCEAN)
ax.add_feature(cartopy.feature.LAND, edgecolor='black')
ax.add_feature(cartopy.feature.LAKES, edgecolor='black')

#### Adding Sample Data to the Map ####-----------------------------------------------
# Sample Data Points. Draw a line from NY to Honolulu
'''
new_york = dict(lon=-74.0060, lat=40.7128)
honolulu = dict(lon=-157.8583, lat=21.3609)
lons = [new_york['lon'], honolulu['lon']]
lats = [new_york['lat'], honolulu['lat']]

ax.plot(lons, lats, label='Equirectangular straight line');
ax.plot(lons, lats, label='Great Circle', transform=ccrs.Geodetic());
# ^^^ `transform` arg tells Cartopy what coordinate system data is defined in
'''
ax.legend();
ax.set_global(); # <-- This zooms it out

## Adding Our Raster Data to the Map ##----------------------------------------
ax.contourf(lon, lat, data, transform=ccrs.PlateCarree());
# ^^^ I guess... this works?

