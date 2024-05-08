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

## Setting up a map of Earth: ##-----------------------------------------------
plt.axes(projection=ccrs.PlateCarree())
plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()

ax.add_feature(cartopy.feature.OCEAN)
ax.add_feature(cartopy.feature.LAND, edgecolor='black')
ax.add_feature(cartopy.feature.LAKES, edgecolor='black')

#### Adding Data to the Map ####-----------------------------------------------
# Sample Data Points. Draw a line from NY to Honolulu
new_york = dict(lon=-74.0060, lat=40.7128)
honolulu = dict(lon=-157.8583, lat=21.3609)
lons = [new_york['lon'], honolulu['lon']]
lats = [new_york['lat'], honolulu['lat']]

ax.plot(lons, lats, label='Equirectangular straight line');
ax.plot(lons, lats, label='Great Circle', transform=ccrs.Geodetic());
# ^^^ `transform` arg tells Cartopy what coordinate system data is defined in
ax.legend();
ax.set_global(); # <-- This zooms it out
