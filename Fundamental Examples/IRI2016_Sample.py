# -*- coding: utf-8 -*-
"""
Pulling Data from IRI 2016

Just... some sample data!

We'll be using the following:
    https://github.com/space-physics/iri2016
    
Note:
    - This Python wrapper of IRI2016 uses the build-on-run technique
    - On the first run or `iri2016.IRI()`, the Fortran code is built
    - You need a Fortran compiler to run this. According to the README:
        - Linux: `apt install gfortran`
        - Max: `brew install gcc`
        - Windows: MSYS2
            - https://www.scivision.dev/install-msys2-windows/
            - Specifically, installing MSYS2. Then, within MSYS2, 
            `pacman -S mingw-w64-ucrt-x86_64-gcc`
            `pacman -S mingw-w64-x86_64-gcc-fortran`
            locate CMake directory. Suppose it's at `/mingw64/bin/cmake.exe`
            `nano ~/.bashrc`
            `export PATH=/mingw64/bin/cmake.exe`
            Then, go to System Properties in Windows, add `C:\msys64\mingw64\bin` to Path
            Restart Spyder, everything should work from there.
    
https://github.com/space-physics/iri2016/blob/main/Examples/Plasma%20Variables.ipynb
"""

import iri2016.profile as iri
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

time_start_stop = (datetime(2020, 6, 1, 0, 0, 0), datetime(2020, 6, 2));
time_step = timedelta(minutes=30);
alt_km_range = (100, 500, 10.);
glat = 65;
glon = -147.5;

# Info on the timeprofile() method here:
# https://github.com/space-physics/iri2016/blob/main/src/iri2016/profile.py
sim = iri.timeprofile(time_start_stop, time_step, alt_km_range, glat, glon);
print(sim);

ax = plt.figure().gca();
ax.pcolormesh(sim.time, sim.alt_km, sim.ne.T, shading="nearest");
ax.set_title("Number Density");
ax.set_xlabel("time [UTC]");
ax.set_ylabel("altitude [km]");