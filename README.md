# HF Prop Sim
 Modular HF Simulator

## Setting Up
The simulator was developed in Spyder, with Python, on a Windows machine.

(1) Setting up the Spyder Environment
    - Install Spyder
        - https://docs.spyder-ide.org/current/installation.html
        - This is a good start. But, we still need more libraries.
    - Install Miniconda
        - https://docs.anaconda.com/free/miniconda/miniconda-install/
        - This'll let us install the necessary libraries.
    - Initialize the Conda environment used for development, in the Git Repo
        - https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
        - Install necessary libraries
    - Point the Spyder Python Interpreter to the Conda environment interpreter
        - Tools > Preferences > Python Interpreter > Use the following Python Interpreter:
        - `C:\Users\...\HF-Prop-Sim\envs\python.exe`
        - **From this point, `EarthMap_PlusData.py` should be able to run**
    - Set up Fortran Compiler for IRI2016
        - https://www.scivision.dev/install-msys2-windows/
        - After installing, from within the MSYS2 prompt:
            - `pacman -S mingw-w64-ucrt-x86_64-gcc`
            - `pacman -S mingw-w64-x86_64-gcc-fortran`
            - Locate CMake directory. Suppose it's at `/mingw64/bin/cmake.exe`
            - `nano ~/.bashrc`
            - `export PATH=/mingw64/bin/cmake.exe`
        - Then, go to System Properties in Windows.
            - Add `C:\msys64\mingw64\bin` to Path
        - Restart Spyder
        **From this point, `IRI2016_Sample.py` should be able to run**