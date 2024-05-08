The conda environment is set up here for ease of access
However, since the environment is outside of the default envs folder:
1) Conda can no longer find your environment with the `--name` flag. 
   You'll generally need to pass the `--prefix` flag along with the environment's full
   path to find the environment.
2) Specifying an install path when creating your conda environments makes it
   so that your cmd prompt is now prefixed with active environment's absolute 
   path rather than the environment's name.

Navigate to the main directory, then:
`conda activate ./envs`

Note To Self: Using miniconda3 terminal.
To install libraries for the Python install, activate the Conda environment, then:
`conda install {package}`