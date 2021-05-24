# pyCLARIS
This is a package containing routines to analyze CLARIS data as well as applied analyses. Top-level script files to reproduce analyses are in the Analyses folder.

The package can be downloaded and used fairly simply using Anaconda
and pip, with steps as follows:

1. Create a new environment named claris_env and put Python v3.6 in it:

    `conda create -n claris_env python=3.6`
    
2. Install the PDAL package from the conda-forge channel:

    `conda install -c conda-forge pdal python-pdal gdal`
    
3. Install the rest of the package:

    `pip install git+https://github.com/conlin-matt/pyCLARIS`


