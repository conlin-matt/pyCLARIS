# pyCLARIS

## Background ##

pyCLARIS is a package to facilitate and complete analysis of CLARIS data. pyCLARIS was built with the specific goal of enabling analysis of point clouds with tens of millions of points quickly and effeciently on a regular laptop. pyCLARIS contains the following functionality:
- A function to transform a .las point cloud file into FRF coordinates, and re-save it in those coordinates
- A point cloud manager class which contains the following methods:
    - view the point cloud
    - manually select and save points from the point cloud
    - grid (quickly) the point cloud into a DSM or an interpolated image
    - extract cross-shore transects from the point cloud
    - load the point cloud into memory as a variable. 
- A function to automatically extract contiguous regions of elevation change from two DSMs (i.e. pre- and post-storm DSMs)
- A scarpManager class to analyze beach scarp changes. Work in progress.

pyCLARIS also contains top-level script files in the Analysis subfolders to reproduce analyses completed by Matt Conlin.

## Installation ##
pyCLARIS can be downloaded and used fairly simply using Anaconda and pip, with steps as follows:

1. Create a new environment named claris_env and put Python v3.6 in it:

    `conda create -n claris_env python=3.6`  
    `conda activate claris_env`
    
2. Install the PDAL package from the conda-forge channel:

    `conda install -c conda-forge pdal python-pdal gdal`
    
3. Install the package and most of its dependencies:

    `pip install git+https://github.com/conlin-matt/pyCLARIS`
    
4. Install the last package and prevent it from upgrading anything that would cause issues:

    `pip install --upgrade --upgrade-strategy only-if-needed pybeach`
    

## Usage ##
```python
import numpy as np
from pyCLARIS import pyCLARISAnalysis as claris

xx = np.arange(50,120,1)
yy = np.arange(0,1000,1)
        
# Create FRF point cloud if it hasn't already been made
claris.createFRFLas("/direc_to_pre-storm_files",croper='frf')   

# Create the PC object #
pc = claris.pcManager("/direc_to_pre-storm_files"+'/FRF.las') 

# View the PC #
pc.viewPC('rgb')

# Grid the PC into a DSM and a colored image #
dsm = pc.gridPC_Slocum(xx,yy,z_val='z',function='min')
im = pc.gridPC_Slocum(xx,yy,z_val='rgb',function='mean')

# Pull xshore transects from the PC #
transects = pc.createTransects(dy=5)   
```


