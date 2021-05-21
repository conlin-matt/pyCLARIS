# Standard library imports #
import os

# 3rd party inputs #
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from scipy.interpolate import interp1d

# Project imports
import pyCLARIS.pyCLARISAnalysis as claris


#############################################
# Examine and quantify changes to scarp at southern end of FRF property #
#############################################

# Manually extract scarp positions #
direcs = ["C:/Users/conli/Documents/FRF/Data/20210318_FRF5km",'C:/Users/conli/Documents/FRF/Data/20210323_FRF5km_postNEer']
scarps = []
for direc in direcs:
    
    # Create FRF point cloud if it hasn't already been made
    if not os.path.exists(direc+'/FRF.las'):
        claris.createFRFLas(direc)   

    # Create the PC object #
    pc = claris.pcManager(direc+'/FRF.las')
    scarp = pc.grabFeatures(1,'rgb')

    scarps.append(scarp)

# Calculate changes in x and z scarp position #
yys = np.arange(min(scarps[0][0][:,1]),max(scarps[0][0][:,1]),0.1)
xis = []
zis = []
for scarp in scarps:
    f = interp1d(scarp[0][:,1],scarp[0][:,0],bounds_error=False,fill_value=np.nan)
    xi = f(yys)

    f = interp1d(scarp[0][:,1],scarp[0][:,2],bounds_error=False,fill_value=np.nan)
    zi = f(yys)

    xis.append(xi)
    zis.append(zi)

# Plot #
plotScarpPositions()





def plotScarpPositions():    

    fig = plt.figure(figsize=(4.5,2.5))
    ax = plt.axes([.1,.15,.8,.75],projection='3d')
    plt.rcParams.update({'font.size': 8})

    ax.plot3D(scarps[0][0][:,0],scarps[0][0][:,1],scarps[0][0][:,2],'k.')
    ax.plot3D(scarps[1][0][:,0],scarps[1][0][:,1],scarps[1][0][:,2],'.',color='grey')
    ax.view_init(40,-70)

    ax.set_yticks(np.arange(50,250,50))

    ax.set_xlabel('FRF X (m)')
    ax.set_ylabel('FRF Y (m)')
    ax.set_zlabel('Elev. (m)')

    fig.legend(['Pre-storm','Post-storm'],handletextpad=.05,fontsize=8)

    fig.show()
    plt.savefig('C:/Users/conli/Documents/FRF/Analyses/March2021Storm/scarpPositionsFig.png',dpi=350)

