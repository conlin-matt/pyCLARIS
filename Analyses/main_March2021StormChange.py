
# Standard library imports #
import os

# 3rd party inputs #
import matplotlib.pyplot as plt
import numpy as np

# Project imports
import pyCLARIS.pyCLARISAnalysis as claris

#############################################
# March 2021 storm surface-change analysis #
#############################################
dsms = list()
Is = list()
T = list()
direcs = ["C:/Users/conli/Documents/FRF/Data/20210318_FRF5km",'C:/Users/conli/Documents/FRF/Data/20210323_FRF5km_postNEer']
for direc in direcs:
    
    # Create FRF point cloud if it hasn't already been made
    if not os.path.exists(direc+'/FRF.las'):
        claris.createFRFLas(direc)   

    # Create the PC object #
    pc = claris.pcManager(direc+'/FRF.las') 

    # Grid the PC into a DSM and a DIM (digital intensity model)#
    xx = np.arange(50,120,1)
    yy = np.arange(0,1000,1)
    dsm,I = pc.gridPC_Slocum(xx,yy,function='min')

    # Pull xshore transects from the PC #
    transects = pc.createTransects(dy=5)

    # Save the DSM and transects #
    dsms.append(dsm)
    Is.append(I)
    T.append(transects)



# Make figures #
plotBothSurfacesAndChangeFig()
plotChangeFig()
plotChangeFigWithTransects()






def plotBothSurfacesAndChangeFig():
    fig = plt.figure(figsize=(6.5,4))
    ax = fig.subplots(3)
    plt.rcParams.update({'font.size': 8})
        
    count = -1
    titles = ['March 18 (pre-storm)','March 21 (post-storm)']
    for axis in ax[0:-1]:
        count+=1
        h = axis.pcolor(yy,xx,np.transpose(dsms[count]),vmin=0,vmax=8,cmap='gist_earth')
        axis.set_ylim(50,125)
        axis.invert_yaxis()
        axis.set_title(titles[count],fontsize=8,pad=0.1)
        axis.set_xticks([])

    h2 = ax[2].pcolor(yy,xx,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')
    ax[2].set_ylim(50,125)
    ax[2].invert_yaxis()
    ax[2].set_xlabel('FRF X (m)')
    ax[2].set_ylabel('FRF Y (m)')
    ax[2].set_title('Elevation change',fontsize=8,pad=0.1)

    cbax1 = plt.axes([.91,.55,.05,.25])
    cbax1.axis('off')
    fig.colorbar(h,ax=cbax1,label='Elevation (m)',orientation='vertical',fraction=1)

    cbax2 = plt.axes([.91,.1,.05,.25])
    cbax2.axis('off') 
    fig.colorbar(h2,label='$\Delta z$ (m)',orientation='vertical',fraction=1,format='%.1f')

    fig.show()
    plt.savefig('C:/Users/conli/Documents/FRF/Analyses/March2021Storm/bothSurfacesAndChangeFig.png',dpi=350)

def plotChangeFig():
    
    fig = plt.figure(figsize=(6.5,2.5))
    ax = plt.axes([.1,.15,.8,.75])
    plt.rcParams.update({'font.size': 8})

    h = ax.pcolor(yy,xx,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')
    ax.set_ylim(50,125)
    ax.invert_yaxis()
    ax.set_xlabel('FRF X (m)')
    ax.set_ylabel('FRF Y (m)')
    ax.text(50,120,'N')
    ax.arrow(70,118,20,0,width=1.5,head_width=5,color='k')

    cbax = plt.axes([.91,.15,.05,.75])
    cbax.axis('off') 
    fig.colorbar(h,label='$\Delta z$ (m)',orientation='vertical',fraction=1,format='%.1f')

    fig.show()
    plt.savefig('C:/Users/conli/Documents/FRF/Analyses/March2021Storm/changeFig.png',dpi=350)


def plotChangeFigWithTransects():

    plt.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(6.5,4))
    ax1 = plt.axes([.1,.55,.8,.4])
    ax2 = plt.axes([.1,.1,.35,.3])
    ax3 = plt.axes([.55,.1,.35,.3])

    h = ax1.pcolor(yy,xx,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')
    ax1.set_ylim(50,125)
    ax1.invert_yaxis()
    ax1.set_xlabel('FRF X (m)')
    ax1.set_ylabel('FRF Y (m)')
    ax1.text(50,120,'N')
    ax1.arrow(70,118,20,0,width=1.5,head_width=5,color='k')

    cbax = plt.axes([.91,.55,.05,.4])
    cbax.axis('off') 
    fig.colorbar(h,label='$\Delta z$ (m)',orientation='vertical',fraction=1,format='%.1f')

    hh1 = ax2.plot(T[0][60][1]['X'],T[0][60][1]['Z'],'k')
    hh2 = ax2.plot(T[1][60][1]['X'],T[1][60][1]['Z'],'--',color='gray')
    ax1.plot(T[1][60][1]['Y'],T[1][60][1]['X'],'k')
    ax1.text(T[1][60][1]['Y'][0],T[1][60][1]['X'][0]-2,"q")
    ax1.text(T[1][60][1]['Y'][0],T[1][60][1]['X'][-1]+4,"q'")
    ax2.text(T[1][60][1]['X'][0],8.9,"q")
    ax2.text(T[1][60][1]['X'][-1],8.9,"q'")
    ax2.set_xlabel('x-shore dist (m)')
    ax2.set_ylabel('Elevation (m)')
    ax2.legend([hh1[0],hh2[0]],['pre-storm','post-storm'])
    ax2.set_ylim(0,8.5)

    ax3.plot(T[0][196][1]['X'],T[0][196][1]['Z'],'k')
    ax3.plot(T[1][196][1]['X'],T[1][196][1]['Z'],'--',color='gray')
    ax1.plot(T[1][196][1]['Y'],T[1][196][1]['X'],'k')
    ax1.text(T[1][196][1]['Y'][0],T[1][196][1]['X'][0]-1,"r")
    ax1.text(T[1][196][1]['Y'][0],T[1][196][1]['X'][-1]+6,"r'")
    ax3.text(T[1][196][1]['X'][0],8.9,"r")
    ax3.text(T[1][196][1]['X'][-1],8.9,"r'")
    ax3.set_xlabel('x-shore dist (m)')
    ax3.set_ylabel('Elevation (m)')
    ax3.set_ylim(0,8.5)
    
    fig.show()
    plt.savefig('C:/Users/conli/Documents/FRF/Analyses/March2021Storm/changeFigWithTransects.png',dpi=350)
    


    
