# Standard library imports #
import os

# 3rd party inputs #
import numpy as np
import pickle

# Project imports
from pyCLARIS import pyCLARISAnalysis as claris
from pyCLARIS import coastalGeoUtils as utils


pc_aft = claris.pcManager('C:/Users/conli/Documents/FRF/Data/20210323_FRF5km_postNEer/FRF.las')  
pc_bef = claris.pcManager('C:/Users/conli/Documents/FRF/Data/20210318_FRF5km/FRF.las')

#############################################
# Visual HWL rgb #
#############################################
##if not os.path.exists('C:/Users/conli/Documents/frf/pyCLARIS_project/Analyses/HWLComp/data/HWL_rgb.pkl'):
##HWL_rgb = pc_aft.grabFeatures(1,'rgb')
##HWL_rgb = HWL_rgb[0]
##else:
##    f = open('C:/Users/conli/Documents/frf/pyCLARIS_project/Analyses/HWLComp/data/HWL_rgb.pkl','rb')
##    HWL_rgb = pickle.load(f)

#############################################
# R2% #
#############################################
T = pc_bef.createTransects(dy=5)
hydro = utils.hydroLab(20210318,20210323,8651370,44056)
HWL_R2 = np.empty([0,3])
for t in T:
    beta = utils.calcTransectSlope(t[1]['X'],t[1]['Z'])
    TWL = hydro.calcTWL(beta)
    hwl_pos = hydro.transectTWLIntercept(max(TWL[:,1]),t[1]['X'],t[1]['Z'])
    HWL_R2 = np.vstack([HWL_R2,np.hstack([hwl_pos,t[1]['Y'][0],max(TWL[:,1])])])


#############################################
# Plot on image and difference figure
#############################################
xx1 = np.arange(50,120,0.25)
yy1 = np.arange(0,1000,0.25)
im = pc_bef.gridPC_Slocum(xx1,yy1,z_val='rgb',function='mean')
xx2 = np.arange(50,120,1)
yy2 = np.arange(0,1000,1)
dsm_bef = pc_bef.gridPC_Slocum(xx2,yy2,z_val='z',function='min')
dsm_aft = pc_aft.gridPC_Slocum(xx2,yy2,z_val='z',function='min')


def plotR2HWL():
    plt.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(6.5,4))
    ax1 = plt.axes([.12,.5,.83,.45])
    ax2 = plt.axes([.12,.1,.37,.3])
    ax3 = plt.axes([.58,.1,.37,.3])
     
    ax1.imshow(np.flipud(np.transpose(im,(1,0,2))),extent=(min(yy1),max(yy1),min(xx1),max(xx1)),aspect='auto')
    ax1.plot(HWL_R2[:,1],HWL_R2[:,0],'b.')
    ax1.invert_yaxis()
    ax1.plot((120,240),(75,75),'k')
    ax1.plot((120,240),(90,90),'k')
    ax1.plot((120,120),(75,90),'k')
    ax1.plot((240,240),(75,90),'k')
    ax1.plot((875,1000),(55,55),'k')
    ax1.plot((875,1000),(70,70),'k')
    ax1.plot((875,875),(55,70),'k')
    ax1.plot((1000,1000),(55,70),'k')
    ax1.set_xlim(0,1000)
    ax1.set_xlabel('FRF Y (m)',labelpad=0.1)
    ax1.set_ylabel('FRF X (m)')

    ax2.imshow(np.flipud(np.transpose(im,(1,0,2))),extent=(min(yy),max(yy),min(xx),max(xx)),aspect='auto')
    ax2.plot(HWL_R2[:,1],HWL_R2[:,0],'b.')
    ax2.invert_yaxis()
    ax2.set_xlim(120,240)
    ax2.set_ylim(90,75)

    ax3.imshow(np.flipud(np.transpose(im,(1,0,2))),extent=(min(yy),max(yy),min(xx),max(xx)),aspect='auto')
    ax3.plot(HWL_R2[:,1],HWL_R2[:,0],'b.')
    ax3.invert_yaxis()
    ax3.set_xlim(875,1000)
    ax3.set_ylim([70,55])
    ax3.yaxis.tick_right()

    con1 = ConnectionPatch(xyA=(180,90), xyB=(120,75), coordsA="data", coordsB="data",
                      axesA=ax1, axesB=ax2, color="k")
    con2 = ConnectionPatch(xyA=(180,90), xyB=(240,75), coordsA="data", coordsB="data",
                      axesA=ax1, axesB=ax2, color="k")
    con3 = ConnectionPatch(xyA=(937.5,70), xyB=(875,55), coordsA="data", coordsB="data",
                      axesA=ax1, axesB=ax3, color="k")
    con4 = ConnectionPatch(xyA=(937.5,70), xyB=(1000,55), coordsA="data", coordsB="data",
                      axesA=ax1, axesB=ax3, color="k")    
    ax1.add_artist(con1)
    ax1.add_artist(con2)
    ax1.add_artist(con3)
    ax1.add_artist(con4)
    
    fig.show()
    plt.savefig('C:/Users/conli/Documents/FRF/pyCLARIS_project/Analyses/HWLComp/figs/R2CalcFig.png',dpi=350)

    
def plotR2HWLOnDoD():
    
    fig = plt.figure(figsize=(6.5,2.5))
    ax = plt.axes([.1,.15,.8,.75])
    plt.rcParams.update({'font.size': 8})

    h = ax.pcolor(yy2,xx2,np.transpose(dsm_aft-dsm_bef),vmin=-1,vmax=1,cmap='seismic_r')
    ax.set_ylim(50,125)
    ax.invert_yaxis()
    ax.set_xlabel('FRF X (m)')
    ax.set_ylabel('FRF Y (m)')
    ax.text(50,120,'N')
    ax.arrow(70,118,20,0,width=1.5,head_width=5,color='k')
    ax.plot(HWL_R2[:,1],HWL_R2[:,0],'k.',markersize=1)

    cbax = plt.axes([.91,.15,.05,.75])
    cbax.axis('off') 
    fig.colorbar(h,label='$\Delta z$ (m)',orientation='vertical',fraction=1,format='%.1f')

    fig.show()
    plt.savefig('C:/Users/conli/Documents/FRF/pyCLARIS_project/Analyses/HWLComp/figs/R2HWLOnDoD.png',dpi=350)

    
def plotHWLs():
    plt.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(6.5,4))
    ax1 = plt.axes([.12,.5,.83,.45])
    ax2 = plt.axes([.12,.1,.83,.3])
     
    ax1.imshow(np.flipud(np.transpose(im,(1,0,2))),extent=(min(yy),max(yy),min(xx),max(xx)),aspect='auto')
    ax1.plot(HWL_rgb[:,1],HWL_rgb[:,0],'r.')
    ax1.plot(HWL_R2[:,1],HWL_R2[:,0],'b.')
    ax1.invert_yaxis()

    ax2.plot(HWL_rgb[:,1],HWL_rgb[:,2],'r.')
    ax2.plot(HWL_R2[:,1],HWL_R2[:,2],'b.')
    
    fig.show()





