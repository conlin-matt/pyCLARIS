import os

# 3rd party inputs #
import numpy as np
import pickle

# Project imports
from pyCLARIS import pyCLARISAnalysis as claris
from pyCLARIS import coastalGeoUtils as utils


pc_aft = claris.pcManager('/Users/frfuser/Documents/pyCLARIS_project/data/20210323_FRF5km_postNEer/FRF.las')  
pc_bef = claris.pcManager('/Users/frfuser/Documents/pyCLARIS_project/data/20210318_FRF5km/FRF.las')


#############################################
# R2% #
#############################################
f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/data/T_5km.pkl','rb');T = pickle.load(f);T_bef = T[0]#T = pc_bef.createTransects(dy=5)
hydro = utils.hydroLab(20210318,20210323,8651370,44056)
HWL_R2 = np.empty([0,3])
print('Calculating TWLs')
for t in T_bef:
    beta = utils.calcTransectSlope(t['X'],t['Z'])
    TWL = hydro.calcTWL(beta)
    hwl_pos = hydro.transectTWLIntercept(max(TWL[:,1]),t['X'],t['Z'])
    HWL_R2 = np.vstack([HWL_R2,np.hstack([hwl_pos,t['Y'][0],max(TWL[:,1])])])


#############################################
# Plot on image and difference figure
#############################################
print('Creating dsms and ims')
xx1 = np.arange(-200,120,0.5)#xx1 = np.arange(50,120,0.25)
yy1 = np.arange(-1000,4000,0.5)#yy1 = np.arange(0,1000,0.25)
im = pc_bef.gridPC_Slocum(xx1,yy1,z_val='rgb',function='mean')
xx2 = np.arange(-200,120,1)#xx2 = np.arange(50,120,1)
yy2 = np.arange(-1000,4000,1)#yy2 = np.arange(0,1000,1)
dsm_bef = pc_bef.gridPC_Slocum(xx2,yy2,z_val='z',function='min')
dsm_aft = pc_aft.gridPC_Slocum(xx2,yy2,z_val='z',function='min')
dsms = [dsm_bef,dsm_aft]






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
    ax.set_xlim(-1000,4000)
    ax.invert_yaxis()
    ax.set_xlabel('FRF Y (m)')
    ax.set_ylabel('FRF X (m)')
    ax.plot(HWL_R2[:,1],HWL_R2[:,0],'k.',markersize=1)

    cbax = plt.axes([.91,.15,.05,.75])
    cbax.axis('off') 
    fig.colorbar(h,label='$\Delta z$ (m)',orientation='vertical',fraction=1,format='%.1f')

    fig.show()
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/HWLAnalysis/figs/R2HWLOnDoD.png',dpi=350)

    
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


def plotDEMsWithRegionsAndTransectsAndHWL_5km():

    plt.rcParams.update({'font.size': 8})

    fig = plt.figure(figsize=(10,7.5))   
    ax1 = plt.axes([.1,.85,.85,.12]);ax1.set_xticklabels([]);ax1.set_title('Pre-storm DSM',fontweight='normal',loc='left',pad=0,fontsize=8)
    ax2 = plt.axes([.1,.70,.85,.12]);ax2.set_xticklabels([]);ax2.set_title('Post-storm DSM',fontweight='normal',loc='left',pad=0,fontsize=8)
    ax3 = plt.axes([.1,.55,.85,.12]);ax3.set_title('Difference',fontweight='normal',loc='left',pad=0,fontsize=8)

    axt1_4 = plt.axes([.1,.1,.12,.1]);axt1_4.set_ylim(0,7.5)
    axt1_3 = plt.axes([.1,.2,.12,.1]);axt1_3.set_xticklabels([]);axt1_3.set_yticklabels([]);axt1_3.set_ylim(0,9)
    axt1_2 = plt.axes([.1,.3,.12,.1]);axt1_2.set_xticklabels([]);axt1_2.set_yticklabels([]);axt1_2.set_ylim(0,9)
    axt1_1 = plt.axes([.1,.4,.12,.1]);axt1_1.set_xticklabels([]);axt1_1.set_yticklabels([]);axt1_1.set_ylim(0,9)

    axt2_4 = plt.axes([.25,.1,.12,.1]);axt2_4.set_ylim(0,7.5);axt2_4.set_yticklabels([])
    axt2_3 = plt.axes([.25,.2,.12,.1]);axt2_3.set_xticklabels([]);axt2_3.set_yticklabels([]);axt2_3.set_ylim(0,9)
    axt2_2 = plt.axes([.25,.3,.12,.1]);axt2_2.set_xticklabels([]);axt2_2.set_yticklabels([]);axt2_2.set_ylim(0,9)
    axt2_1 = plt.axes([.25,.4,.12,.1]);axt2_1.set_xticklabels([]);axt2_1.set_yticklabels([]);axt2_1.set_ylim(0,9)  

    axt3_4 = plt.axes([.4,.1,.12,.1]);axt3_4.set_ylim(0,7.5)
    axt3_3 = plt.axes([.4,.2,.12,.1]);axt3_3.set_xticklabels([]);axt3_3.set_yticklabels([]);axt3_3.set_ylim(0,9)
    axt3_2 = plt.axes([.4,.3,.12,.1]);axt3_2.set_xticklabels([]);axt3_2.set_yticklabels([]);axt3_2.set_ylim(0,9)
    axt3_1 = plt.axes([.4,.4,.12,.1]);axt3_1.set_xticklabels([]);axt3_1.set_yticklabels([]);axt3_1.set_ylim(0,9)

    axt4_4 = plt.axes([.55,.1,.12,.1]);axt4_4.set_ylim(0,7.5)
    axt4_3 = plt.axes([.55,.2,.12,.1]);axt4_3.set_xticklabels([]);axt4_3.set_yticklabels([]);axt4_3.set_ylim(0,9)
    axt4_2 = plt.axes([.55,.3,.12,.1]);axt4_2.set_xticklabels([]);axt4_2.set_yticklabels([]);axt4_2.set_ylim(0,9)
    axt4_1 = plt.axes([.55,.4,.12,.1]);axt4_1.set_xticklabels([]);axt4_1.set_yticklabels([]);axt4_1.set_ylim(0,9)

    axt5_4 = plt.axes([.7,.1,.12,.1]);axt5_4.set_ylim(0,7.5)
    axt5_3 = plt.axes([.7,.2,.12,.1]);axt5_3.set_xticklabels([]);axt5_3.set_yticklabels([]);axt5_3.set_ylim(0,9)
    axt5_2 = plt.axes([.7,.3,.12,.1]);axt5_2.set_xticklabels([]);axt5_2.set_yticklabels([]);axt5_2.set_ylim(0,9)
    axt5_1 = plt.axes([.7,.4,.12,.1]);axt5_1.set_xticklabels([]);axt5_1.set_yticklabels([]);axt5_1.set_ylim(0,9)

    axt6_4 = plt.axes([.85,.1,.12,.1]);axt6_4.set_ylim(0,7.5)
    axt6_3 = plt.axes([.85,.2,.12,.1]);axt6_3.set_xticklabels([]);axt6_3.set_yticklabels([]);axt6_3.set_ylim(0,9)
    axt6_2 = plt.axes([.85,.3,.12,.1]);axt6_2.set_xticklabels([]);axt6_2.set_yticklabels([]);axt6_2.set_ylim(0,9)
    axt6_1 = plt.axes([.85,.4,.12,.1]);axt6_1.set_xticklabels([]);axt6_1.set_yticklabels([]);axt6_1.set_ylim(0,9) 

    h = ax1.contourf(yy2,xx2,np.transpose(dsms[0]),levels=[0,1,2,3,4,5,6,7,8],vmin=0,vmax=8,cmap='gist_earth')
    ax1.set_ylim(-75,125)
    ax1.set_xlim(-1000,4000)
    ax1.invert_yaxis()

    ax2.contourf(yy2,xx2,np.transpose(dsms[1]),levels=[0,1,2,3,4,5,6,7,8],vmin=0,vmax=8,cmap='gist_earth')
    ax2.set_ylim(-75,125)
    ax2.set_xlim(-1000,4000)
    ax2.invert_yaxis()

    ax3.pcolor(yy2,xx2,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')    
    ax3.set_ylim(-75,125)
    ax3.set_xlim(-1000,4000)
    ax3.invert_yaxis()
##    for region in regions_agg:
##        ax3.plot(region[:,1],region[:,0],'c')
##    for region in regions_deg:
##        ax3.plot(region[:,1],region[:,0],'m')


    T_nums = [ [185,186,191,196],[350,355,360,363],[475,480,485,490],[650,655,660,665],[780,785,790,795],[920,925,930,935] ]
    for i in range(0,6):
        for t in range(0,len(T_nums[i])):
            axx = eval('axt'+str(i+1)+'_'+str(t+1))
            axx.plot(T[0][T_nums[i][t]]['X'],T[0][T_nums[i][t]]['Z'],'k-')
            axx.plot(T[1][T_nums[i][t]]['X'],T[1][T_nums[i][t]]['Z'],'--',color='grey')
            axx.plot((HWL_R2[T_nums[i][t],0],HWL_R2[T_nums[i][t],0]),axx.get_ylim(),'b',linewidth=0.75)
            axx.set_ylim(0,9)
            ax3.plot(T[0][T_nums[i][t]]['Y'],T[0][T_nums[i][t]]['X'],'k',linewidth=0.5)            


    fig.show()
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/HWLAnalysis/figs/changesOverview_5km_withHWL.png',dpi=350)    
