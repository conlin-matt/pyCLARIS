
# Standard library imports #
import os

# 3rd party inputs #
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pickle
import numpy as np
from scipy.interpolate import interp1d

# Project imports
import pyCLARIS.pyCLARISAnalysis as claris

#############################################
# Surface change analysis #
#############################################

analysisLen = 5
createdVars = True

if analysisLen == 1:

    dsms = list()
    ims = list()
    T = list()
    direcs = ["/Users/frfuser/Documents/pyCLARIS_project/data/20210318_FRF5km","/Users/frfuser/Documents/pyCLARIS_project/data/20210323_FRF5km_postNEer"]
    for direc in direcs:
        
        # Create FRF point cloud if it hasn't already been made
        if not os.path.exists(direc+'/FRF.las'):
            claris.createFRFLas(direc,croper='frf')   

        # Create the PC object #
        pc = claris.pcManager(direc+'/FRF.las') 

        # Grid the PC into a DSM and a colored image #
        xx1 = np.arange(50,120,1)
        yy1 = np.arange(0,1000,1)
        print('Making DSM')
        dsm = pc.gridPC_Slocum(xx1,yy1,z_val='z',function='min')
        xx2 = np.arange(50,120,0.25)#xx2 = np.arange(-200,120,0.5)
        yy2 = np.arange(0,1000,0.25)#yy2 = np.arange(-1000,4000,0.5)
        print('Making image')
        im = pc.gridPC_Slocum(xx2,yy2,z_val='rgb',function='mean')

        # Pull xshore transects from the PC #
        print('Making transects')
        transects = pc.createTransects(dy=5)

        # Save the DSM and transects #
        dsms.append(dsm)
        ims.append(im)
        T.append(transects)


    print('Extracting regions')
    # Extract areas of change and plot on pre-storm DEM and RGB image #
    ##regions_agg,regions_deg = claris.extractChangeAreas(xx1,yy1,dsms[0],dsms[1],thresh=0.1)
    
elif analysisLen == 5:

    if createdVars:
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/data/dsms_5km.pkl','rb')
        dsms = pickle.load(f)
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/data/ims_5km.pkl','rb')
        ims = pickle.load(f)
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/data/T_5km.pkl','rb')
        T = pickle.load(f)
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/data/regions_agg_5km.pkl','rb')
        regions_agg = pickle.load(f)
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/data/regions_deg_5km.pkl','rb')
        regions_deg = pickle.load(f)        
        xx1 = np.arange(-200,120,1)
        yy1 = np.arange(-1000,4000,1)
        xx2 = np.arange(-200,120,0.5)
        yy2 = np.arange(-1000,4000,0.5)
        

    else:
        dsms = list()
        ims = list()
        T = list()
        direcs = ["/Users/frfuser/Documents/pyCLARIS_project/data/20210318_FRF5km","/Users/frfuser/Documents/pyCLARIS_project/data/20210323_FRF5km_postNEer"]
        for direc in direcs:
            
            # Create FRF point cloud if it hasn't already been made
            if not os.path.exists(direc+'/FRF.las'):
                claris.createFRFLas(direc,croper='5km')   

            # Create the PC object #
            pc = claris.pcManager(direc+'/FRF.las') 

            # Grid the PC into a DSM and a colored image #
            xx1 = np.arange(-200,120,1)
            yy1 = np.arange(-1000,4000,1)
            print('Making DSM')
            dsm = pc.gridPC_Slocum(xx1,yy1,z_val='z',function='min')
            xx2 = np.arange(-200,120,0.5)
            yy2 = np.arange(-1000,4000,0.5)
            print('Making image')
            im = pc.gridPC_Slocum(xx2,yy2,z_val='rgb',function='mean')

            # Pull xshore transects from the PC #
            print('Making transects')
            transects = pc.createTransects(dy=5)

            # Save the DSM and transects #
            dsms.append(dsm)
            ims.append(im)
            T.append(transects)


            print('Extracting regions')
            # Extract areas of change and plot on pre-storm DEM and RGB image #
            regions_agg,regions_deg = claris.extractChangeAreas(xx1,yy1,dsms[0],dsms[1],thresh=0.1)


 
### Make figures #
##plotBothSurfacesAndChangeFig()
##plotChangeFig()
##plotChangeFigWithTransects()
##plotChangeFigWithRegions()
##plotRegions()











def plotBothSurfacesAndChangeFig():
    fig = plt.figure(figsize=(6.5,4))
    ax = fig.subplots(3)
    plt.rcParams.update({'font.size': 8})
        
    count = -1
    titles = ['March 18 (pre-storm)','March 21 (post-storm)']
    for axis in ax[0:-1]:
        count+=1
        h = axis.pcolor(yy1,xx1,np.transpose(dsms[count]),vmin=0,vmax=8,cmap='gist_earth')
        axis.set_ylim(-75,125)
        axis.set_xlim(-1000,4000)
        axis.invert_yaxis()
        axis.set_title(titles[count],fontsize=8,pad=0.1)
        axis.set_xticks([])

    h2 = ax[2].pcolor(yy1,xx1,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')
    ax[2].set_ylim(-75,125)
    ax[2].set_xlim(-1000,4000)
    ax[2].invert_yaxis()
    ax[2].set_xlabel('FRF Y (m)')
    ax[2].set_ylabel('FRF X (m)')
    ax[2].set_title('Elevation change',fontsize=8,pad=0.1)

    cbax1 = plt.axes([.91,.55,.05,.25])
    cbax1.axis('off')
    fig.colorbar(h,ax=cbax1,label='Elevation (m)',orientation='vertical',fraction=1)

    cbax2 = plt.axes([.91,.1,.05,.25])
    cbax2.axis('off') 
    fig.colorbar(h2,label='$\Delta z$ (m)',orientation='vertical',fraction=1,format='%.1f')

    fig.show()
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/figs/bothSurfacesAndChangeFig.png',dpi=350)

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
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/figs/changeFig.png',dpi=350)


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
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/figs/changeFigWithTransects.png',dpi=350)
    
def plotChangeFigWithRegions():
    
    fig = plt.figure(figsize=(6.5,2.5))
    ax = plt.axes([.1,.15,.8,.75])
    plt.rcParams.update({'font.size': 8})

    h = ax.pcolor(yy1,xx1,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')
    for region in regions_agg:
        ax.plot(region[:,1],region[:,0],'c')
    for region in regions_deg:
        ax.plot(region[:,1],region[:,0],'m')
    ax.set_ylim(-75,125)#ax.set_ylim(50,125)
    ax.invert_yaxis()
    ax.set_xlabel('FRF Y (m)')
    ax.set_ylabel('FRF X (m)')
##    ax.text(-910,120,'N')
##    ax.arrow(-800,118,200,0,width=3,head_width=5,color='k')

    cbax = plt.axes([.91,.15,.05,.75])
    cbax.axis('off') 
    fig.colorbar(h,label='$\Delta z$ (m)',orientation='vertical',fraction=1,format='%.1f')

    fig.show()
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/figs/changeFigWithRegions.png',dpi=350)


def plotRegions():

    plt.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(6.5,8))
    ax1 = plt.axes([.1,.69,.8,.22])
    ax2 = plt.axes([.1,.42,.8,.22])
    ax3 = plt.axes([.1,.15,.8,.22])

    ax1.imshow(np.flipud(np.transpose(ims[0],(1,0,2))),extent=(min(yy2),max(yy2),min(xx2),max(xx2)),aspect='auto')
    ax1.invert_yaxis()
    for region in regions_agg:
        ax1.plot(region[:,1],region[:,0],'c')
    for region in regions_deg:
        ax1.plot(region[:,1],region[:,0],'m')

    h = ax2.contourf(yy1,xx1,np.transpose(dsms[0]),levels=[0,1,2,3,4,5,6,7,8],vmin=0,vmax=8,cmap='gist_earth')
##    plt.contour(yy,xx,np.transpose(dsms[0]),levels=8,colors='grey')
    ax2.set_ylim(50,125)
    ax2.invert_yaxis()
    for region in regions_agg:
        ax2.plot(region[:,1],region[:,0],'c')
    for region in regions_deg:
        ax2.plot(region[:,1],region[:,0],'m')

    
    ax3.contourf(yy,xx,np.transpose(dsms[1]),levels=[0,1,2,3,4,5,6,7,8],vmin=0,vmax=8,cmap='gist_earth')
##    plt.contour(yy,xx,np.transpose(dsms[0]),levels=8,colors='grey')
    ax3.set_ylim(50,125)
    ax3.invert_yaxis()

    cbax = plt.axes([.91,.25,.05,.35])
    cbax.axis('off') 
    fig.colorbar(h,label='Elev. (m)',orientation='vertical',fraction=1,format='%.1f')

    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax1.set_title('Pre-storm image',fontsize=8,fontweight='normal',loc='left',pad=.1)
    ax2.set_title('Pre-storm DTM',fontsize=8,fontweight='normal',loc='left',pad=.1)
    ax3.set_title('Post-storm DTM',fontsize=8,fontweight='normal',loc='left',pad=.1)
    ax3.set_xlabel('FRF Y (m)')
    ax3.set_ylabel('FRF X (m)')
    ax3.text(50,120,'N')
    ax3.arrow(70,118,20,0,width=1.5,head_width=5,color='k')

    fig.show()
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/figs/regionsFigOverview.png',dpi=350)


def plotDEMsWithRegionsAndTransects():

    plt.rcParams.update({'font.size': 8})

    fig = plt.figure(figsize=(6.5,9))   
    ax1 = plt.axes([.1,.78,.85,.19]);ax1.set_xticklabels([]);ax1.set_title('Pre-storm DSM',fontweight='normal',loc='left',pad=0,fontsize=8)
    ax2 = plt.axes([.1,.55,.85,.19]);ax2.set_title('Post-storm DSM',fontweight='normal',loc='left',pad=0,fontsize=8)

    axt1_4 = plt.axes([.1,.1,.25,.1]);axt1_4.set_ylim(0,7.5)
    axt1_3 = plt.axes([.1,.2,.25,.1]);axt1_3.set_xticklabels([]);axt1_3.set_yticklabels([]);axt1_3.set_ylim(0,9)
    axt1_2 = plt.axes([.1,.3,.25,.1]);axt1_2.set_xticklabels([]);axt1_2.set_yticklabels([]);axt1_2.set_ylim(0,9)
    axt1_1 = plt.axes([.1,.4,.25,.1]);axt1_1.set_xticklabels([]);axt1_1.set_yticklabels([]);axt1_1.set_ylim(0,9)

    axt2_4 = plt.axes([.4,.1,.25,.1]);axt2_4.set_ylim(0,7.5);axt2_4.set_yticklabels([])
    axt2_3 = plt.axes([.4,.2,.25,.1]);axt2_3.set_xticklabels([]);axt2_3.set_yticklabels([]);axt2_3.set_ylim(0,9)
    axt2_2 = plt.axes([.4,.3,.25,.1]);axt2_2.set_xticklabels([]);axt2_2.set_yticklabels([]);axt2_2.set_ylim(0,9)
    axt2_1 = plt.axes([.4,.4,.25,.1]);axt2_1.set_xticklabels([]);axt2_1.set_yticklabels([]);axt2_1.set_ylim(0,9)  

    axt3_4 = plt.axes([.7,.1,.25,.1]);axt3_4.set_ylim(0,7.5)
    axt3_3 = plt.axes([.7,.2,.25,.1]);axt3_3.set_xticklabels([]);axt3_3.set_yticklabels([]);axt3_3.set_ylim(0,9)
    axt3_2 = plt.axes([.7,.3,.25,.1]);axt3_2.set_xticklabels([]);axt3_2.set_yticklabels([]);axt3_2.set_ylim(0,9)
    axt3_1 = plt.axes([.7,.4,.25,.1]);axt3_1.set_xticklabels([]);axt3_1.set_yticklabels([]);axt3_1.set_ylim(0,9)  

    h = ax1.contourf(yy1,xx1,np.transpose(dsms[0]),levels=[0,1,2,3,4,5,6,7,8],vmin=0,vmax=8,cmap='gist_earth')
    ax1.set_ylim(50,125)
    ax1.invert_yaxis()
##    for region in regions_agg:
##        ax1.plot(region[:,1],region[:,0],'c')
##    for region in regions_deg:
##        ax1.plot(region[:,1],region[:,0],'m')

    ax2.contourf(yy1,xx1,np.transpose(dsms[1]),levels=[0,1,2,3,4,5,6,7,8],vmin=0,vmax=8,cmap='gist_earth')
    ax2.set_ylim(50,125)
    ax2.invert_yaxis()

    T_nums = [25,26,31,36]
    for t in range(0,len(T_nums)):
        axx = eval('axt1_'+str(t+1))
        axx.plot(T[0][T_nums[t]]['X'],T[0][T_nums[t]]['Z'],'k-')
        axx.plot(T[1][T_nums[t]]['X'],T[1][T_nums[t]]['Z'],'--',color='grey')
        ax1.plot(T[0][T_nums[t]]['Y'],T[0][T_nums[t]]['X'],'k',linewidth=0.5)

    T_nums = [81,91,101,111]
    for t in range(0,len(T_nums)):
        axx = eval('axt2_'+str(t+1))
        axx.plot(T[0][T_nums[t]]['X'],T[0][T_nums[t]]['Z'],'k-')
        axx.plot(T[1][T_nums[t]]['X'],T[1][T_nums[t]]['Z'],'--',color='grey')
        ax1.plot(T[0][T_nums[t]]['Y'],T[0][T_nums[t]]['X'],'k',linewidth=0.5)

    T_nums = [161,171,181,194]
    for t in range(0,len(T_nums)):
        axx = eval('axt3_'+str(t+1))
        axx.plot(T[0][T_nums[t]]['X'],T[0][T_nums[t]]['Z'],'k-')
        axx.plot(T[1][T_nums[t]]['X'],T[1][T_nums[t]]['Z'],'--',color='grey')
        ax1.plot(T[0][T_nums[t]]['Y'],T[0][T_nums[t]]['X'],'k',linewidth=0.5)

    fig.show()
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/figs/changesOverview1.png',dpi=350)



def plotDEMsWithRegionsAndTransects_5km():

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

    h = ax1.contourf(yy1,xx1,np.transpose(dsms[0]),levels=[0,1,2,3,4,5,6,7,8],vmin=0,vmax=8,cmap='gist_earth')
    ax1.set_ylim(-75,125)
    ax1.set_xlim(-1000,4000)
    ax1.invert_yaxis()

    ax2.contourf(yy1,xx1,np.transpose(dsms[1]),levels=[0,1,2,3,4,5,6,7,8],vmin=0,vmax=8,cmap='gist_earth')
    ax2.set_ylim(-75,125)
    ax2.set_xlim(-1000,4000)
    ax2.invert_yaxis()

    ax3.pcolor(yy1,xx1,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')    
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
            axx.set_ylim(0,9)
            ax3.plot(T[0][T_nums[i][t]]['Y'],T[0][T_nums[i][t]]['X'],'k',linewidth=0.5)            


    fig.show()
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/figs/changesOverview_5km.png',dpi=350)






