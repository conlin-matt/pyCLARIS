# Standard library imports #
import os
os.chdir('/Users/frfuser/Documents/pyCLARIS_project')


# 3rd party inputs #
import datetime
import math
from matplotlib import cm
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
from scipy.interpolate import griddata

# Project imports
import pyCLARIS.pyCLARISAnalysis as claris
import pyCLARIS.coastalGeoUtils as utils


# DEMGridSize = '0.5' # '1' or '0.5' #

# files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if 'scarpResults' in i])

# dates = np.unique([int(i[0:8]) for i in files])
# dates_post = np.unique([int(i[9:17]) for i in files])
# numDates = len(dates)

# # Put the mean and std results in a DataFrame #
# results = pd.DataFrame(columns=['Date','Scarp','Method','u','sigma'])
# for d in range(0,numDates):
#     for method in ['manual_pc','_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m','pd','mc','rr','ml']:
#         try:
#             file = [i for i in files if str(dates[d])+'-' in i and method in i][0]
#             f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
#             scarpResults = pickle.load(f) 
#             BT = scarpResults.BT
#             for scarp in range(0,len(BT)):
#                 u = np.nanmean(BT[scarp])
#                 sigma = np.nanstd(BT[scarp])
                
#                 dat = [dates[d],scarp+1,method,u,sigma]
#                 results = results.append({'Date':dates[d],'Scarp':scarp+1,'Method':method,'u':u,'sigma':sigma},ignore_index=True)
#         except:
#             pass
            
            
def plotExampleProfiles_UpAndDown():
    
    method = 'mc'
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if 'scarpResults' in i])
    files_use = [i for i in files if method in i]
    file_up = files_use[0] # Choose your storm #
    file_down = files_use[2]

    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file_up,'rb')
    scarpResults = pickle.load(f)
    T_scarps_up = scarpResults.T_scarps
    BT_up = scarpResults.BT
    toes_up = scarpResults.scarpToes
    
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file_down,'rb')
    scarpResults = pickle.load(f)
    T_scarps_down = scarpResults.T_scarps
    BT_down = scarpResults.BT
    toes_down = scarpResults.scarpToes    
    
    
    fig = plt.figure(figsize=(8.9,9))
    plt.rc('axes', labelsize=24)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=22)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
    plt.rc('legend', fontsize=24)    # legend fontsize
    
    ax_up = plt.axes([0.1,0.63,0.85,0.35])
    ax_down = plt.axes([0.1,0.15,0.85,0.35])
   
    ax_up.set_xlabel('Alongshore dist. (m)')
    ax_up.set_ylabel('Elevation (m)')
    ax_down.set_xlabel('Alongshore dist. (m)')
    ax_down.set_ylabel('Elevation (m)')
    
    ax_up.plot(T_scarps_up[0][0][15]['X'],T_scarps_up[0][0][15]['Z'],'k',linewidth=5,label='Pre-storm')
    ax_up.plot(T_scarps_up[0][1][15]['X'],T_scarps_up[0][1][15]['Z'],'grey',linewidth=5,label='Post-storm')
    ax_up.plot(toes_up[0][0][15][0],toes_up[0][0][15][2],'ko',markersize=12)    
    ax_up.plot(toes_up[0][1][15][0],toes_up[0][1][15][2],'o',color='grey',markersize=12)  
    ax_up.legend(loc='upper right')
    
    ax_down.plot(T_scarps_down[0][0][99]['X'],T_scarps_down[0][0][99]['Z'],'k',linewidth=5,label='Pre-storm')
    ax_down.plot(T_scarps_down[0][1][99]['X'],T_scarps_down[0][1][99]['Z'],'grey',linewidth=5,label='Post-storm')
    ax_down.plot(toes_down[0][0][99][0],toes_down[0][0][99][2],'ko',markersize=12)    
    ax_down.plot(toes_down[0][1][99][0],toes_down[0][1][99][2],'o',color='grey',markersize=12)  
    ax_down.legend(loc='upper right')
    
    for ii in range(0,len(T_scarps_down[0][0])):
        fig,ax = plt.subplots(1)
        ax.plot(T_scarps_down[0][0][ii]['X'],T_scarps_down[0][0][ii]['Z'],'k',linewidth=5,label='Pre-storm')
        ax.plot(T_scarps_down[0][1][ii]['X'],T_scarps_down[0][1][ii]['Z'],'grey',linewidth=5,label='Post-storm')
        ax.plot(toes_down[0][0][ii][0],toes_down[0][0][ii][2],'ko',markersize=12)    
        ax.plot(toes_down[0][1][ii][0],toes_down[0][1][ii][2],'o',color='grey',markersize=12)  


            
def plotExampleProfiles_Longshore():

    
    method = 'mc'
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '5km_T' in i])
    file = files[4]

    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
    T = pickle.load(f)
    mgr = claris.scarpManager(thresh_vertChange=0.5,thresh_longshoreContinuity=100,thresh_minimumElev=1.5,thresh_slope_after=35)
    toes = mgr.idDuneToe(T,method='mc_supervised')
    
    # Estimate berm width #
    widths = []
    heights = []
    for i in range(0,len(T[0])):
        xx = T[0][i]['X']
        zz = T[0][i]['Z']
        dzdx = np.divide(np.diff(zz),np.diff(xx))
        def smooth(x,window_len=7,window='hanning'):
            s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
            #print(len(s))
            if window == 'flat': #moving average
                w=np.ones(window_len,'d')
            else:
                w=eval('np.'+window+'(window_len)')
                            
            y=np.convolve(w/w.sum(),s,mode='valid')
            return y[round(window_len/2-1):-round(window_len/2)]

        dzdx_smooth = smooth(dzdx,window_len=25)
        try:
            i_SignChange = np.where(np.diff(np.sign(dzdx_smooth).astype(int))==2)[0][-1] # Find the most offshore point where the slope changes from negetive to postive (i.e. the landward end of the berm #
        except:
            widths.append(np.nan)
            heights.append(np.nan)
        else:       
            x_land = xx[i_SignChange]
            x_sea = utils.transectElevIntercept(1,xx,zz)

            bermWidth = x_sea-x_land
            bermHeight = max(zz[xx>x_land])-zz[xx==x_land][0]

            widths.append(bermWidth)
            heights.append(bermHeight)
        
##        fig,ax = plt.subplots(2,1,sharex=True)
##        ax[0].plot(xx,zz)
##        ax[1].plot(xx[1:len(xx)],smooth(dzdx,window_len=25))
##        ax[0].plot(xx[i_SignChange],zz[i_SignChange],'ro')
##        ax[0].set_title('width = '+str(round(bermWidth,1))+', height = '+str(round(bermHeight,1)))
##        plt.show()
##        plt.pause(2)
##        plt.close('all')


    numProfs=5
    yLocs = [350,775,1500,2300,3220]
    def makeFig():
        
        fig = plt.figure(figsize=(6.5,4))
        ax0 = plt.axes([0.1,0.6,0.8,0.35])
        
        axs1 = plt.axes([0.1,0.15,0.15,0.3])
        axs2 = plt.axes([0.27,0.15,0.15,0.3],sharey=axs1)
        axs3 = plt.axes([0.44,0.15,0.15,0.3],sharey=axs1)
        axs4 = plt.axes([0.61,0.15,0.15,0.3],sharey=axs1)
        axs5 = plt.axes([0.78,0.15,0.15,0.3],sharey=axs1)

        ax0.spines['left'].set_color('b')
        ax0.xaxis.label.set_color('b')
        ax0.tick_params(axis='y', colors='b')
        ax0.set_ylabel('Berm height (m)',color='b')
        ax0.set_xlabel('Alongshore (m)',color='k')
        
        ax0.plot(np.arange(0,4000,5),heights,'bs',markersize=1.5)
        yl=ax0.get_ylim()
        r1 = Rectangle((2180,yl[0]),2840-2180,np.diff(yl),facecolor='k',alpha=0.4)
        ax0.add_artist(r1)
        ax0.set_ylim(yl)
        r2 = Rectangle((0,yl[0]),450,np.diff(yl),facecolor='k',alpha=0.4)
        ax0.add_artist(r2)
        ax0.set_ylim(yl)
        for y in yLocs:
            ax0.plot((y,y),yl,'k--')
            ax0.set_ylim(yl)

        ax1 = ax0.twinx()   
        ax1.plot(np.arange(0,4000,5),widths,'ks',markersize=1.5)
        ax1.set_ylabel('Berm width (m)')

        axs = [axs1,axs2,axs3,axs4,axs5]
        for i in range(0,len(axs)):
            tt = T[0]
            yy = [ti['Y'][0] for ti in tt]
            iUse = int(np.where(np.array(yy)==yLocs[i])[0])

            h1 = axs[i].plot(T[0][iUse]['X'],T[0][iUse]['Z'],'k',linewidth=1)
            h2 = axs[i].plot(T[1][iUse]['X'],T[1][iUse]['Z'],'grey',linewidth=1)
            xl = axs[i].get_xlim()
            axs[i].set_xticks([round((xl[0]+10)/10)*10,round((xl[1]-10)/10)*10])
            axs[i].set_xlim(xl)
            if i==0:
                axs[i].set_xlabel('xShore (m)')
                axs[i].set_ylabel('Elev (m)')
            else:
                axs[i].set_yticklabels([])
        fig.legend([h1[0],h2[0]],['Pre','Post'],ncol=2,loc='lower right')

        plt.show()



        def makeManyProfsFig():

            import matplotlib
            cmap = matplotlib.cm.get_cmap('jet')
             
            ys1 = np.arange(1200,1500,25)
            ys2 = np.arange(1700,2100,30)
            ys3 = np.arange(2200,2800,40)

            fig,ax = plt.subplots(1,figsize=(6.5,4))
            hh = []
            for s in [ys1,ys2,ys3]:
                if s[0]==ys1[0]:
                    cols = cmap(np.linspace(0,0.2,len(s)))
                elif s[0]==ys2[0]:
                    cols = cmap(np.linspace(0.4,0.6,len(s)))
                elif s[0]==ys3[0]:
                    cols = cmap(np.linspace(0.8,1,len(s)))
        

                for i in range(0,len(s)):
                    loc=s[i]
                    tt = T[0]
                    yy = [ti['Y'][0] for ti in tt]
                    iUse = int(np.where(np.array(yy)==loc)[0])

                    xx = T[0][iUse]['X']
                    zz = T[0][iUse]['Z']
                    iRef = np.where(abs(zz-3)==min(abs(zz-3)))[0][0]
                    xxRel = xx-xx[iRef]
                    if i==0:
                        h= ax.plot(xxRel,zz,color=cols[i])
##                        h = ax.plot((xx-min(xx))/(max(xx)-min(xx)),zz,color=cols[i])
##                        h = ax.plot(xx,zz,color=cols[i])
                        hh.append(h[0])
                    else:
                        ax.plot(xxRel,zz,color=cols[i])
##                        ax.plot((xx-min(xx))/(max(xx)-min(xx)),zz,color=cols[i])
##                        ax.plot(xx,zz,color=cols[i])
            fig.legend(hh,['Y=1200-1500\nSmall berm, no upper beach change','Y=1700-2100\nBig berm, no upper beach change','Y=2200-2800\nBig berm, upper beach change'],loc='upper right')
            ax.set_xlabel('Cross-shore (m)')
            ax.set_ylabel('Elevation (m)')
            plt.show()


        def makeManyProfsFig_scatters():
             
#            ys1 = np.arange(1000,2000,50)
##            ys2 = np.arange(1700,2100,30)
            ys3 = np.arange(1000,2000,5)

            fig = plt.figure(figsize=(3.25,4))
            ax = plt.axes([0.12,0.11,0.8,0.6])
##            cbax1 = plt.axes([0.12,0.9,0.4,0.02])
##            cbax2 = plt.axes([0.12,0.85,0.4,0.02])
            cbax3 = plt.axes([0.12,0.8,0.4,0.02])

            ax.text(0.01,0.01,'a',fontsize=8,fontweight='bold',transform=ax.transAxes)
##            cbax1.text(1.02,0.5,'Y=1000-2000',ha='left',va='center',fontsize=8,transform=cbax1.transAxes,fontweight='bold')
##            cbax1.text(0.5,1.1,'$\Delta Z (m)$',ha='center',va='bottom',fontsize=8,transform=cbax1.transAxes)
##            cbax2.text(1.02,0.5,'Y=1700-2100',ha='left',va='center',fontsize=8,transform=cbax2.transAxes,fontweight='bold')
            cbax3.text(1.02,0.5,'Y=2200-3200',ha='left',va='center',fontsize=8,transform=cbax3.transAxes,fontweight='bold')
            
            hh = []
            it=-1
            for s,c in zip([ys3],[cbax3]):
                it+=1

##                if s[0]==ys1[0]:
##                    cols = 'seismic_r'
##                elif s[0]==ys2[0]:
##                    cols = 'PuOr'
##                elif s[0]==ys3[0]:
##                cols = 'PuOr'

                xi = np.arange(-100,0,0.25)
                ziPre_all = np.empty([0,len(xi)])
                dz_all = np.empty([0,len(xi)])
                for i in range(0,len(s)):
                    loc=s[i]
                    tt = T[0]
                    yy = [ti['Y'][0] for ti in tt]
                    iUse = int(np.where(np.array(yy)==loc)[0])

                    xx = T[0][iUse]['X']
                    zzPre = T[0][iUse]['Z']
                    xxPost = T[1][iUse]['X']
                    zzPost = T[1][iUse]['Z']
                    zzPosti = np.interp(xx,xxPost,zzPost)
                    iRef = np.where(abs(zzPre-5)==min(abs(zzPre-5)))[0][0]
                    xxRel = xx-xx[iRef]

                    
                    ziPre = np.interp(xi,xx-max(xx),zzPre)
                    ziPost = np.interp(xi,xxPost-max(xxPost),zzPost)
                    dz = np.array(ziPost-ziPre)

                    ziPre_all = np.vstack([ziPre_all,np.array(ziPre).reshape(1,-1)])
                    dz_all = np.vstack([dz_all,np.array(dz).reshape(1,-1)])

                zi_mean = np.nanmean(ziPre_all,axis=0)
                zi_std = np.nanstd(ziPre_all,axis=0)
                dz_mean = np.nanmean(dz_all,axis=0)


                    

                    h = ax.scatter(xxRel,zzPre,2,zzPosti-zzPre,cmap=cols,vmin=-1,vmax=1)
                    if c==cbax3:
                        plt.colorbar(h,c,orientation='horizontal',ticks=[-1,0,1])
                    else:
                        cb = plt.colorbar(h,c,orientation='horizontal',ticks=None)
                        cb.set_ticks([])

                ax.set_xlabel('Relative cross-shore (m)')
                ax.set_ylabel('Elevation (m)')
                ax.set_ylim(0,4)
##                ax.set_xlim(0,45)
                ax.set_yticks([0,1,2,3,4])

##            plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/TeddyBerms.png',dpi=400)
            plt.show()
            
                    
        def makeManyProfsFig_meanProf():
             
            from matplotlib.colors import ListedColormap

            def createRYBcmap():
                m = 256
                m1=m*0.5
                r = np.arange(0,m1)/(m1-1)
                g = np.arange(0,m1)/(m1-1)
                b = np.arange(0,m1)/(m1-1)
                r = np.vstack([r.reshape(-1,1),np.ones([len(r),1])])
                g = np.vstack([g.reshape(-1,1),np.flipud(g.reshape(-1,1))])
                b = np.vstack([b.reshape(-1,1),np.ones([len(b),1])])
                b = np.flipud(b)
                b = np.flipud(np.linspace(0,1,len(r))).reshape(-1,1)
                
                c = np.hstack([r,g,b,np.ones([len(r),1])])
                return np.flipud(c)
                
            c = createRYBcmap()
            ryb = ListedColormap(c)

            ys1 = np.arange(1000,2005,5)
            ys3 = np.arange(2200,3205,5)

            fig,ax = plt.subplots(2,1,figsize=(3,5),sharex=True,sharey=True)
            cbax = plt.axes([0.2,0.94,0.6,0.02])

            ax[0].text(0.01,0.01,'a (Y=1000-2000 m)',fontsize=8,fontweight='bold',transform=ax[0].transAxes)
            ax[1].text(0.01,0.01,'b (Y=2200-3200 m)',fontsize=8,fontweight='bold',transform=ax[1].transAxes)
            cbax.text(0.5,1.1,'$\Delta Z (m)$',ha='center',va='bottom',fontsize=8,transform=cbax.transAxes)
            
            hh = []
            it=-1
            for s,axx in zip([ys1,ys3],ax):
                it+=1

                xi = np.arange(-100,0,0.25)
                ziPre_all = np.empty([0,len(xi)])
                dz_all = np.empty([0,len(xi)])
                slope_all = []
                b_all = []
                for i in range(0,len(s)):
                    loc=s[i]
                    tt = T[0]
                    yy = [ti['Y'][0] for ti in tt]
                    iUse = int(np.where(np.array(yy)==loc)[0])

                    xx = T[0][iUse]['X']
                    zzPre = T[0][iUse]['Z']
                    xxPost = T[1][iUse]['X']
                    zzPost = T[1][iUse]['Z']
                    zzPosti = np.interp(xx,xxPost,zzPost)
                    iRef = np.where(abs(zzPre-5)==min(abs(zzPre-5)))[0][0]
                    xxRel = xx-xx[iRef]

                    
                    ziPre = np.interp(xi,xx-max(xx),zzPre)
                    ziPost = np.interp(xi,xxPost-max(xxPost),zzPost)
                    dz = np.array(ziPost-ziPre)

                    ziPre_all = np.vstack([ziPre_all,np.array(ziPre).reshape(1,-1)])
                    dz_all = np.vstack([dz_all,np.array(dz).reshape(1,-1)])

                    xxSlope = xi[np.logical_and(xi>=-15,xi<=-10)]
                    zzSlope = ziPre[np.logical_and(xi>=-15,xi<=-10)]
                    bb = utils.linearRegression(xxSlope,zzSlope)[0]
                    b = bb[0]
                    slope = bb[1]
                    slope_all.append(slope)
                    b_all.append(b)
                    

                zi_mean = np.nanmean(ziPre_all,axis=0)
                zi_std = np.nanstd(ziPre_all,axis=0)
                dz_mean = np.nanmean(dz_all,axis=0)
                slope_mean = np.mean(slope_all)
                slope_std = np.std(slope_all)
                b_mean = np.mean(b_all)
                print('mean='+str(round(slope_mean,5))+', std='+str(round(slope_std,5)))

              
                h = axx.scatter(xi,zi_mean,5,dz_mean,cmap=ryb,vmin=-1,vmax=1,zorder=10)
                axx.plot(xi,zi_mean+zi_std,'k',linewidth=1,zorder=11)
                axx.plot(xi,zi_mean-zi_std,'k',linewidth=1,zorder=12)
                axx.plot((-15,-15),(0,5),'k--',linewidth=1,zorder=1)
                axx.plot((-10,-10),(0,5),'k--',linewidth=1,zorder=2)
                axx.plot(np.arange(-30,-5,0.25),(np.arange(-30,-5,0.25)*slope_mean)+b_mean,'k',linewidth=1,zorder=8)
                if axx == ax[0]:
                    axx.text(-38,1.4,'m = '+str(round(slope_mean,3))+'$\pm$ '+ str(round(slope_std,3)),color='k')
                else:
                    axx.text(-38,1.1,'m = '+str(round(slope_mean,3))+'$\pm$ '+ str(round(slope_std,3)),color='k')                    
                axx.set_xlim(-40,2)
                axx.set_ylim(0,5)
                axx.set_xticklabels([0,10,20,30,40])
                if axx == ax[1]:
                    axx.set_xlabel('Cross-shore (m)')
                axx.set_ylabel('Elev. (m)')

                plt.colorbar(h,cbax,orientation='horizontal',ticks=[-1,0,1])
                
                
            plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/TeddyBerms.png',dpi=400)
            plt.show()


                    

     
    
def plotProfilesThroughTime():
    
    import matplotlib
    
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '5km_T' in i])
    files = [i for i in files if '2011' not in i and '20170918' not in i and '202003' not in i] # Remove Irene and Jose #
    
    y_locs = [200,1900] # Alongshore location of profile to plot #
    
    cmap = matplotlib.cm.get_cmap('nipy_spectral')
    cols = cmap([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])

    
    fig = plt.figure(figsize=(6.5,4))
    ax1 = plt.axes([0.1,0.1,0.4,0.5])
    ax2 = plt.axes([0.55,0.1,0.4,0.5],sharey=ax1)

    hh = []
    iterr=-1
    for file in files:
        iterr+=1
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        T = pickle.load(f)
        
        if file==files[0]: # Find the transect to use #
            tt = T[0]
            yy = [ti['Y'][0] for ti in tt]
            iUse1 = int(np.where(np.array(yy)==y_locs[0])[0])
            iUse2 = int(np.where(np.array(yy)==y_locs[1])[0])

        
        hh.append(ax1.plot(T[0][iUse1]['X'],T[0][iUse1]['Z'],color=cols[iterr]))
        hh.append(ax1.plot(T[1][iUse1]['X'],T[1][iUse1]['Z'],'--',color=cols[iterr]))
        ax2.plot(T[0][iUse2]['X'],T[0][iUse2]['Z'],color=cols[iterr])
        ax2.plot(T[1][iUse2]['X'],T[1][iUse2]['Z'],'--',color=cols[iterr])

    hh = [i[0] for i in hh]
    hh = hh[0:len(hh):2]+hh[1:len(hh):2]
    names = ['pre 2013NE','post 2013 NE',
                  'pre Maria','post Maria',
                  'pre Riley','post Riley',
                  'pre Dorian','post Dorian',
                  'pre Humberto','post Humberto',
                  'pre 2019NE1','post 2019NE1',
                  'pre 2019NE2','post 2019NE2',
                  'pre Teddy','post Teddy',
                  'pre 2021NE','post 2021NE']
    names = names[0:len(names):2]+names[1:len(names):2]  
    fig.legend(hh,names,ncol=2,loc='upper center',columnspacing=0.5,handletextpad=0.3)
    ax1.set_xlabel('Cross-shore (m)')
    ax2.set_xlabel('Cross-shore (m)')
    ax1.set_ylabel('Elev (m, NAVD88)')
    ax2.set_ylabel('')

    ax1.text(0.02, 0.1, 'a (Y = '+str(y_locs[0])+')',ha='left', va='top',transform=ax1.transAxes,fontweight='bold')
    ax2.text(0.02, 0.1, 'b (Y = '+str(y_locs[1])+')',ha='left', va='top',transform=ax2.transAxes,fontweight='bold')

    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/ProfilesThroughTime.png',dpi=400)
    plt.show()
        
        
def plotDoDs():
    
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '5km_dsms' in i])
    files = [i for i in files if '2011' not in i and '20170918' not in i] # Remove Irene and Jose #   

    fig,ax = plt.subplots(9,1,figsize=(6.5,9))    
    plt.subplots_adjust(top=0.94,bottom=0.1,left=0.1,right=0.95,hspace=0.1)
    cbax = plt.axes([0.325,0.965,0.4,0.007])
    names = ['i (2013NE)','h (Maria)','g (Riley)','f (Dorian)','e (Humberto)','d (2019NE1)','c (2019NE2)','b (Teddy)','a (2021NE)']
    
    for i in range(0,len(files)):        
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+files[i],'rb')
        dsms = pickle.load(f)  
        xx = np.arange(-200,150,0.5)
        yy = np.arange(0,4000,0.5)
        h = ax[-i-1].pcolor(yy,xx,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')
        ax[-i-1].contour(yy,xx,np.transpose(dsms[0]),[0.36],colors=['k'])
        ax[-i-1].contour(yy,xx,np.transpose(dsms[0]),[4],colors=['k'])
        ax[-i-1].set_ylim(-80,120)
        ax[-i-1].set_yticks([60,110])
        ax[-i-1].invert_yaxis()
        ax[-i-1].set_xticks([0,1000,2000,3000,4000])
        ax[-i-1].text(50,-20,names[i],fontsize=8,fontweight='bold',ha='left',va='bottom')
        if i!=0:
            ax[-i-1].set_xticklabels([])
            ax[-i-1].set_yticklabels([])
        else:
            ax[-i-1].set_xlabel('Alongshore (m)',fontsize=10)
            ax[-i-1].set_ylabel('Cross-shore (m)',fontsize=10)
            ard = Arrow(3800,60,0,25,width=50,color='k')
            aru = Arrow(3800,60,0,-25,width=50,color='k')
            ax[-i-1].add_patch(ard)
            ax[-i-1].add_artist(aru)            
            ax[-i-1].text(3800,85,'offshore',ha='center',va='top',fontsize=8)
            ax[-i-1].text(3800,35,'onshore',ha='center',va='bottom',fontsize=8)
            arl = Arrow(3200,95,-100,0,width=25,color='k')
            arr = Arrow(3200,95,100,0,width=25,color='k')
            ax[-i-1].add_patch(arl)
            ax[-i-1].add_artist(arr)  
            ax[-i-1].text(3075,100,'south',va='center',ha='right',fontsize=8)
            ax[-i-1].text(3325,100,'north',va='center',ha='left',fontsize=8)
            
        cb = plt.colorbar(h,cbax,orientation='horizontal')
        cb.ax.set_title('$\Delta Z $(m)',fontsize=10)
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/DoDs.png',dpi=400)
    plt.show()
            
        

def plotBTValsAll_Hovmoller():
    '''
    Plots every computed BT value as a function of space and time #

    Returns
    -------
    None.

    '''
    
    import datetime
    from matplotlib.colors import ListedColormap
    import matplotlib.dates as mdates
    myFmt = mdates.DateFormatter('%Y')

    def createRYBcmap():
        m = 256
        m1=m*0.5
        r = np.arange(0,m1)/(m1-1)
        g = np.arange(0,m1)/(m1-1)
        b = np.arange(0,m1)/(m1-1)
        r = np.vstack([r.reshape(-1,1),np.ones([len(r),1])])
        g = np.vstack([g.reshape(-1,1),np.flipud(g.reshape(-1,1))])
        b = np.vstack([b.reshape(-1,1),np.ones([len(b),1])])
        b = np.flipud(b)
        b = np.flipud(np.linspace(0,1,len(r))).reshape(-1,1)
        
        c = np.hstack([r,g,b,np.ones([len(r),1])])
        return np.flipud(c)
        
    c = createRYBcmap()
    ryb = ListedColormap(c)
    
    yy = np.array([datetime.datetime(2011,1,1)+datetime.timedelta(days=i) for i in range(0,(365*11)+20,100)])
    xx = np.arange(-1000,4100,100)
    
    # Plot the storms with labels #
    fig = plt.figure(figsize=(4,6))
    plt.rcParams.update({'font.size': 8})
    ax1 = plt.axes([0.1,0.08,0.72,0.25])
    ax1.set_xlim(0,4000)
    ax1.set_ylim(datetime.datetime(2013,1,1),datetime.datetime(2014,12,31))
    ax1.spines['top'].set_visible(False)
    ax1.set_yticks([datetime.datetime(2013,1,1),datetime.datetime(2014,1,1)])
    ax1.yaxis.set_major_formatter(myFmt)
    ax2 = plt.axes([0.1,0.33,0.72,0.62])
    ax2.set_xlim(0,4000)
    ax2.set_ylim(datetime.datetime(2017,1,1),datetime.datetime(2021,12,13)) 
    ax2.spines['bottom'].set_visible(False)
    ax2.set_xticks([])
    ax2.set_yticks([datetime.datetime(2017,1,1),datetime.datetime(2018,1,1),datetime.datetime(2019,1,1),datetime.datetime(2020,1,1),datetime.datetime(2021,1,1)])
    ax2.yaxis.set_major_formatter(myFmt)
    
    # # AGU poster format #
    # fig = plt.figure(figsize=(6.2,7.5))
    # plt.rc('axes', labelsize=22)    # fontsize of the x and y labels
    # plt.rc('xtick', labelsize=22)    # fontsize of the tick labels
    # plt.rc('ytick', labelsize=22)    # fontsize of the tick labels
    # plt.rc('legend', fontsize=16)     
    # ax1 = plt.axes([0.18,0.12,0.6,0.25])
    # ax1.set_xlim(-1000,4000)
    # ax1.set_ylim(datetime.datetime(2013,1,1),datetime.datetime(2014,12,31))
    # ax1.spines['top'].set_visible(False)
    # ax1.set_yticks([datetime.datetime(2013,1,1),datetime.datetime(2014,1,1)])
    # ax1.yaxis.set_major_formatter(myFmt)
    # ax2 = plt.axes([0.18,0.37,0.6,0.62])
    # ax2.set_xlim(-1000,4000)
    # ax2.set_ylim(datetime.datetime(2017,1,1),datetime.datetime(2021,12,13)) 
    # ax2.spines['bottom'].set_visible(False)
    # ax2.set_xticks([])
    # ax2.set_yticks([datetime.datetime(2017,1,1),datetime.datetime(2018,1,1),datetime.datetime(2019,1,1),datetime.datetime(2020,1,1),datetime.datetime(2021,1,1)])
    # ax2.yaxis.set_major_formatter(myFmt)
    
    d = .015  # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (.7-d, .7+d), **kwargs)        # top-left diagonal
    ax1.plot((1 - d, 1 + d), (.7-d, .7+d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (-.09-d, -.11+d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (-.09-d, -.11+d), **kwargs)  # bottom-right diagonal

    
    for ax in [ax1,ax2]:
        dates_storms = [datetime.datetime(2013,3,13),datetime.datetime(2017,9,25),
                        datetime.datetime(2018,3,6),datetime.datetime(2019,9,7),
                        datetime.datetime(2019,9,17),datetime.datetime(2019,10,12),
                        datetime.datetime(2019,11,16),datetime.datetime(2020,3,7,12),
                        datetime.datetime(2020,9,17,12),datetime.datetime(2021,3,20,12)]
        for i in range(0,len(dates_storms)):
            ax.plot((min(xx),max(xx)),(dates_storms[i],dates_storms[i]),'k--',linewidth=0.5)
            
    ax1.text(4050,dates_storms[0],'NorEaster',va='center',ha='left',fontweight='bold')
    ax2.text(4050,dates_storms[1],'Maria',va='center')
    ax2.text(4050,dates_storms[2],'Riley',va='center',ha='left',fontweight='bold')
    ax2.text(4050,dates_storms[3]-datetime.timedelta(days=50),'Dorian',va='center',ha='left',fontweight='bold')
    ax2.text(4050,dates_storms[4]-datetime.timedelta(days=20),'Humberto',va='center',ha='left')
    ax2.text(4050,dates_storms[5]+datetime.timedelta(days=5),'NorEaster',va='center',ha='left',fontweight='bold')
    ax2.text(4050,dates_storms[6]+datetime.timedelta(days=15),'NorEaster',va='center',ha='left',fontweight='bold')
    ax2.text(4050,dates_storms[7],'NorEaster',va='center',ha='left')   
    ax2.text(4050,dates_storms[8],'Teddy',va='center',ha='left',fontweight='bold')
    ax2.text(4050,dates_storms[9],'NorEaster',va='center',ha='left',fontweight='bold')
    ax1.set_xlabel('FRF Y (m)')
    
    # # AGU poster format #
    # ax1.text(4001,dates_storms[0],'NE',va='center',ha='left',fontweight='bold',fontsize=18)
    # ax2.text(4001,dates_storms[1],'Maria',va='center',fontsize=18)
    # ax2.text(4001,dates_storms[2],'Riley',va='center',ha='left',fontweight='bold',fontsize=18)
    # ax2.text(4001,dates_storms[3]-datetime.timedelta(days=120),'Dorian',va='center',ha='left',fontweight='bold',fontsize=18)
    # ax2.text(4001,dates_storms[4]-datetime.timedelta(days=60),'Humberto',va='center',ha='left',fontsize=18)
    # ax2.text(4001,dates_storms[5]-datetime.timedelta(days=15),'NE',va='center',ha='left',fontweight='bold',fontsize=18)
    # ax2.text(4001,dates_storms[6]+datetime.timedelta(days=25),'NE',va='center',ha='left',fontweight='bold',fontsize=18)
    # ax2.text(4001,dates_storms[7],'NE',va='center',ha='left',fontsize=18)   
    # ax2.text(4001,dates_storms[8],'Teddy',va='center',ha='left',fontweight='bold',fontsize=18)
    # ax2.text(4001,dates_storms[9],'NE',va='center',ha='left',fontweight='bold',fontsize=18)
    # ax1.set_xlabel('alongshore (m)')
    
    # Plot the data #
    method = 'mc'
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if 'scarpResults' in i])
    files_use = [i for i in files if method in i]
    c=-1
    for file in files_use:
        c+=1
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        scarpResults = pickle.load(f)
        T_scarps = scarpResults.T_scarps
        BT = scarpResults.BT
        forwardFlag = scarpResults.forwardFlag
        
        date_pre = datetime.datetime(int(file[0:4]),int(file[4:6]),int(file[6:8]))
        date_post = datetime.datetime(int(file[9:13]),int(file[13:15]),int(file[15:17]))
        date_plot = date_pre+((date_post-date_pre)/2)       
        # Move 2019 data to help visualization #
        if c==2:
            date_plot = date_plot-datetime.timedelta(days=40)
        elif c==4:
            date_plot = date_plot+datetime.timedelta(days=40)
        
        
        for region in range(0,len(BT)):
            locs = [T_scarps[region][0][i]['Y'][0] for i in range(0,len(T_scarps[region][0]))]
            times = np.tile(date_plot,len(locs))
            BT_vals = BT[region]
            BT_vals[forwardFlag[region]] = np.nan # Do not plot those where dToe advanced #
            
            for ax in [ax1,ax2]:
                a = ax.scatter(locs,times,200,BT_vals,marker='|',cmap=ryb,vmin=-.2,vmax=.2,alpha=1)
                # a = ax.scatter(locs,times,600,BT_vals,marker='|',cmap=ryb,vmin=-.2,vmax=.2,alpha=1)


    cbax = plt.axes([0.33,0.3,0.3,0.01])
    # cbax = plt.axes([0.35,0.42,0.3,0.01])    
    cbax.set_xticks([])
    cbax.set_yticks([])
    plt.colorbar(a,cbax,orientation='horizontal',label='$B_T$')    
    
    

def compareDtoeMethods():

    method1 = 'ml'
    method2 = 'mc'
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if 'scarpResults' in i])
    files_use1 = [i for i in files if method1 in i]
    files_use2 = [i for i in files if method2 in i]
        
    bigArray1 = np.empty([0,7])
    bigArray2 = np.empty([0,7])
    for i in range(0,max(len(files_use1),len(files_use1))):
        file1 = files_use1[i]
        dateTest = file1[0:8]
        file2 = [i for i in files_use2 if dateTest in i]
        
        if len(file2)>0:
            file2 = file2[0]
            f1 = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file1,'rb') 
            scarpResults1 = pickle.load(f1)
            f2 = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file2,'rb') 
            scarpResults2 = pickle.load(f2)
            
            # Create an array of date,y-value,pre/post storm toe x and z, and BT for both methods #
            c=-1
            for scarpResults in [scarpResults1,scarpResults2]:
                c+=1
                for scarp in range(0,len(scarpResults.T_scarps)):
                    y_vals = [scarpResults.T_scarps[scarp][0][ii]['Y'][0] for ii in range(0,len(scarpResults.T_scarps[scarp][0]))]
                    toes_pre_x = [scarpResults.scarpToes[scarp][0][ii][0] for ii in range(0,len(scarpResults.T_scarps[scarp][0]))]
                    toes_pre_z = [scarpResults.scarpToes[scarp][0][ii][2] for ii in range(0,len(scarpResults.T_scarps[scarp][0]))]
                    toes_post_x = [scarpResults.scarpToes[scarp][1][ii][0] for ii in range(0,len(scarpResults.T_scarps[scarp][0]))]
                    toes_post_z = [scarpResults.scarpToes[scarp][1][ii][2] for ii in range(0,len(scarpResults.T_scarps[scarp][0]))]   
                    BT_vals = scarpResults.BT[scarp] 
                    date_vals = np.tile(int(file1[0:8]),(len(BT_vals)))
                    ar = np.transpose(np.vstack([date_vals.reshape(1,-1),
                                    np.array(y_vals).reshape(1,-1),
                                    np.array(toes_pre_x).reshape(1,-1),
                                    np.array(toes_pre_z).reshape(1,-1),
                                    np.array(toes_post_x).reshape(1,-1),
                                    np.array(toes_post_z).reshape(1,-1),
                                    np.array(BT_vals).reshape(1,-1)
                                    ]))
                    if c==1:
                        forwardFlag = scarpResults.forwardFlag[scarp]
                        ar_new = np.delete(ar,forwardFlag,axis=0)
                        ar = ar_new
                    if c==0:
                        bigArray1 = np.vstack([bigArray1,ar])
                    elif c==1:
                        bigArray2 = np.vstack([bigArray2,ar])
    
    # Find where transects have been extracted in same place for same storm for both methods and extract BT values #
    compArray_BT = np.empty([0,3])
    for i in range(0,len(bigArray1)):
        test = any(np.logical_and(bigArray2[:,0]==bigArray1[i,0],bigArray2[:,1]==bigArray1[i,1]))
        if test:
            bigArray2_rowToKeep = bigArray2[np.logical_and(bigArray2[:,0]==bigArray1[i,0],bigArray2[:,1]==bigArray1[i,1]),:]
            compArray_BT = np.vstack([compArray_BT,np.hstack([bigArray1[i,0],bigArray1[i,6],bigArray2_rowToKeep[0,6]])])

    # Plot a 1:1 comparison of the methods #
    fig,ax = plt.subplots(1)
    dates = np.unique(compArray_BT[:,0])
    colors = ['r','g','b','k','c','m','y',[0.5,0.8,0.3]]
    h_all = []
    for d in range(0,len(dates)):
        i = np.where(compArray_BT[:,0]==dates[d])
        h = ax.plot(compArray_BT[i,1],compArray_BT[i,2],'.',markersize=2,color=colors[d])
        h_all.append(h[0])
    ax.axis('equal')
    ax.plot(ax.get_xlim(),ax.get_xlim(),'k--')
    ax.plot((0,0),ax.get_ylim(),'k',linewidth=1)
    ax.plot(ax.get_xlim(),(0,0),'k',linewidth=1)
    ax.set_xlabel(method1)
    ax.set_ylabel(method2)
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.legend(h_all,np.unique(compArray_BT[:,0]).astype(int),markerscale=2)
    ax.legend(h_all,['2013NE','Riley','Dorian','2019NE1','2019NE2','Teddy','2021NE'],markerscale=2)
    propGoodQuad = len(np.where(np.sign(compArray_BT[:,1])==np.sign(compArray_BT[:,2]))[0])/len(compArray_BT)
    propManualNegAutoPos = len(np.where(np.logical_and(np.sign(compArray_BT[:,1])==-1,np.sign(compArray_BT[:,2])==1))[0]) / len(np.where(np.sign(compArray_BT[:,1])!=np.sign(compArray_BT[:,2]))[0])
    
    # # AGU 2021 poster figure #
    # fig = plt.figure(figsize=(6.4,4))
    # plt.rc('axes', labelsize=22)    # fontsize of the x and y labels
    # plt.rc('xtick', labelsize=22)    # fontsize of the tick labels
    # plt.rc('ytick', labelsize=22)    # fontsize of the tick labels
    # plt.rc('legend', fontsize=16) 
    # ax = plt.axes([0.2,0.2,0.75,0.75])
    # dates = np.unique(compArray_BT[:,0])
    # colors = ['r','g','b','k','c','m','y',[0.5,0.8,0.3]]
    # h_all = []
    # for d in range(0,len(dates)):
    #     i = np.where(compArray_BT[:,0]==dates[d])
    #     h = ax.plot(compArray_BT[i,1],compArray_BT[i,2],'.',markersize=4,color=colors[d])
    #     h_all.append(h[0])
    # ax.axis('equal')
    # ax.plot(ax.get_xlim(),ax.get_xlim(),'k--')
    # ax.plot((0,0),ax.get_ylim(),'k',linewidth=1)
    # ax.plot(ax.get_xlim(),(0,0),'k',linewidth=1)
    # ax.set_xlabel('$B_T$ manual')
    # ax.set_ylabel('$B_T$ auto')
    # ax.set_xlim(-1.2,1)
    # ax.set_ylim(-.9,.8)
    # ax.legend(h_all,['2013NE','Riley','Dorian','2019NE1','2019NE2','Teddy','2021NE'],markerscale=2,loc='lower right',handletextpad=0.001,borderpad=0.1) 

    
def examineAndPlot_VarsVsBT_AllWays():   
        
    method = 'mc'
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if 'scarpResults' in i])
    files_use = [i for i in files if method in i]
    results_byStorm = pd.DataFrame(columns=['Name','results'])
    for file in files_use:
        
        results = pd.DataFrame(columns=['BT','Bf','$Z_{Toe}$ (m)','Post-storm toe elev (m)','$Bf_{post}$','dx (m)','dx (m, contour)',
                                'dz (m)','$\B_{Dune}$','Bf/Bd ',
                                'F (m)','W (m)','t_i (hr)','$max(H_s)$ (m)','\sumE (MJh/m^2)','TWL_ref (m)'
                                '$\Delta V_B$ ($m^3/m$)','$\Delta V_D$ ($m^3/m$)','Dune Height (m)'])
        bdate = int(file[0:8])
        edate = int(file[9:17])
        
        if '2013'  not in str(bdate):
            ws = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/frf_17m_waves.xlsx'
            hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=17.8)  
        else:  # Data from 17 m buoy missing for first few months of 2013. Use data from 26 m buoy and transform to 17 m with LWT #
            ws='/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/frf_26m_waves.xlsx'
            hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=26)  
            waves = hydro.waves                      
            for i in range(0,len(waves)):
                dat = waves.iloc[i]
                d = 26
                H = dat['wvht (m)']
                T = dat['Tp (s)']
                theta = dat['MWD (degT)']
                
                if 0<theta<180:               
                    d_vec = np.arange(-d,-16,1) # Vector of water depths #
                    # Wave parameters at buoy #  
                    k0 = utils.newtRaph(T,d)
                    C0 = ((2*np.pi)/k0)/T
                    n0 = 0.5*(1+((2*k0*d)/np.sinh(2*k0*d)))
                    alpha0 = 90-theta
                    Hh = H
                    # Transform #
                    for h in d_vec[1:len(d_vec)]:
                        k = utils.newtRaph(T,-h)
                        C = ((2*np.pi)/k)/T
                        n = 0.5*(1+((2*k*-h)/np.sinh(2*k*-h))) 
                        alpha = np.degrees(math.asin((C/C0)*np.sin(np.radians(alpha0))))
                        Hs = math.sqrt((n0*C0)/(n*C))*math.sqrt(math.cos(math.radians(alpha0))/math.cos(math.radians(alpha)))*Hh
                        
                        k0 = k
                        C0 = C
                        n0 = n
                        alpha0 = alpha
                        Hh = Hs
                        
                    waves.replace({'wvht (m)':waves.iloc[i]['wvht (m)'],'MWD (degT)':waves.iloc[i]['MWD (degT)']},
                                          {'wvht (m)':Hh,'MWD (degT)':360-(270+alpha)},inplace=True)
                else:
                    waves.replace({'wvht (m)':waves.iloc[i]['wvht (m)'],'MWD (degT)':waves.iloc[i]['MWD (degT)']},
                                          {'wvht (m)':np.nan,'MWD (degT)':np.nan},inplace=True)     
                    
            hydro.waves = waves
        

        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        scarpResults=pickle.load(f)
        Bf = scarpResults.Bf
        BT = scarpResults.BT
        scarpToes = scarpResults.scarpToes
        T_scarps = scarpResults.T_scarps
        forwardFlag = scarpResults.forwardFlag
        
        for scarp in range(0,len(T_scarps)):
            for ii in range(0,len(T_scarps[scarp][0])):
                
                if ii not in forwardFlag[scarp]:
                
                    BT_here = BT[scarp][ii]
                    Bf_here = Bf[scarp][ii]
                    
                    
                    # Calc toe elev pre and post #
                    toe_pre = scarpToes[scarp][0][ii][2]
                    toe_post = scarpToes[scarp][1][ii][2]
                    
                    # Calc post-storm Bf #
                    x_beta = T_scarps[scarp][1][ii]['X'][T_scarps[scarp][1][ii]['X']>=scarpToes[scarp][1][ii][0]]
                    z_beta = T_scarps[scarp][1][ii]['Z'][T_scarps[scarp][1][ii]['X']>=scarpToes[scarp][1][ii][0]]
                    try:
                        Bf_post = abs((z_beta[-1]-z_beta[0])/(x_beta[-1]-x_beta[0]))
                    except:
                        Bf_post = np.nan
                    
                    # Calc max wave height #
                    maxHs = float(max(hydro.waves['wvht (m)']))
                    
                    # Calc cumulative wave energy between surveys #
                    t_sec=[]
                    for t in range(0,len(hydro.waves)):
                        t_sec.append((datetime.datetime(int(hydro.waves['yr'][t]),int(hydro.waves['mo'][t]),int(hydro.waves['day'][t]),int(hydro.waves['hr'][t]),int(hydro.waves['mm'][t]))-
                                  datetime.datetime(int(hydro.waves['yr'][0]),int(hydro.waves['mo'][0]),int(hydro.waves['day'][0]),int(hydro.waves['hr'][0]),int(hydro.waves['mm'][0]))).total_seconds())
                    Hs = np.array(hydro.waves['wvht (m)']).astype(float)
                    cumE = np.trapz((1/16)*1025*9.81*Hs[~np.isnan(Hs)]**2,np.array(t_sec)[~np.isnan(Hs)])/3600/1000000 #units are MJh/m^2 #
                    
                    # Calc TWL freeboard #
                    twl = hydro.calcTWL(Bf_here,exceedance=16)
                    twl_ref = hydro.calcTWL(0.1,exceedance=2)
                    F = max(twl[:,1])-toe_pre
                    
                    # Calc toe dx #
                    dx = scarpToes[scarp][1][ii][0] - scarpToes[scarp][0][ii][0]
                    
                    # Calc toe dx using change in contour position #
                    x_pre = scarpToes[scarp][0][ii][0]
                    # x_post = utils.transectElevIntercept(toe_pre, T_scarps[scarp][1][ii]['X'], T_scarps[scarp][1][ii]['Z'])
                    dx_contour = np.nan#x_post-x_pre
                    
                    # Calc toe dz #
                    dz = scarpToes[scarp][1][ii][2] - scarpToes[scarp][0][ii][2]
                    
                    # Calc pre-storm slope from toe to 5 m #
                    t = T_scarps[scarp][0][ii]
                    # fig,ax = plt.subplots(1)
                    # ax.plot(t['X'],t['Z'])
                    t = t[np.logical_and(t['X']<=scarpToes[scarp][0][ii][0],t['Z']<=5,t['Z']-scarpToes[scarp][0][ii][2]>=0)]
                    # t = t[np.logical_and(t['X']<=scarpToes[scarp][0][ii][0],t['Z']-scarpToes[scarp][0][ii][2]<=0.5,t['Z']-scarpToes[scarp][0][ii][2]>=0)]
                    m = -(t['Z'][0]-t['Z'][-1])/(t['X'][0]-t['X'][-1])
                    # inter,m = utils.linearRegression(t['X'],t['Z'])
                    slope = abs(m)
                    # ax.plot(t['X'],t['Z'],'c')
                    # ax.plot(t['X'],(t['X']*m)+inter,'r')
                    # ax.set_title(str(m))
                    # fig.show()
                    # plt.pause(2)
                    # plt.close('all')
                    
                    # Calc slope ratio #
                    slopeRatio = Bf_here/slope
                    
                    # Calc pre-storm beach width (toe to shoreline)
                    x_toe = scarpToes[scarp][0][ii][0]
                    try:
                        x_sl = utils.transectElevIntercept(0.36,T_scarps[scarp][0][ii]['X'],T_scarps[scarp][0][ii]['Z'])
                    except IndexError: # Error when transects does not reach MHW #
                        x_sl = np.nan
                    width = x_sl-x_toe
                    
                    # Calc impact duration #
                    impact_dur_i = np.where(twl[:,1]>scarpToes[scarp][0][ii][2])[0]
                    if len(impact_dur_i)==0:
                        impact_dur=0
                    elif len(impact_dur_i)==1:
                        impact_dur=0.5
                    else:
                        impact_dur = (len(impact_dur_i)-1)/2 # Assuming sampling is half hour #
                    
                    # Calc pre- and post-storm volume (beach) #
                    x_start = scarpToes[scarp][0][ii][0]
                    x_end = min((max(T_scarps[scarp][0][ii]['X']),max(T_scarps[scarp][1][ii]['X'])))
                    vol_beach1 = []
                    for d in range(0,2):
                        i1 = np.logical_and(T_scarps[scarp][d][ii]['X']>=x_start,T_scarps[scarp][d][ii]['X']<=x_end)
                        xx = T_scarps[scarp][d][ii]['X'][i1]
                        zz = T_scarps[scarp][d][ii]['Z'][i1]
                        vol = np.trapz(zz[~np.isnan(zz)],xx[~np.isnan(zz)])
                        vol_beach1.append(vol)
                    vol_beach = vol_beach1[1]-vol_beach1[0]
                    
                    
                    # Calc pre- and post-storm volume (dune) #
                # for ii in range(0,len(T_scarps[scarp][0])):                
                    x_end = scarpToes[scarp][0][ii][0]
                    x_start = max((min(T_scarps[scarp][0][ii]['X']),min(T_scarps[scarp][1][ii]['X'])))
                    vol_dune1 = []
                    for d in range(0,2):
                        i1 = np.logical_and(T_scarps[scarp][d][ii]['X']>=x_start,T_scarps[scarp][d][ii]['X']<=x_end)
                        xx = T_scarps[scarp][d][ii]['X'][i1]
                        zz = T_scarps[scarp][d][ii]['Z'][i1]
                        vol = np.trapz(zz[~np.isnan(zz)],xx[~np.isnan(zz)])
                        vol_dune1.append(vol)
                    vol_dune = vol_dune1[1]-vol_dune1[0]
                   
                    
                # fix,ax = plt.subplots(1)
                # ax.plot(T_scarps[scarp][0][ii]['X'],T_scarps[scarp][0][ii]['Z'],'k')
                # ax.plot(T_scarps[scarp][1][ii]['X'],T_scarps[scarp][1][ii]['Z'],'grey')
                # ax.plot(scarpToes[scarp][0][ii][0],scarpToes[scarp][0][ii][2],'ko')
                # ax.plot(scarpToes[scarp][1][ii][0],scarpToes[scarp][1][ii][2],'o',color='grey')
                # ax.plot((x_start,x_start),ax.get_ylim(),'r')
                # ax.plot((x_end,x_end),ax.get_ylim(),'r')
                # ax.set_title(vol_dune)
                # plt.show()
                # plt.pause(3)
                # plt.close('all')
                    
                    # fig,ax = plt.subplots(1)
                    # ax.plot(T_scarps[scarp][0][ii]['X'],T_scarps[scarp][0][ii]['Z'],'k')
                    # ax.plot(scarpToes[scarp][0][ii][0],scarpToes[scarp][0][ii][2],'k.')
                    # ax.plot(T_scarps[scarp][0][ii]['X'][np.where(abs(T_scarps[scarp][0][ii]['Z']-6) == min(abs(T_scarps[scarp][0][ii]['Z']-6)))],
                    #            T_scarps[scarp][0][ii]['Z'][np.where(abs(T_scarps[scarp][0][ii]['Z']-6) == min(abs(T_scarps[scarp][0][ii]['Z']-6)))],'r.')
                    # ax.plot(T_scarps[scarp][1][ii]['X'],T_scarps[scarp][1][ii]['Z'],'grey')
                    # ax.plot(scarpToes[scarp][1][ii][0],scarpToes[scarp][1][ii][2],'.',color='grey')
                    # ax.plot(T_scarps[scarp][1][ii]['X'][np.where(abs(T_scarps[scarp][1][ii]['Z']-6) == min(abs(T_scarps[scarp][1][ii]['Z']-6)))],
                    #            T_scarps[scarp][1][ii]['Z'][np.where(abs(T_scarps[scarp][1][ii]['Z']-6) == min(abs(T_scarps[scarp][1][ii]['Z']-6)))],'b.')
                    # fig.suptitle(vol_dune)
                    # plt.pause(5)
                    # plt.close('all')
                    
                    # Calc pre-storm dune height #
                    dHeight = max(T_scarps[scarp][0][ii]['Z'])-toe_pre
                    
             
                    results = results.append({'BT':BT_here,
                                              'Bf':Bf_here,
                                              '$Z_{Toe}$ (m)':toe_pre,
                                              'Post-storm toe elev (m)':toe_post,
                                              '$Bf_{post}$':Bf_post,
                                              'dx (m)':dx,
                                              'dx (m, contour)':dx_contour,
                                              'dz (m)':dz,
                                              '$B_{Dune}$':slope,
                                              'Bf/Bd ':slopeRatio,
                                              'F (m)':F,
                                              'W (m)':width,
                                              't_i (hr)':impact_dur,
                                              '$max(H_s)$ (m)':maxHs,
                                              '\sumE (MJh/m^2)':cumE,
                                              'TWL_ref (m)':twl_ref,
                                              '$\Delta V_B$ ($m^3/m$)':vol_beach,
                                              '$\Delta V_D$ ($m^3/m$)':vol_dune,
                                              'Dune Height (m)':dHeight},
                                             ignore_index=True)      
                
        if bdate==20130306:
            name = '2013NE'
        elif bdate==20180303:
            name = 'Riley'
        elif bdate==20190904:
            name = 'Dorian'
        elif bdate==20191011:
            name = '2019NE1'
        elif bdate==20191114:
            name = '2019NE2'
        elif bdate==20200910:
            name = 'Teddy'
        elif bdate==20210318:
            name = '2021NE'
            
        results_byStorm = results_byStorm.append({'Name':name,'results':results},ignore_index=True)

    
    names = list(results_byStorm['Name'])
    BT_all1 = []
    ii = -1
    vals = []
    iUse = []
    for i in names:
        ii+=1
        try:
            vals = results_byStorm.loc[results_byStorm['Name'] == i]['results'][ii]['BT']
            iUse1 = np.where(abs(vals)<=1)[0]
            vals = vals[iUse1]
        except:
            ii-=1
        else:
            BT_all1.append(vals)
            iUse.append(iUse1)
    BT_all = [item for sublist in BT_all1 for item in sublist]



    # Example up vs. down profiles #
    file_up = files_use[1]
    file_down = files_use[2]
     
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file_up,'rb')
    scarpResults=pickle.load(f)
    T_up = scarpResults.T_scarps
    toes_up = scarpResults.scarpToes
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file_down,'rb')
    scarpResults=pickle.load(f)
    T_down = scarpResults.T_scarps
    toes_down = scarpResults.scarpToes
      
    use_up = 5
    use_down = 100#40
    
    dv_ud = []
    for ud in ['up','down']:
        T = eval('T_'+ud)
        x_end = eval('toes_'+ud+'[0][0][use_'+ud+'][0]')
        x_start = max((min(T[0][0][eval('use_'+ud)]['X']),min(T[0][1][eval('use_'+ud)]['X'])))        
        vol_dune1 = []
        for d in range(0,2):
            i1 = np.logical_and(T[0][d][eval('use_'+ud)]['X']>=x_start,T[0][d][eval('use_'+ud)]['X']<=x_end)
            xx = T[0][d][eval('use_'+ud)]['X'][i1]
            zz = T[0][d][eval('use_'+ud)]['Z'][i1]
            vol = np.trapz(zz[~np.isnan(zz)],xx[~np.isnan(zz)])
            vol_dune1.append(vol)
        dv = vol_dune1[1]-vol_dune1[0]
        dv_ud.append(dv)
            
    
    dv_up = dv_ud[0]
    dv_down = dv_ud[1]
    
    fig,ax = plt.subplots(2,1,figsize=(4,4))
    h1 = ax[0].plot(T_up[0][0][use_up]['X'],T_up[0][0][use_up]['Z'],'k',linewidth=2)
    ax[0].plot(toes_up[0][0][use_up][0],toes_up[0][0][use_up][2],'ko')
    h2 = ax[0].plot(T_up[0][1][use_up]['X'],T_up[0][1][use_up]['Z'],'grey',linewidth=2)
    ax[0].plot(toes_up[0][1][use_up][0],toes_up[0][1][use_up][2],'o',color='grey')
    ax[1].plot(T_down[0][0][use_down]['X'],T_down[0][0][use_down]['Z'],'k',linewidth=2)
    ax[1].plot(toes_down[0][0][use_down][0],toes_down[0][0][use_down][2],'ko')
    ax[1].plot(T_down[0][1][use_down]['X'],T_down[0][1][use_down]['Z'],'grey',linewidth=2)
    ax[1].plot(toes_down[0][1][use_down][0],toes_down[0][1][use_down][2],'o',color='grey')  
    ax[1].set_xlabel('Cross-shore (m)')
    ax[1].set_ylabel('Elev. (m)')
    ax[0].set_ylabel('Elev. (m)')
    fig.legend([h1[0],h2[0]],['Pre-storm','Post-storm'],loc='upper center',ncol=2)
    
    ax[0].text(.6, .9,'Riley',ha='left', va='top',transform=ax[0].transAxes,fontsize=8,fontweight='bold')
    ax[0].text(.6, .8,'Y = '+str(int(T_up[0][0][use_up]['Y'][0])),ha='left', va='top',transform=ax[0].transAxes,fontsize=8,fontweight='bold')
    ax[0].text(.6, .7,'$\Delta V_D$ = '+str(round(dv_up,2))+' $m^3/m$',ha='left', va='top',transform=ax[0].transAxes,fontsize=8,fontweight='bold')
    ax[1].text(.6, .9,'Dorian',ha='left', va='top',transform=ax[1].transAxes,fontsize=8,fontweight='bold')
    ax[1].text(.6, .8,'Y = '+str(int(T_down[0][0][use_down]['Y'][0])),ha='left', va='top',transform=ax[1].transAxes,fontsize=8,fontweight='bold')
    ax[1].text(.6, .7,'$\Delta V_D$ = '+str(round(dv_down,2))+' $m^3/m$',ha='left', va='top',transform=ax[1].transAxes,fontsize=8,fontweight='bold')
    ax[0].text(0.03, 0.1, 'a', transform=ax[0].transAxes,fontsize=8, fontweight='bold', va='top',color='k')
    ax[1].text(0.03, 0.1, 'b', transform=ax[1].transAxes,fontsize=8, fontweight='bold', va='top',color='k')    
    
    # Histogram of a parameter for each storm #
    varsToPlot = ['$Z_{Toe}$ (m)'] 
    fig,ax = plt.subplots(3,3,figsize=(6,5))         
    for colorVal in varsToPlot:             
        # Do the binned plot and overall plot #
        vals_all1 =  []
        ii = -1
        for i in names:
            ii+=1
            try:
                vals = results_byStorm.loc[results_byStorm['Name'] == i]['results'][ii][colorVal]
            except:
                ii-=1
            else:
                vals_all1.append(vals)               
    axx = []
    for row in range(0,len(ax[:,0])):
        for col in range(0,len(ax[0,:])):
            axx.append(ax[row,col])
    for s in range(0,len(vals_all1)):       
        axx[s].hist(vals_all1[s])
        axx[s].set_title(results_byStorm['Name'][s],fontsize=12)
        if s==6:
            axx[s].set_xlabel(colorVal)
            axx[s].set_ylabel('Occurences')
    axx[7].set_visible(False)
    axx[8].set_visible(False)
   
    
    # Bar plot of a parameter for each storm # 
    colorVal = 'BT'
    fig,axxx = plt.subplots(1,figsize=(6.5,2))
    for ii in range(0,len(results_byStorm)):
        name = results_byStorm.loc[ii]['Name']
        results = results_byStorm.loc[ii]['results']
        val_mean = np.nanmedian(results[colorVal])
        val_25 = np.nanpercentile(results[colorVal],25)
        val_75 = np.nanpercentile(results[colorVal],75)
        
        axxx.add_artist((Rectangle((ii-0.4,0),0.8,val_mean,edgecolor='k')))
        axxx.plot((ii,ii),(val_mean,val_25),'k',linewidth=1)
        axxx.plot((ii,ii),(val_mean,val_75),'k',linewidth=1)        
        # axxx.plot(ii,val_mean,'k.',markersize=0.1)
    axxx.plot((-0.5,6.5),(0,0),'k')
    axxx.set_xlim(-0.5,6.5)
    axxx.set_xticks([0,1,2,3,4,5,6])
    axxx.set_xticklabels(names)
    axxx.set_ylabel(colorVal)
    # axxx.invert_yaxis()
    
    # Dune retreat vs. storm params #
    x1 = '\sumE (MJh/m^2)'
    x2 = 'TWL_ref (m)'
    y2 = '$\Delta V_D$ ($m^3/m$)'
    
    E = []
    twl = []
    L = []
    Vd = []
    Vd_err = []
    for ii in range(0,len(results_byStorm)):
        results = results_byStorm.loc[ii]['results']
        E.append(np.nanmedian(results[x1]))
        twl.append(max(results[x2][0][:,1]))
        L.append(len(results[x2])*5)
        Vd.append(np.nanmedian(results[y2])); Vd_err.append((np.nanpercentile(results[y2],25),np.nanpercentile(results[y2],75)))
        
    fig,ax = plt.subplots(2,2,figsize=(4,3))
    for i in range(0,len(E)):
        ax[0][0].plot(E[i],L[i],'s',markersize=16)
        ax[1][0].plot(E[i],Vd[i],'s',markersize=16),ax[1][0].invert_yaxis()
        ax[0][1].plot(twl[i],L[i],'s',markersize=16)
        ax[1][1].plot(twl[i],Vd[i],'s',markersize=16),ax[1][1].invert_yaxis()
    
    params,yhat,r2,p_values = utils.linearRegression(E,L)
    params,yhat,r2,p_values = utils.linearRegression(twl,L)
    params,yhat,r2,p_values = utils.linearRegression(E,Vd)
    params,yhat,r2,p_values = utils.linearRegression(twl,Vd)
      
                
    # Val vs. BT all scatter and storm-average #        
    varsToPlot = 'F (m)'  
    fig,ax = plt.subplots(1,2,figsize=(4,1.75))
    # ax = np.hstack([ax.reshape(-1,1),np.vstack([1,1,1])])
    plt.rcParams.update({'font.size': 8})
    plt.rc('axes', titlesize=8) 
    plt.rc('axes', labelsize=8)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    
    BT_mean_all = []
    val_mean_all = []
    BT_all = []
    val_all = []
    for ii in range(0,len(results_byStorm)):
        name = results_byStorm.loc[ii]['Name']
        results = results_byStorm.loc[ii]['results']
        BT_mean = np.nanmedian(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]])
        BT_25 = np.nanpercentile(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]],25)
        BT_75 = np.nanpercentile(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]],75)
        val_mean = np.nanmedian(np.array(results[varsToPlot]).astype(float)[np.where(abs(results['BT'])<1)[0]])
        val_25 = np.nanpercentile(np.array(results[varsToPlot]).astype(float)[np.where(abs(results['BT'])<1)[0]],25)
        val_75 = np.nanpercentile(np.array(results[varsToPlot]).astype(float)[np.where(abs(results['BT'])<1)[0]],75)
        
        alpha_scaled = 0.8#((n-20)/(100-20))
        
        ax[0].plot(np.array(results[varsToPlot])[np.where(abs(results['BT'])<1)[0]],
                   np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]],
                   '.',markersize=1)
        ax[0].set_xlabel(varsToPlot)
        ax[0].set_ylabel('BT')
        
        ax[1].plot( (val_mean,val_25),(BT_mean,BT_mean),'k',linewidth=1 )
        ax[1].plot( (val_mean,val_75),(BT_mean,BT_mean),'k',linewidth=1 )
        ax[1].plot( (val_mean,val_mean),(BT_mean,BT_25),'k',linewidth=1 )
        ax[1].plot( (val_mean,val_mean),(BT_mean,BT_75),'k',linewidth=1 )
        ax[1].plot(val_mean,BT_mean,'s',alpha=alpha_scaled,markersize=14,label=results_byStorm.loc[ii]['Name'])
        ax[1].set_xlabel(varsToPlot)
        
        val_mean_all.append(val_mean)
        BT_mean_all.append(BT_mean)
        val_all.append(np.array(results[varsToPlot])[np.where(abs(results['BT'])<1)[0]])
        BT_all.append(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]])
    
    ax[0].text(0.03, 0.1, 'a', transform=ax[0].transAxes,fontsize=8, fontweight='bold', va='top',color='k')
    ax[1].text(0.03, 0.1, 'b', transform=ax[1].transAxes,fontsize=8, fontweight='bold', va='top',color='k')
    ax[0].set_ylim(-1,0.5)
    # ax[0].set_xlim(-25,25)
    # ax[1].set_xlim(-11,0)
    
    fig.legend(loc='upper left',bbox_to_anchor=[0.78,0.91],frameon=False,handletextpad=0.1,fontsize=8,ncol=1,markerscale=0.5)

    val_all = np.hstack(val_all)
    BT_all = np.hstack(BT_all) 
    val_all = val_all[np.where(np.isnan(val_all)==0)[0]]
    BT_all = BT_all[np.where(np.isnan(val_all)==0)[0]]
    b,yhat,r2,p_values = utils.linearRegression(val_all,BT_all)
    
    b,yhat,r2,p_values = utils.linearRegression(np.array(val_mean_all),np.array(BT_mean_all))
    
    
    def power_law(x, a, b, c):
        return (a*np.power(x,b))+c
    def poly(x,a,b,c):
        return (a*x**2)+(b*x)+c
    def linear(x,m,b):
        return (m*x)+b
    
    pars, cov = curve_fit(f=power_law, xdata=np.array(val_mean_all)+2, ydata=np.array(BT_mean_all), p0=[0,-5,0], bounds=(-np.inf, np.inf))
    yhat_power =  (pars[0]*np.power(np.array(val_mean_all)+2,pars[1]))+pars[2]
    pars, cov = curve_fit(f=linear, xdata=np.array(val_mean_all)+2, ydata=np.array(BT_mean_all), p0=[0,0], bounds=(-np.inf, np.inf))
    yhat_linear = (pars[0]*(np.array(val_mean_all)+2))+pars[1]
    pars, cov = curve_fit(f=linear, xdata=np.array(val_mean_all)[0:-1]+2, ydata=np.array(BT_mean_all)[0:-1], p0=[0,0], bounds=(-np.inf, np.inf))
    yhat_linear_sub = (pars[0]*(np.array(val_mean_all)[0:-1]+2))+pars[1]
    
    fig,ax = plt.subplots(1)
    ax.plot(np.array(val_mean_all),np.array(BT_mean_all),'ks',markersize=16)
    ax.plot(np.array(val_mean_all)[np.argsort(np.array(val_mean_all))],yhat_power[np.argsort(np.array(val_mean_all))],'r',label='$y=ax^b+c$')
    ax.plot(np.array(val_mean_all)[np.argsort(np.array(val_mean_all))],yhat_linear[np.argsort(np.array(val_mean_all))],'b',label='$y=mx+b$')
    ax.plot(np.array(val_mean_all)[0:-1][np.argsort(np.array(val_mean_all)[0:-1])],yhat_linear_sub[np.argsort(np.array(val_mean_all)[0:-1])],'b--',label='$y=mx+b$')
    ax.set_xlabel('F (m)')
    ax.set_ylabel('BT')
    ax.legend()
            
    # F and ti vs. BT combined plot #
    varsToPlot1 = 'F (m)'
    varsToPlot2 = 't_i (hr)'
    
    BT_mean_all = []
    val_mean_all1 = []
    val_mean_all2 = []
    val_all1 = []
    val_all2 = []
    BT_all = []
    for ii in range(0,len(results_byStorm)):
        name = results_byStorm.loc[ii]['Name']
        results = results_byStorm.loc[ii]['results']
        BT_mean = np.nanmedian(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]])
        BT_25 = np.nanpercentile(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]],25)
        BT_75 = np.nanpercentile(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]],75)
        val_mean1 = np.nanmedian(np.array(results[varsToPlot1]).astype(float)[np.where(abs(results['BT'])<1)[0]])
        val_251 = np.nanpercentile(np.array(results[varsToPlot1]).astype(float)[np.where(abs(results['BT'])<1)[0]],25)
        val_751 = np.nanpercentile(np.array(results[varsToPlot1]).astype(float)[np.where(abs(results['BT'])<1)[0]],75)
        val_mean2 = np.nanmedian(np.array(results[varsToPlot2]).astype(float)[np.where(abs(results['BT'])<1)[0]])
        val_252 = np.nanpercentile(np.array(results[varsToPlot2]).astype(float)[np.where(abs(results['BT'])<1)[0]],25)
        val_752 = np.nanpercentile(np.array(results[varsToPlot2]).astype(float)[np.where(abs(results['BT'])<1)[0]],75)        
        
        val_mean_all1.append(val_mean1)
        val_mean_all2.append(val_mean2)
        BT_mean_all.append(BT_mean)
        
        val_all1.append(np.array(results[varsToPlot1])[np.where(abs(results['BT'])<1)[0]])
        val_all2.append(np.array(results[varsToPlot2])[np.where(abs(results['BT'])<1)[0]])
        BT_all.append(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]])
        
    val_all1s = np.hstack(val_all1)
    val_all2s = np.hstack(val_all2)
    BT_alls = np.hstack(BT_all)
        
    fig = plt.figure(figsize=(4,3))
    ax = plt.axes([0.15,0.15,0.65,0.75])
    cbax = plt.axes([0.82,0.3,0.03,0.4])
    h = ax.scatter(np.array(val_all1s),np.array(val_all2s),2,np.array(BT_alls),vmin=-.3,vmax=.3,cmap='turbo_r')
    ax.set_xlabel(varsToPlot1)
    ax.set_ylabel(varsToPlot2)
    plt.colorbar(h,cbax)
    ax.text(1.025,.75,'BT',transform=ax.transAxes,fontsize=10)
        
        
        
        
        
    # Vals vs BT plots #
    varsToPlot = ['t_i (hr)']  
    fig,ax = plt.subplots(3,len(varsToPlot),figsize=(3,5))
    if len(varsToPlot)==1:
        ax = np.hstack([ax.reshape(-1,1),np.vstack([1,1,1])])
    plt.rcParams.update({'font.size': 8})
    plt.rc('axes', titlesize=8) 
    plt.rc('axes', labelsize=8)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)

    iii = -1
    for colorVal in varsToPlot:
        iii+=1
              
        # Do the binned plot and overall plot #
        vals_all1 =  []
        ii = -1
        for i in names:
            ii+=1
            try:
                vals = results_byStorm.loc[results_byStorm['Name'] == i]['results'][ii][colorVal]
                vals = vals[iUse[ii]]
            except:
                ii-=1
            else:
                vals_all1.append(vals)
        vals_all = [item for sublist in vals_all1 for item in sublist]
        
        vals_all = np.array(vals_all)
        BT_all = np.array(BT_all)
        vals_all[np.abs(BT_all)>1] = np.nan
        BT_all[np.abs(BT_all)>1] = np.nan
        
        for s in range(0,len(vals_all1)):
            ax[1,iii].plot(vals_all1[s],BT_all1[s],'.',markersize=2)
        ax[1,iii].set_xlabel(colorVal)
        ax[1,iii].set_ylabel('BT')
        ax[1,iii].set_ylim(-1,1)
    

        
        
        
        # from sklearn.linear_model import LinearRegression
        # from scipy import stats
        # reg = LinearRegression().fit(np.array(vals_all)[~np.isnan(vals_all)].reshape(-1,1),
        #                              np.array(BT_all)[~np.isnan(vals_all)].reshape(-1,1))
        # R2 = reg.score(np.array(vals_all)[~np.isnan(vals_all)].reshape(-1,1),
        #                              np.array(BT_all)[~np.isnan(vals_all)].reshape(-1,1))
        # coef = reg.coef_
        # inter = reg.intercept_
        # # Calculate p-values. Method taken from top answer to this SO question https://stackoverflow.com/questions/27928275/find-p-value-significance-in-scikit-learn-linearregression #
        # x = np.array(vals_all)[~np.isnan(vals_all)].reshape(-1,1)
        # y = np.array(BT_all)[~np.isnan(vals_all)].reshape(-1,1)
        # params = np.append(inter[0],coef[0])
        # predictions = reg.predict(x)
        # newX = np.append(np.ones((len(x),1)), x.reshape(-1,1), axis=1)
        # MSE = (sum((y-predictions)**2))/(len(newX)-len(newX[0]))
        # var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
        # sd_b = np.sqrt(var_b)
        # ts_b = params/ sd_b
        # p_values =[2*(1-stats.t.cdf(np.abs(i),(len(newX)-len(newX[0])))) for i in ts_b] 
        
        
        # if p_values[1]>+0.01:
        #     t = ax[1,iii].text(min(ax[1,iii].get_xlim()),-0.5,'$R^2=$'+str(round(R2,2))+'\n$p=$'+str(round(p_values[1],3)),ha='left')
        # else:
        #     a = p_values[1]
        #     t = ax[1,iii].text(min(ax[1,iii].get_xlim()),-0.5,'$R^2=$'+str(round(R2,2))+'\n'+r"$p = {0:0.2e}$".format(a),ha='left')
        # t.set_bbox(dict(facecolor='red', alpha=0.3))
    
        
    
    
 
        
        
        
        
        
        from matplotlib.patches import Rectangle
        binsize = 0.1
        binedges = np.arange(-1,0.5,binsize)
        iBins = np.digitize(np.array(BT_all),binedges)
        iPos = np.where(np.array(BT_all)>0)[0]
        iNeg = np.where(np.array(BT_all)<0)[0]
               
        val_pos_mean = np.nanmean(np.array(vals_all)[iPos].astype(float))
        val_neg_mean = np.nanmean(np.array(vals_all)[iNeg].astype(float))
        
        # stat,p_val = stats.ks_2samp(np.array(vals_all)[iPos][~np.isnan(np.array(vals_all)[iPos])],np.array(vals_all)[iNeg][~np.isnan(np.array(vals_all)[iNeg])])
         
        # Manual KS test to check #
        # import random
        # binedges = np.linspace(min(vals_all),max(vals_all),100)
        # H_pos,X_pos = np.histogram( np.array(vals_all)[iPos][~np.isnan(np.array(vals_all)[iPos])], bins = binedges, density = True )
        # H_neg,X_neg = np.histogram( np.array(vals_all)[iNeg][~np.isnan(np.array(vals_all)[iNeg])], bins = binedges, density = True )
        # dx = X_pos[1] - X_pos[0]
        # F_pos = np.cumsum(H_pos)*dx
        # F_neg = np.cumsum(H_neg)*dx
        # ks_stat = max(abs(F_pos-F_neg))
        
        # N_pos = len(np.array(vals_all)[iPos][~np.isnan(np.array(vals_all)[iPos])])
        # N_neg = len(np.array(vals_all)[iNeg][~np.isnan(np.array(vals_all)[iNeg])])
        # vals_all_nonan = np.array(vals_all)[~np.isnan(vals_all)]
        # ks_stat_dist = []
        # for it in range(0,1000):
        #     i_pos = random.sample(range(0,len(vals_all_nonan)),N_pos)
        #     i_neg = [ii for ii in np.arange(0,len(vals_all_nonan)) if ii not in i_pos]
            
        #     H_pos1,X_pos1 = np.histogram( vals_all_nonan[i_pos], bins = binedges, density = True )
        #     H_neg1,X_neg1 = np.histogram( vals_all_nonan[i_neg], bins = binedges, density = True )
        #     dx1 = X_pos1[1] - X_pos1[0]
        #     F_pos1 = np.cumsum(H_pos1)*dx
        #     F_neg1 = np.cumsum(H_neg1)*dx
        #     ks_stat1 = max(abs(F_pos1-F_neg1))
        #     ks_stat_dist.append(ks_stat1)
        
  
        
        for binn in range(1,len(binedges)):
            val_mean = np.nanmean(np.array(vals_all)[np.where(iBins==binn)[0]].astype(float))
            val_std = np.nanstd(np.array(vals_all)[np.where(iBins==binn)[0]].astype(float))
            val_25 = np.nanpercentile(np.array(vals_all)[np.where(iBins==binn)[0]].astype(float),25)
            val_75 = np.nanpercentile(np.array(vals_all)[np.where(iBins==binn)[0]].astype(float),75)
           
            r = Rectangle((0,binedges[binn-1]),val_mean,binedges[binn]-binedges[binn-1],
                          edgecolor='k',facecolor='grey')
            ax[0,iii].add_patch(r)
            ax[0,iii].plot((val_mean,val_25),(binedges[binn-1]+((binedges[binn]-binedges[binn-1])/2),binedges[binn-1]+((binedges[binn]-binedges[binn-1])/2)),'k')
            ax[0,iii].plot((val_mean,val_75),(binedges[binn-1]+((binedges[binn]-binedges[binn-1])/2),binedges[binn-1]+((binedges[binn]-binedges[binn-1])/2)),'k')
            ax[0,iii].set_ylim(min(binedges),max(binedges))
            ax[0,iii].set_xlabel(colorVal)
            ax[0,iii].set_ylabel('BT')
            # ax[0,iii].text(0.033,((binedges[binn-1]-binedges[binn])/2)+binedges[binn],'n='+str(len(np.where(iBins==binn)[0])),va='center',ha='center')
            if colorVal==varsToPlot[-1]:
                ax[0,iii].text(65,((binedges[binn-1]-binedges[binn])/2)+binedges[binn],'n='+str(len(np.where(iBins==binn)[0])),va='center',ha='left',fontsize=7)
                
        ax[0,iii].set_xlim(ax[1,iii].get_xlim())
        m1 = ax[0,iii].plot(val_pos_mean,min(binedges),'b*',alpha=1)
        m2 = ax[0,iii].plot(val_neg_mean,min(binedges),'r*',alpha=1)
        m1[0].set_clip_on(False)
        m2[0].set_clip_on(False)
        fig.legend([m1[0],m2[0]],['Mean for (+)BT','Mean for (-)BT'],loc='upper center',frameon=False,ncol=2,handletextpad=0.1)
        
        # a = p_val
        # t = ax[0,iii].text(min(ax[0,iii].get_xlim()),-0.5,r"$p = {0:0.2e}$".format(a),ha='left')
        # t.set_bbox(dict(facecolor='red', alpha=0.3))
        
        
        # Do the storm-average plot #   
        BT_mean_all = []
        val_mean_all = []
        for ii in range(0,len(results_byStorm)):
            name = results_byStorm.loc[ii]['Name']
            results = results_byStorm.loc[ii]['results']
            BT_mean = np.nanmedian(results['BT'][np.where(results['BT']!=-np.inf)[0]])
            BT_25 = np.nanpercentile(results['BT'][np.where(results['BT']!=-np.inf)[0]],25)
            BT_75 = np.nanpercentile(results['BT'][np.where(results['BT']!=-np.inf)[0]],75)
            val_mean = np.nanmedian(results[colorVal])
            val_25 = np.nanpercentile(results[colorVal],25)
            val_75 = np.nanpercentile(results[colorVal],75)
            
            n = len(results[colorVal])
            alpha_scaled = 0.8#((n-20)/(100-20))
            # if alpha_scaled<0.1:
            #     alpha_scaled=0.1
            # print(alpha_scaled)
            
            ax[2,iii].plot( (val_mean,val_25),(BT_mean,BT_mean),'k',linewidth=1 )
            ax[2,iii].plot( (val_mean,val_75),(BT_mean,BT_mean),'k',linewidth=1 )
            ax[2,iii].plot( (val_mean,val_mean),(BT_mean,BT_25),'k',linewidth=1 )
            ax[2,iii].plot( (val_mean,val_mean),(BT_mean,BT_75),'k',linewidth=1 )
            ax[2,iii].plot(val_mean,BT_mean,'s',alpha=alpha_scaled,markersize=14,label=results_byStorm.loc[ii]['Name'])
            ax[2,iii].set_xlabel(colorVal)
            ax[2,iii].set_ylabel('BT')
            # if iii==0:
            # ax[0,iii].set_xlim(ax[1,iii].get_xlim())
            # else:
            ax[2,iii].set_xlim(ax[1,iii].get_xlim())
            
            val_mean_all.append(val_mean)
            BT_mean_all.append(BT_mean)
            
        # coefs,yhat,r2,p_values = utils.linearRegression(val_mean_all,BT_mean_all)
        
        # reg = LinearRegression().fit(np.array(val_mean_all).reshape(-1,1),np.array(BT_mean_all).reshape(-1,1))
        # R2 = reg.score(np.array(val_mean_all).reshape(-1,1),np.array(BT_mean_all).reshape(-1,1))
        # coef = reg.coef_
        # inter = reg.intercept_
        # # Calculate p-values. Method taken from top answer to this SO question https://stackoverflow.com/questions/27928275/find-p-value-significance-in-scikit-learn-linearregression #
        # x = np.array(val_mean_all).reshape(-1,1)
        # y = np.array(BT_mean_all).reshape(-1,1)
        # params = np.append(inter[0],coef[0])
        # predictions = reg.predict(x)
        # newX = np.append(np.ones((len(x),1)), x.reshape(-1,1), axis=1)
        # MSE = (sum((y-predictions)**2))/(len(newX)-len(newX[0]))
        # var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
        # sd_b = np.sqrt(var_b)
        # ts_b = params/ sd_b
        # p_values =[2*(1-stats.t.cdf(np.abs(i),(len(newX)-len(newX[0])))) for i in ts_b] 
        
        
        # if p_values[1]>+0.01:
        #     t = ax[2,iii].text(min(ax[1,iii].get_xlim()),-0.5,'$R^2=$'+str(round(R2,2))+'\n$p=$'+str(round(p_values[1],3)),ha='left')
        # else:
        #     a = p_values[1]
        #     t = ax[2,iii].text(min(ax[1,iii].get_xlim()),-0.5,'$R^2=$'+str(round(R2,2))+'\n'+r"$p = {0:0.2e}$".format(a),ha='left')
        # t.set_bbox(dict(facecolor='red', alpha=0.3))
            
        if colorVal==varsToPlot[0]:
            fig.legend(loc='lower center',frameon=False,handletextpad=0.1,fontsize=8,ncol=3)
            
        
        # Scatters for each storm #        
        fig,ax = plt.subplots(3,3,figsize=(6,5))         
        # fig,ax = plt.subplots(3,3,figsize=(8.3,7.8)) 
        # plt.rc('axes', labelsize=22)    # fontsize of the x and y labels
        # plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
        # plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
        # plt.rc('legend', fontsize=16) 

        axx = []
        for row in range(0,len(ax[:,0])):
            for col in range(0,len(ax[0,:])):
                axx.append(ax[row,col])
        for s in range(0,len(vals_all1)):
            iPos = np.where(np.sign(BT_all1[s])==1)[0]
            iNeg = np.where(np.sign(BT_all1[s])==-1)[0]
            axx[s].plot(vals_all1[s][iPos],BT_all1[s][iPos],'b.',markersize=2)
            axx[s].plot(vals_all1[s][iNeg],BT_all1[s][iNeg],'r.',markersize=2)            
            axx[s].set_title(results_byStorm['Name'][s],fontsize=12)
            if s==6:
                axx[s].set_xlabel(colorVal)
                axx[s].set_ylabel('$B_T$')
            params_pos,yhat_pos,r2_pos,p_values_pos = utils.linearRegression(vals_all1[s][iPos][np.logical_and(~np.isnan(vals_all1[s][iPos]),~np.isinf(BT_all1[s][iPos]))],
                                                              BT_all1[s][iPos][np.logical_and(~np.isnan(vals_all1[s][iPos]),~np.isinf(BT_all1[s][iPos]))])
            plt.text(.02, .34, 'm = '+str(round(params_pos[1],1)),ha='left', va='center',transform=axx[s].transAxes,color='b',fontsize=8)
            plt.text(.02, .22, '$R^2$ = '+str(round(r2_pos,1)),ha='left', va='center',transform=axx[s].transAxes,color='b',fontsize=8)
            plt.text(.02, .1, 'p = '+str(round(p_values_pos[1],2)),ha='left', va='center',transform=axx[s].transAxes,color='b',fontsize=8)
            
            try:
                params_neg,yhat_neg,r2_neg,p_values_neg = utils.linearRegression(vals_all1[s][iNeg][np.logical_and(~np.isnan(vals_all1[s][iNeg]),~np.isinf(BT_all1[s][iNeg]))],
                                                              BT_all1[s][iNeg][np.logical_and(~np.isnan(vals_all1[s][iNeg]),~np.isinf(BT_all1[s][iNeg]))])
            except:
                pass
            else:       
                plt.text(.98, .34, 'm = '+str(round(params_neg[1],1)),ha='right', va='center',transform=axx[s].transAxes,color='r',fontsize=8)
                plt.text(.98, .22, '$R^2$ = '+str(round(r2_neg,1)),ha='right', va='center',transform=axx[s].transAxes,color='r',fontsize=8)
                plt.text(.98, .1, 'p = '+str(round(p_values_neg[1],2)),ha='right', va='center',transform=axx[s].transAxes,color='r',fontsize=8)
                    
        for ii in range(0,len(results_byStorm)):
            name = results_byStorm.loc[ii]['Name']
            results = results_byStorm.loc[ii]['results']
            BT_mean = np.nanmedian(results['BT'][np.where(results['BT']!=-np.inf)[0]])
            BT_25 = np.nanpercentile(results['BT'][np.where(results['BT']!=-np.inf)[0]],25)
            BT_75 = np.nanpercentile(results['BT'][np.where(results['BT']!=-np.inf)[0]],75)
            val_mean = np.nanmedian(results[colorVal])
            val_25 = np.nanpercentile(results[colorVal],25)
            val_75 = np.nanpercentile(results[colorVal],75)
            
            n = len(results[colorVal])
            alpha_scaled = 0.8#((n-20)/(100-20))
            # if alpha_scaled<0.1:
            #     alpha_scaled=0.1
            # print(alpha_scaled)
            axx[-2].plot( (val_mean,val_25),(BT_mean,BT_mean),'k',linewidth=1 )
            axx[-2].plot( (val_mean,val_75),(BT_mean,BT_mean),'k',linewidth=1 )
            axx[-2].plot( (val_mean,val_mean),(BT_mean,BT_25),'k',linewidth=1 )
            axx[-2].plot( (val_mean,val_mean),(BT_mean,BT_75),'k',linewidth=1 )
            axx[-2].plot(val_mean,BT_mean,'s',alpha=alpha_scaled,markersize=14,label=results_byStorm.loc[ii]['Name'])
            axx[-2].set_xlabel(colorVal)
            
        axx[-1].spines['left'].set_visible(False)
        axx[-1].spines['right'].set_visible(False)
        axx[-1].spines['top'].set_visible(False)
        axx[-1].spines['bottom'].set_visible(False)
        axx[-1].set_xticks([])
        axx[-1].set_yticks([])

        fig.legend(loc='upper left',bbox_to_anchor=[0.67,0.26],frameon=False,handletextpad=0.1,fontsize=8,ncol=2)



def examineAndPlot_TWLVsBT():
    
    method = 'mc'
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if 'scarpResults' in i])
    files_use = [i for i in files if method in i]
    results_byStorm = pd.DataFrame(columns=['Name','results'])
    for file in files_use:
        
        results = pd.DataFrame(columns=['BT',
                                        '$TWL_{R2%}$','$TWL_{R7%}$','$TWL_{R16%}$','$TWL_{R50%}$',
                                        '$F_{R2%}$','$F_{R7%}$','$F_{R16%}$','$F_{R50%}$',
                                        '$ti_{R2%}$','$ti_{R7%}$','$ti_{R16%}$','$ti_{R50%}$'])
        bdate = int(file[0:8])
        edate = int(file[9:17])
        
        if '2013'  not in str(bdate):
            ws = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/frf_17m_waves.xlsx'
            hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=17.8)  
        else:  # Data from 17 m buoy missing for first few months of 2013. Use data from 26 m buoy and transform to 17 m with LWT #
            ws='/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/frf_26m_waves.xlsx'
            hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=26)  
            waves = hydro.waves                      
            for i in range(0,len(waves)):
                dat = waves.iloc[i]
                d = 26
                H = dat['wvht (m)']
                T = dat['Tp (s)']
                theta = dat['MWD (degT)']
                
                if 0<theta<180:               
                    d_vec = np.arange(-d,-16,1) # Vector of water depths #
                    # Wave parameters at buoy #  
                    k0 = utils.newtRaph(T,d)
                    C0 = ((2*np.pi)/k0)/T
                    n0 = 0.5*(1+((2*k0*d)/np.sinh(2*k0*d)))
                    alpha0 = 90-theta
                    Hh = H
                    # Transform #
                    for h in d_vec[1:len(d_vec)]:
                        k = utils.newtRaph(T,-h)
                        C = ((2*np.pi)/k)/T
                        n = 0.5*(1+((2*k*-h)/np.sinh(2*k*-h))) 
                        alpha = np.degrees(math.asin((C/C0)*np.sin(np.radians(alpha0))))
                        Hs = math.sqrt((n0*C0)/(n*C))*math.sqrt(math.cos(math.radians(alpha0))/math.cos(math.radians(alpha)))*Hh
                        
                        k0 = k
                        C0 = C
                        n0 = n
                        alpha0 = alpha
                        Hh = Hs
                        
                    waves.replace({'wvht (m)':waves.iloc[i]['wvht (m)'],'MWD (degT)':waves.iloc[i]['MWD (degT)']},
                                          {'wvht (m)':Hh,'MWD (degT)':360-(270+alpha)},inplace=True)
                else:
                    waves.replace({'wvht (m)':waves.iloc[i]['wvht (m)'],'MWD (degT)':waves.iloc[i]['MWD (degT)']},
                                          {'wvht (m)':np.nan,'MWD (degT)':np.nan},inplace=True)     
                    
            hydro.waves = waves
        

        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        scarpResults=pickle.load(f)
        Bf = scarpResults.Bf
        BT = scarpResults.BT
        scarpToes = scarpResults.scarpToes
        T_scarps = scarpResults.T_scarps
        forwardFlag = scarpResults.forwardFlag
        
        for scarp in range(0,len(T_scarps)):
            for ii in range(0,len(T_scarps[scarp][0])):
                
                if ii not in forwardFlag[scarp]:
                                
                    BT_here = BT[scarp][ii]
                    Bf_here = Bf[scarp][ii]
                    toe_pre = scarpToes[scarp][0][ii][2]
                    
                    F_all = []
                    ti_all = []
                    twl_all = []
                    for exc in [2,7,16,50]:
                        # Calc TWL freeboard #
                        twl = hydro.calcTWL(Bf_here,exceedance=exc)
                        F = max(twl[:,1])-toe_pre
                        F_all.append(F)
                        twl_all.append(max(twl[:,1]))
                        
    
                        # Calc impact duration #
                        impact_dur_i = np.where(twl[:,1]>scarpToes[scarp][0][ii][2])[0]
                        if len(impact_dur_i)==0:
                            impact_dur=0
                        elif len(impact_dur_i)==1:
                            impact_dur=0.5
                        else:
                            impact_dur = (len(impact_dur_i)-1)/2 # Assuming sampling is half hour #
                        ti_all.append(impact_dur)
                
                    results = results.append({'BT':BT_here,
                                              '$TWL_{R2%}$':twl_all[0],
                                              '$TWL_{R7%}$':twl_all[1],
                                              '$TWL_{R16%}$':twl_all[2],
                                              '$TWL_{R50%}$':twl_all[3],
                                              '$F_{R2%}$':F_all[0],
                                              '$F_{R7%}$':F_all[1],
                                              '$F_{R16%}$':F_all[2],
                                              '$F_{R50%}$':F_all[3],
                                              '$ti_{R2%}$':ti_all[0],
                                              '$ti_{R7%}$':ti_all[1],
                                              '$ti_{R16%}$':ti_all[2],
                                              '$ti_{R50%}$':ti_all[3]
                                              },
                                             ignore_index=True)      
                
        if bdate==20130306:
            name = '2013NE'
        elif bdate==20180303:
            name = 'Riley'
        elif bdate==20190904:
            name = 'Dorian'
        elif bdate==20191011:
            name = '2019NE1'
        elif bdate==20191114:
            name = '2019NE2'
        elif bdate==20200910:
            name = 'Teddy'
        elif bdate==20210318:
            name = '2021NE'
            
        results_byStorm = results_byStorm.append({'Name':name,'results':results},ignore_index=True)


    # Val vs. BT all scatter and storm-average #        
    BT_mean = []
    F_mean = []
    ti_mean = []
    BT_all = []
    F_all = []
    ti_all = []
    for exc in [2,7,16]:
    
        BT_mean1 = []
        F_mean1 = []
        ti_mean1 = []
        BT_all1 = []
        F_all1 = []
        ti_all1 = []
        for ii in range(0,len(results_byStorm)):
            name = results_byStorm.loc[ii]['Name']
            results = results_byStorm.loc[ii]['results']
            BT = np.nanmedian(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]])
            BT_25 = np.nanpercentile(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]],25)
            BT_75 = np.nanpercentile(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]],75)
            F = np.nanmedian(np.array(results['$F_{R'+str(exc)+'%}$'])[np.where(abs(results['BT'])<1)[0]])
            F_25 = np.nanpercentile(np.array(results['$F_{R'+str(exc)+'%}$'])[np.where(abs(results['BT'])<1)[0]],25)
            F_75 = np.nanpercentile(np.array(results['$F_{R'+str(exc)+'%}$'])[np.where(abs(results['BT'])<1)[0]],75)
            ti = np.nanmedian(np.array(results['$ti_{R'+str(exc)+'%}$'])[np.where(abs(results['BT'])<1)[0]])
            ti_25 = np.nanpercentile(np.array(results['$ti_{R'+str(exc)+'%}$'])[np.where(abs(results['BT'])<1)[0]],25)
            ti_75 = np.nanpercentile(np.array(results['$ti_{R'+str(exc)+'%}$'])[np.where(abs(results['BT'])<1)[0]],75) 
            
            BT_mean1.append([BT,(BT_25,BT_75)])
            F_mean1.append([F,(F_25,F_75)])
            ti_mean1.append([ti,(ti_25,ti_75)])
            
            BT_all1.append(np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]])
            F_all1.append(np.array(results['$F_{R'+str(exc)+'%}$'])[np.where(abs(results['BT'])<1)[0]])
            ti_all1.append(np.array(results['$ti_{R'+str(exc)+'%}$'])[np.where(abs(results['BT'])<1)[0]])
        
        BT_mean = BT_mean1
        F_mean.append(F_mean1)
        ti_mean.append(ti_mean1)
        BT_all.append(BT_all1)
        F_all.append(F_all1)
        ti_all.append(ti_all1)
        
    fig,ax = plt.subplots(3,2,figsize=(6.5,5))
    excs=[2,7,16]
    for i in range(0,len(excs)):
        [ax[i][0].plot(F_all[i][ii],BT_all[i][ii],'.',markersize=2) for ii in range(0,len(BT_all[0]))]
        hh = []
        for ii in range(0,len(F_mean[i])):
            ax[i][1].plot((F_mean[i][ii][0],F_mean[i][ii][0]),(BT_mean[ii][0],BT_mean[ii][1][0]),'k',linewidth=1)
            ax[i][1].plot((F_mean[i][ii][0],F_mean[i][ii][0]),(BT_mean[ii][0],BT_mean[ii][1][1]),'k',linewidth=1)
            ax[i][1].plot((F_mean[i][ii][0],F_mean[i][ii][1][0]),(BT_mean[ii][0],BT_mean[ii][0]),'k',linewidth=1)
            ax[i][1].plot((F_mean[i][ii][0],F_mean[i][ii][1][1]),(BT_mean[ii][0],BT_mean[ii][0]),'k',linewidth=1)
            hh.append(ax[i][1].plot(F_mean[i][ii][0],BT_mean[ii][0],'s',alpha=0.8,markersize=14))

            def power_law(x, a, b, c):
                return (a*np.power(x,b))+c

            xx = np.hstack([F_all[i][ii] for ii in range(0,len(BT_all[0]))])
            yy = np.hstack([BT_all[i][ii] for ii in range(0,len(BT_all[0]))])
            pars, cov = curve_fit(f=power_law, 
                              xdata=xx[~np.isnan(xx)]+5, 
                              ydata=yy[~np.isnan(xx)],
                              p0=[0,-5,0], 
                              bounds=(-np.inf, np.inf),
                              maxfev=10000)
            yhat_all =  (pars[0]*np.power(xx[~np.isnan(xx)]+5,pars[1]))+pars[2]

            sst = np.sum((yy[~np.isnan(xx)]-np.mean(yy[~np.isnan(xx)]))**2)
            ssr = np.sum((yhat_all-np.mean(yy[~np.isnan(xx)]))**2)
            r2_all = ssr/sst

            xx = np.array([F_mean[i][ii][0] for ii in range(0,len(F_mean[0]))])
            yy = np.array([BT_mean[ii][0] for ii in range(0,len(BT_mean))])
            pars, cov = curve_fit(f=power_law, 
                              xdata=xx[~np.isnan(xx)]+5, 
                              ydata=yy[~np.isnan(xx)],
                              p0=[0,-5,0], 
                              bounds=(-np.inf, np.inf),
                              maxfev=10000)
            yhat_mean =  (pars[0]*np.power(xx[~np.isnan(xx)]+5,pars[1]))+pars[2]

            sst = np.sum((yy[~np.isnan(xx)]-np.mean(yy[~np.isnan(xx)]))**2)
            ssr = np.sum((yhat_mean-np.mean(yy[~np.isnan(xx)]))**2)
            r2_mean = ssr/sst


            print((r2_all,r2_mean))



        ax[i,0].set_ylabel('$B_T$')
        ax[i,1].set_ylabel('$B_T$')
    fig.legend([ii[0] for ii in hh],['2013NE','Riley','Dorian','2019NE1','2019NE2','Teddy','2021NE'],loc='center right',
               frameon=False,handletextpad=0.1,fontsize=8,ncol=1,markerscale=0.5)
    ax[2,0].set_xlabel('$F_T$ (m)')
    ax[2,1].set_xlabel('$F_T$ (m)')
    ax[0,0].text(-.65,.5,'R2%',fontweight='bold',fontsize=10,transform=ax[0,0].transAxes)
    ax[1,0].text(-.65,.5,'R7%',fontweight='bold',fontsize=10,transform=ax[1,0].transAxes)
    ax[2,0].text(-.65,.5,'R16%',fontweight='bold',fontsize=10,transform=ax[2,0].transAxes)
    ax[0,0].text(0.01,0.9,'a',fontsize=8,fontweight='bold',transform=ax[0][0].transAxes)
    ax[0,1].text(0.01,0.9,'b',fontsize=8,fontweight='bold',transform=ax[0][1].transAxes)
    ax[1,0].text(0.01,0.9,'c',fontsize=8,fontweight='bold',transform=ax[1][0].transAxes)
    ax[1,1].text(0.01,0.9,'d',fontsize=8,fontweight='bold',transform=ax[1][1].transAxes)            
    ax[2,0].text(0.01,0.9,'e',fontsize=8,fontweight='bold',transform=ax[2][0].transAxes)
    ax[2,1].text(0.01,0.9,'f',fontsize=8,fontweight='bold',transform=ax[2][1].transAxes)                
    
    
    
    # TWL vs. BT this study and literature studies #
    BT_BL = [0.21,0.13]
    F_BL = [0.73,0.58]
    BT_SP = [0.05,0.11,-0.7,0.01,0.05]          
    F_SP = [1.95,2.07,0.13,1.52,2.11]
    BT_OL = [0.03,-0.01,-0.02]
    F_OL = [1.9,1,1.6]
    
    BT_all = []
    F_all = []
    for ii in range(0,len(results_byStorm)):
        name = results_byStorm.loc[ii]['Name']
        results = results_byStorm.loc[ii]['results']
        BT = np.array(results['BT'])[np.where(abs(results['BT'])<1)[0]]
        F = np.array(results['$F_{R2%}$'])[np.where(abs(results['BT'])<1)[0]]
        BT_all.append(BT)
        F_all.append(F)
    BT = np.hstack(BT_all)
    F = np.hstack(F_all)
    
    BT_all = np.hstack([BT,BT_BL,BT_SP,BT_OL])
    F_all = np.hstack([F,F_BL,F_SP,F_OL])
    
    def power_law(x, a, b, c):
        return (a*np.power(x,b))+c
    pars, cov = curve_fit(f=power_law, 
                      xdata=F_all[~np.isnan(F_all)]+5, 
                      ydata=BT_all[~np.isnan(F_all)],
                      p0=[0,-5,0], 
                      bounds=(-np.inf, np.inf),
                      maxfev=10000)
    yhat_power =  (pars[0]*np.power(F_all[~np.isnan(F_all)]+5,pars[1]))+pars[2]

    sst = np.sum((BT_all[~np.isnan(F_all)]-np.mean(BT_all[~np.isnan(F_all)]))**2)
    ssr = np.sum((yhat_power-np.mean(BT_all[~np.isnan(F_all)]))**2)
    r2_all_power = ssr/sst
    print(r2_all_power)
    
    fig = plt.figure(figsize=(4,3.5))
    ax = plt.axes([0.18,0.15,0.77,0.8])
    ax.plot(F,BT,'k.',markersize=1,label='This study')
    ax.plot(F_SP,BT_SP,'mo',markersize=4,label='Splinter & Palmsten (2012)')
    ax.plot(F_BL,BT_BL,'co',markersize=4,label='Bonte & Levoy (2015)')
    ax.plot(F_OL,BT_OL,'yo',markersize=4,label='Overbeck et al. (2017)')
    ax.plot(F_all[~np.isnan(F_all)][np.argsort(F_all[~np.isnan(F_all)])],yhat_power[np.argsort(F_all[~np.isnan(F_all)])],'r',linewidth=2,label='$y = '+str(round(pars[0],2))+'x^{'+str(round(pars[1],2))+'}+ '+str(round(pars[2],2))+'$')
    ax.legend(fontsize=8,loc='lower right',borderpad=0.1,handletextpad=0.1)
    ax.set_xlabel('$F_T$ (m)')
    ax.set_ylabel('$B_T$')
    
    
    
    
    # R2% vs. R50% #
    twlr = []
    for ii in range(0,len(results_byStorm)):
        name = results_byStorm.loc[ii]['Name']
        results = results_byStorm.loc[ii]['results']
        twl2 = np.array(results['$TWL_{R2%}$'])[np.where(abs(results['BT'])<1)[0]]
        twl50 = np.array(results['$TWL_{R50%}$'])[np.where(abs(results['BT'])<1)[0]]
        twlr.append(twl50/twl2)
    fig,ax = plt.subplots(1,figsize=(6.5,2))
    for i in range(0,len(twlr)):
        R = Rectangle((i+1-0.4,0),0.8,np.nanmedian(twlr[i]))
        ax.add_patch(R)
    ax.set_xlim(0.5,7.5)
    ax.set_xticklabels(['','2013NE','Riley','Dorian','2019NE1','2019NE2','Teddy','2021NE'])
    ax.set_ylabel()
 
    
        

def plotDuneChangeSummary():
    
    import datetime
    from matplotlib.colors import ListedColormap
    import matplotlib.dates as mdates
    myFmt = mdates.DateFormatter('%Y')

    def createRYBcmap():
        m = 256
        m1=m*0.5
        r = np.arange(0,m1)/(m1-1)
        g = np.arange(0,m1)/(m1-1)
        b = np.arange(0,m1)/(m1-1)
        r = np.vstack([r.reshape(-1,1),np.ones([len(r),1])])
        g = np.vstack([g.reshape(-1,1),np.flipud(g.reshape(-1,1))])
        b = np.vstack([b.reshape(-1,1),np.ones([len(b),1])])
        b = np.flipud(b)
        b = np.flipud(np.linspace(0,1,len(r))).reshape(-1,1)
        
        c = np.hstack([r,g,b,np.ones([len(r),1])])
        return np.flipud(c)
        
    c = createRYBcmap()
    ryb = ListedColormap(c)
    
    
    method = 'mc'
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if 'scarpResults' in i])
    files_use = [i for i in files if method in i]
    results_byStorm = pd.DataFrame(columns=['Name','results'])
    for file in files_use:
        
        results = pd.DataFrame(columns=['Y','BT','dx (m)','$\Delta V_D$ ($m^3/m$)'])
        bdate = int(file[0:8])
        edate = int(file[9:17])
     
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        scarpResults=pickle.load(f)
        Bf = scarpResults.Bf
        BT = scarpResults.BT
        scarpToes = scarpResults.scarpToes
        T_scarps = scarpResults.T_scarps
        forwardFlag = scarpResults.forwardFlag
        
        for scarp in range(0,len(T_scarps)):
            for ii in range(0,len(T_scarps[scarp][0])):
                
                if ii not in forwardFlag[scarp]:
                
                    BT_here = BT[scarp][ii]
                    
                    y_here = T_scarps[scarp][0][ii]['Y'][0]
                                      
                    # Calc toe dx #
                    dx = scarpToes[scarp][1][ii][0] - scarpToes[scarp][0][ii][0]
                               
                    # Calc pre- and post-storm volume (dune) #
                    x_end = scarpToes[scarp][0][ii][0]
                    x_start = max((min(T_scarps[scarp][0][ii]['X']),min(T_scarps[scarp][1][ii]['X'])))
                    vol_dune1 = []
                    for d in range(0,2):
                        i1 = np.logical_and(T_scarps[scarp][d][ii]['X']>=x_start,T_scarps[scarp][d][ii]['X']<=x_end)
                        xx = T_scarps[scarp][d][ii]['X'][i1]
                        zz = T_scarps[scarp][d][ii]['Z'][i1]
                        vol = np.trapz(zz[~np.isnan(zz)],xx[~np.isnan(zz)])
                        vol_dune1.append(vol)
                    vol_dune = vol_dune1[1]-vol_dune1[0]
                                     
             
                    results = results.append({'Y':y_here,
                                              'BT':BT_here,                                             
                                              'dx (m)':dx,                                              
                                              '$\Delta V_D$ ($m^3/m$)':vol_dune},
                                             ignore_index=True)      
                
        
        
        if bdate==20130306:
            name = '2013NE'
        elif bdate==20180303:
            name = 'Riley'
        elif bdate==20190904:
            name = 'Dorian'
        elif bdate==20191011:
            name = '2019NE1'
        elif bdate==20191114:
            name = '2019NE2'
        elif bdate==20200910:
            name = 'Teddy'
        elif bdate==20210318:
            name = '2021NE'
            
        results_byStorm = results_byStorm.append({'Name':name,'results':results},ignore_index=True)    
    
  
    
    
    yy = np.array([datetime.datetime(2011,1,1)+datetime.timedelta(days=i) for i in range(0,(365*11)+20,100)])
    xx = np.arange(-1000,4100,100)
    
    # Plot the storms with labels #
    fig = plt.figure(figsize=(3.25,6))
    plt.rcParams.update({'font.size': 8})
    ax1 = plt.axes([0.15,0.08,0.72,0.25])#plt.axes([0.1,0.08,0.72,0.25])
    ax1.set_xlim(0,4000)
    ax1.set_ylim(datetime.datetime(2013,1,1),datetime.datetime(2014,12,31))
    ax1.spines['top'].set_visible(False)
    ax1.set_yticks([datetime.datetime(2013,1,1),datetime.datetime(2014,1,1)])
    ax1.yaxis.set_major_formatter(myFmt)
    ax2 = plt.axes([0.15,0.33,0.72,0.62])#plt.axes([0.1,0.33,0.72,0.62])
    ax2.set_xlim(0,4000)
    ax2.set_ylim(datetime.datetime(2017,1,1),datetime.datetime(2021,12,13)) 
    ax2.spines['bottom'].set_visible(False)
    ax2.set_xticks([])
    ax2.set_yticks([datetime.datetime(2017,1,1),datetime.datetime(2018,1,1),datetime.datetime(2019,1,1),datetime.datetime(2020,1,1),datetime.datetime(2021,1,1)])
    ax2.yaxis.set_major_formatter(myFmt)
    
    
    d = .015  # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (.7-d, .7+d), **kwargs)        # top-left diagonal
    ax1.plot((1 - d, 1 + d), (.7-d, .7+d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (-.09-d, -.11+d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (-.09-d, -.11+d), **kwargs)  # bottom-right diagonal

    
    for ax in [ax1,ax2]:
        dates_storms = [datetime.datetime(2013,3,13),datetime.datetime(2017,9,25),
                        datetime.datetime(2018,3,6),datetime.datetime(2019,9,7),
                        datetime.datetime(2019,9,17),datetime.datetime(2019,10,12),
                        datetime.datetime(2019,11,16),
                        datetime.datetime(2020,9,17,12),datetime.datetime(2021,3,20,12)]
        for i in range(0,len(dates_storms)):
            ax.plot((min(xx),max(xx)),(dates_storms[i],dates_storms[i]),'k--',linewidth=0.5)
##            
##    ax1.text(4025,dates_storms[0],'NorEaster',va='center',ha='left',fontweight='bold',fontsize=8)
##    ax2.text(4025,dates_storms[1],'Maria',va='center',fontsize=8)
##    ax2.text(4025,dates_storms[2],'Riley',va='center',ha='left',fontweight='bold',fontsize=8)
##    ax2.text(4025,dates_storms[3]-datetime.timedelta(days=50),'Dorian',va='center',ha='left',fontweight='bold',fontsize=8)
##    ax2.text(4025,dates_storms[4]-datetime.timedelta(days=20),'Humberto',va='center',ha='left',fontsize=8)
##    ax2.text(4025,dates_storms[5]+datetime.timedelta(days=5),'NorEaster',va='center',ha='left',fontweight='bold',fontsize=8)
##    ax2.text(4025,dates_storms[6]+datetime.timedelta(days=15),'NorEaster',va='center',ha='left',fontweight='bold',fontsize=8)
##    ax2.text(4025,dates_storms[7],'Teddy',va='center',ha='left',fontweight='bold',fontsize=8)
##    ax2.text(4025,dates_storms[8],'NorEaster',va='center',ha='left',fontweight='bold',fontsize=8)
    ax1.set_xlabel('Alongshore (m)')
        
    dates_storms_use = [datetime.datetime(2013,3,13),
                        datetime.datetime(2018,3,6),datetime.datetime(2019,9,7),
                        datetime.datetime(2019,10,12),
                        datetime.datetime(2019,11,16),
                        datetime.datetime(2020,9,17,12),datetime.datetime(2021,3,20,12)]
    
    for i in range(0,len(dates_storms_use)):
        date_plot = dates_storms_use[i]
        
        # Move 2019 data to help visualization #
        if i==2:
            date_plot = date_plot-datetime.timedelta(days=40)
        elif i==4:
            date_plot = date_plot+datetime.timedelta(days=40)
            
        locs = np.array(results_byStorm.loc[i]['results']['Y'])
        vals = np.array(results_byStorm.loc[i]['results']['$\Delta V_D$ ($m^3/m$)'])
        times = np.tile(date_plot,len(locs))
            
        for ax in [ax1,ax2]:
            a = ax.scatter(locs,times,200,vals,marker='|',cmap='copper',vmin=-15,vmax=0,alpha=1)


    cbax = plt.axes([0.33,0.3,0.3,0.01])#plt.axes([0.33,0.3,0.3,0.01])
    cbax.set_xticks([])
    cbax.set_yticks([])
    plt.colorbar(a,cbax,orientation='horizontal',label='$\Delta V_D$ ($m^3/m$)',ticks=[-15,-10,-5,0])   
    
    ax2.text(0.02,0.97,'a',transform=ax2.transAxes,fontweight='bold')
##    ax2.set_yticklabels([])
##    ax1.set_yticklabels([])
    
    
    # Plot the bar plots #
    fig,ax = plt.subplots(4,1,figsize=(6.5,3),sharex=True)
    for ii in range(0,len(results_byStorm)):
        name = results_byStorm.loc[ii]['Name']
        results = results_byStorm.loc[ii]['results']
        
        L = len(results)*5
        ax[0].add_artist(Rectangle((ii-0.4,0),0.8,L))
        ax[0].plot(ii,L,'k.',markersize=0.01)
        
        V_mean = np.nanmedian(np.array(results['$\Delta V_D$ ($m^3/m$)'])[abs(results['BT']<1)])
        V_25 = np.nanpercentile(np.array(results['$\Delta V_D$ ($m^3/m$)'])[abs(results['BT']<1)],25)
        V_75 = np.nanpercentile(np.array(results['$\Delta V_D$ ($m^3/m$)'])[abs(results['BT']<1)],75)
        ax[1].add_artist(Rectangle((ii-0.4,0),0.8,V_mean))
        ax[1].plot((ii,ii),(V_mean,V_25),'k',linewidth=1)
        ax[1].plot((ii,ii),(V_mean,V_75),'k',linewidth=1)
        dat = np.array(results['$\Delta V_D$ ($m^3/m$)'])[abs(results['BT']<1)]
        iFliers = np.argsort(dat)[0:11]
##        ax[1].plot(np.tile(ii,len(iFliers)),dat[iFliers],'kx',markersize=2)
##        print(str(V_mean)+','+str(V_75-V_25)+','+str((V_75-V_25)/V_mean))
        
        dx_mean = np.nanmedian(np.array(results['dx (m)'])[abs(results['BT']<1)])
        dx_25 = np.nanpercentile(np.array(results['dx (m)'])[abs(results['BT']<1)],25)
        dx_75 = np.nanpercentile(np.array(results['dx (m)'])[abs(results['BT']<1)],75)  
        ax[2].add_artist(Rectangle((ii-0.4,0),0.8,dx_mean))
        ax[2].plot((ii,ii),(dx_mean,dx_25),'k',linewidth=1)
        ax[2].plot((ii,ii),(dx_mean,dx_75),'k',linewidth=1)

        BT_mean = np.nanmedian(np.array(results['BT'])[abs(results['BT']<1)])
        BT_25 = np.nanpercentile(np.array(results['BT'])[abs(results['BT']<1)],25)
        BT_75 = np.nanpercentile(np.array(results['BT'])[abs(results['BT']<1)],75)         
        ax[3].add_artist(Rectangle((ii-0.4,0),0.8,BT_mean))
        ax[3].plot((ii,ii),(BT_mean,BT_25),'k',linewidth=1)
        ax[3].plot((ii,ii),(BT_mean,BT_75),'k',linewidth=1)
        print(str(BT_mean)+','+str(BT_75-BT_25)+','+str((BT_75-BT_25)/BT_mean))

    ax[0].set_ylim(0,2000)
##    ax[1].set_ylim(-12,0)
    ax[2].set_ylim(-6,0)
    ax[3].set_ylim(-0.27,0.27)
    ax[1].invert_yaxis()
    ax[2].invert_yaxis()
    ax[0].set_ylabel('L (m)')
    ax[1].set_ylabel('$\Delta V_D$ ($m^3/m$)')
    ax[2].set_ylabel('dx (m)')
    ax[3].set_ylabel('$B_T$')
    for i in [0,1,2,3]:
        ax[i].plot((-1,8),(0,0),'k',linewidth=1)
        ax[i].set_xlim(-0.5,6.5)
    fig.align_ylabels()
    ax[3].set_xticklabels(['','2013NE','Riley','Dorian','2019NE1','2019NE2','Teddy','2021NE'])
    ax[0].text(0.005,0.8,'c',transform=ax[0].transAxes,fontweight='bold')
    ax[1].text(0.005,0.8,'d',transform=ax[1].transAxes,fontweight='bold')
    ax[2].text(0.005,0.8,'e',transform=ax[2].transAxes,fontweight='bold')
    ax[3].text(0.005,0.8,'f',transform=ax[3].transAxes,fontweight='bold')
        
    plt.show()  
 

def plotStormParams():   
    
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '5km_T' in i])
    files_use = [i for i in files if '2011' not in i and '20170918' not in i] # Get rid of Irene and Jose #
    results = pd.DataFrame(columns=['Name','max $H_s$ (m)',r'$\sum E\,(\frac{mJh}{m^2})$','$Z_{toe}$ (m)','$TWL_{max}$ (m)','max F (m)'])
    for file in files_use:
        
        bdate = int(file[0:8])
        edate = int(file[9:17])
        if '2013'  not in str(bdate):
            ws = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/frf_17m_waves.xlsx'
            hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=17.8)  
        else:  # Data from 17 m buoy missing for first few months of 2013. Use data from 26 m buoy and transform to 17 m with LWT #
            ws='/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/frf_26m_waves.xlsx'
            hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=26)  
            waves = hydro.waves                      
            for i in range(0,len(waves)):
                dat = waves.iloc[i]
                d = 26
                H = dat['wvht (m)']
                T = dat['Tp (s)']
                theta = dat['MWD (degT)']
                
                if 0<theta<180:               
                    d_vec = np.arange(-d,-16,1) # Vector of water depths #
                    # Wave parameters at buoy #  
                    k0 = utils.newtRaph(T,d)
                    C0 = ((2*np.pi)/k0)/T
                    n0 = 0.5*(1+((2*k0*d)/np.sinh(2*k0*d)))
                    alpha0 = 90-theta
                    Hh = H
                    # Transform #
                    for h in d_vec[1:len(d_vec)]:
                        k = utils.newtRaph(T,-h)
                        C = ((2*np.pi)/k)/T
                        n = 0.5*(1+((2*k*-h)/np.sinh(2*k*-h))) 
                        alpha = np.degrees(math.asin((C/C0)*np.sin(np.radians(alpha0))))
                        Hs = math.sqrt((n0*C0)/(n*C))*math.sqrt(math.cos(math.radians(alpha0))/math.cos(math.radians(alpha)))*Hh
                        
                        k0 = k
                        C0 = C
                        n0 = n
                        alpha0 = alpha
                        Hh = Hs
                        
                    waves.replace({'wvht (m)':waves.iloc[i]['wvht (m)'],'MWD (degT)':waves.iloc[i]['MWD (degT)']},
                                          {'wvht (m)':Hh,'MWD (degT)':360-(270+alpha)},inplace=True)
                else:
                    waves.replace({'wvht (m)':waves.iloc[i]['wvht (m)'],'MWD (degT)':waves.iloc[i]['MWD (degT)']},
                                          {'wvht (m)':np.nan,'MWD (degT)':np.nan},inplace=True)     
                    
            hydro.waves = waves
            
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        T = pickle.load(f)
        
        # Calc max wave height #
        maxHs = float(max(hydro.waves['wvht (m)']))
        
        # Calc cumulative wave energy between surveys #
        t_sec=[]
        for t in range(0,len(hydro.waves)):
            t_sec.append((datetime.datetime(int(hydro.waves['yr'][t]),int(hydro.waves['mo'][t]),int(hydro.waves['day'][t]),int(hydro.waves['hr'][t]),int(hydro.waves['mm'][t]))-
                      datetime.datetime(int(hydro.waves['yr'][0]),int(hydro.waves['mo'][0]),int(hydro.waves['day'][0]),int(hydro.waves['hr'][0]),int(hydro.waves['mm'][0]))).total_seconds())
        Hs = np.array(hydro.waves['wvht (m)']).astype(float)
        cumE = np.trapz((1/16)*1025*9.81*Hs**2,t_sec)/3600/1000000 #units are MJh/m^2 #
        cumEPert = cumE/((t_sec[-1]-t_sec[0])/60/60)
        
        # Calculate TWL and F using slope and toe at every transect #
        mgr = claris.scarpManager()
        toes = mgr.idDuneToe(T, method='mc_supervised')
        Bf = []
        for ii in range(0,len(toes[0])):
        
            x_beta = T[0][ii]['X'][T[0][ii]['X']>=toes[0][ii,0]]                       
            z_beta = T[0][ii]['Z'][T[0][ii]['X']>=toes[0][ii,0]]
            try:
                m = abs((z_beta[-1]-z_beta[0])/(x_beta[-1]-x_beta[0]))
            except:
                m = np.nan
            Bf.append(m)
        
        twl = []
        F =[]
        for Bf_here,toe_here in zip(Bf,toes[0]):
            if not np.isnan(Bf_here):
                twl1 = hydro.calcTWL(Bf_here,exceedance=2)
                F1 = max(twl1[:,1])-toe_here[2]
                twl.append(max(twl1[:,1]))
                F.append(F1)
            else:
                twl.append(np.nan)
                F.append(np.nan)
                
        # Get the name #                       
        if bdate==20130306:
            name = '2013NE'
        elif bdate==20170922:
            name = 'Maria'
        elif bdate==20180303:
            name = 'Riley'
        elif bdate==20190904:
            name = 'Dorian'
        elif bdate==20190910:
            name = 'Humberto'
        elif bdate==20191011:
            name = '2019NE1'
        elif bdate==20191114:
            name = '2019NE2'
        elif bdate==20200305:
            name = '2020NE'
        elif bdate==20200910:
            name = 'Teddy'
        elif bdate==20210318:
            name = '2021NE'
            
        results = results.append({'Name':name,
                       'max $H_s$ (m)':maxHs,
                       r'$\sum E\,(\frac{mJh}{m^2})$':cumE,
                       '$Z_{toe}$ (m)':toes[0][:,2],
                       '$TWL_{max}$ (m)':twl,
                       'max F (m)':F},
                       ignore_index=True)
                
    
    
    # Make the comparison figure (with Freeboard and colors # #
    fig,ax = plt.subplots(3,1,figsize=(6.5,4))
    c = -1
    for val in results.columns[[1,2,4]]:
        c+=1
        ax[c].set_ylabel(val)  
        
        if c!=2:       
            ax[c].bar(np.arange(0,len(results)),results[val],edgecolor='k')
            ax[c].set_xticks(np.arange(0,len(results)))
            ax[c].set_xticklabels([])
        else:
            ax[c].set_xticklabels(results['Name'],fontsize=8,rotation=-15)
            medians = [np.nanmedian(results.iloc[i][val]) for i in range(0,len(results))]
            ptile25 = [np.nanpercentile(results.iloc[i][val],25) for i in range(0,len(results))]
            ptile75 = [np.nanpercentile(results.iloc[i][val],75) for i in range(0,len(results))]
            ax[c].bar(np.arange(0,len(results)),medians,edgecolor='k')
            for i in range(0,len(results)):
                ax[c].plot((i,i),(medians[i],ptile25[i]),'k',linewidth=1)
                ax[c].plot((i,i),(medians[i],ptile75[i]),'k',linewidth=1)
            ax[c].set_xticks(np.arange(0,len(results)))
            xl = ax[c].get_xlim()
            ax[c].plot(xl,(0,0),'k',linewidth=.5)
            ax[c].set_xlim(xl)

    # Make the comparison figure (without Freeboard and colors # #
    fig,ax = plt.subplots(2,1,figsize=(6.5,3))
    c = -1
    for val in results.columns[[1,2]]:
        c+=1
        ax[c].set_ylabel(val)  
        
        if c!=1:       
            ax[c].bar(np.arange(0,len(results)),results[val],edgecolor='k')
            ax[c].set_xticks(np.arange(0,len(results)))
            ax[c].set_xticklabels([])
        else:
            ax[c].set_xticklabels(results['Name'],fontsize=8,rotation=-15)
            ax[c].bar(np.arange(0,len(results)),results[val],edgecolor='k')
            
            ax[c].set_xticks(np.arange(0,len(results)))
            xl = ax[c].get_xlim()
            ax[c].plot(xl,(0,0),'k',linewidth=.5)
            ax[c].set_xlim(xl)
    plt.show()
    
    
    
    
    def runAvg(x,y,window=10):
        xx = x[0]
        x_smooth=[]
        y_smooth=[]
        while xx in x:
            x_start = xx
            x_end = xx+window
            x_window = x[np.logical_and(x>=x_start,x<=x_end)]
            y_window = y[np.logical_and(x>=x_start,x<=x_end)]
            x_smooth.append(np.nanmean(x_window))
            y_smooth.append(np.nanmean(y_window))
            xx=x_end   
        return np.array(x_smooth),np.array(y_smooth)    
    
    names = ['2013NE','Maria','Riley','Dorian','Humberto','2019NE1','2019NE2','Teddy','2021NE']
    labs = ['i','h','g','f','e','d','c','b','a']
    fig,ax = plt.subplots(9,1,figsize=(5,6))
    method = 'mc'
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '5km_T' in i])
    files_use = [i for i in files if '2011' not in i and '20170918' not in i ]
    for i in range(0,len(files_use)):        
        
        storm = names[i]
        xx,ztoe = runAvg(np.arange(0,4000,5),results['$Z_{toe}$ (m)'][np.where(results['Name']==storm)[0][0]],window=50)
        xx,twl = runAvg(np.arange(0,4000,5),np.array(results['$TWL_{max}$ (m)'][np.where(results['Name']==storm)[0][0]]),window=50)
        h1 = ax[-i-1].plot(xx,ztoe,'k')
        h2 = ax[-i-1].plot(xx,twl,'b')
        ax[-i-1].set_ylim(1.5,5.5)
        
        try:
            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+files_use[i][0:17]+'_5km_scarpResults_mc_supervised.pkl','rb')
        except:
            pass
        else:
            scarpResults=pickle.load(f)
            T_scarps = scarpResults.T_scarps
            yl = ax[-i-1].get_ylim()
            difs_all = []
            for ii in range(0,len(T_scarps)):
                y_start = T_scarps[ii][0][0]['Y'][0]
                y_end = T_scarps[ii][0][-1]['Y'][0]
                iUse = np.where(np.logical_and(xx>y_start,xx<y_end))
                difs = (twl[iUse]-ztoe[iUse])/ztoe[iUse]
                difs_all.append(difs)
                r = Rectangle((y_start,yl[0]),y_end-y_start,np.diff(yl),facecolor='grey')
                ax[-i-1].add_artist(r)
        ax[-i-1].set_ylim(yl)
        ax[-i-1].set_xlim(0,4000)
        ax[-i-1].set_xticks([0,1000,2000,3000,4000])
        ax[-i-1].text(4050,3.5,labs[i]+' ('+names[i]+')',va='center',ha='left',fontsize=8,fontweight='bold')
        if -i-1 != -1:
            ax[-i-1].set_xticklabels([])
        ax[-1].set_xlabel('Alongshore (m)')
        ax[4].set_ylabel('Elev (m)')
        fig.legend([h1[0],h2[0]],['$Z_{toe}$','$TWL_{2\%}$'],loc='upper center',ncol=2,frameon=False)
             
    plt.show()  
    
    
    
    
    
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '5km_T' in i])
    files_use_T = [i for i in files if '2011' not in i and '20170918' not in i] # Get rid of Irene and Jose #
    i = -2
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+files_use_T[i],'rb')
    T = pickle.load(f)
    yy = [int(T[0][ii]['Y'][0]) for ii in range(0,len(T[0]))]
    y_want = 2425
    i_Use = np.where(np.array(yy)==y_want)[0][0]
    fig,ax = plt.subplots(1)
    ax.plot(T[0][i_Use]['X'],T[0][i_Use]['Z'],'k',linewidth=2)
    ax.plot(T[1][i_Use]['X'],T[1][i_Use]['Z'],'grey',linewidth=2)
    
    
    
    
    storm='Dorian'
    xx,ztoe = runAvg(np.arange(0,4000,5),results['$Z_{toe}$ (m)'][np.where(results['Name']==storm)[0][0]],window=50)
    xx,twl = runAvg(np.arange(0,4000,5),np.array(results['$TWL_{max}$ (m)'][np.where(results['Name']==storm)[0][0]]),window=50)
    fig,ax = plt.subplots(1)
    ax.plot(xx,ztoe,'k')
    ax.plot(xx,twl,'b')
    
    
    storm='2013NE'
    fig,ax = plt.subplots(2,1)
    ax[0].plot(np.arange(0,4000,5),results['max F (m)'][np.where(results['Name']==storm)[0][0]])
    ax[0].set_ylabel('F (m)')
    ax[1].plot(np.arange(0,4000,5),results['$Z_{toe}$ (m)'][np.where(results['Name']==storm)[0][0]],'k.',label='$Z_{toe}$ (m)')
    ax[1].plot(np.arange(0,4000,5),results['$TWL_{max}$ (m)'][np.where(results['Name']==storm)[0][0]],'b.',label='$TWL_{max}$ (m)')
    ax[1].set_ylabel('Elev (m)')
    ax[1].legend()
    ax[1].set_xlabel('FRF Y (m)')
    
    
    
    # AGU 2021 figure #
    fig = plt.figure(figsize=(9,7.4))
    plt.rc('axes', labelsize=22)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
    plt.rc('legend', fontsize=16) 
    
    ax1 = plt.axes([0.1,0.7,0.85,0.28])
    ax1.set_xticks(np.arange(0,len(results['Name'])))
    ax1.set_xticklabels(results['Name'],fontsize=18,rotation=-17)
    medians = [np.nanmedian(results.iloc[i]['max F (m)']) for i in range(0,len(results))]
    ptile25 = [np.nanpercentile(results.iloc[i]['max F (m)'],25) for i in range(0,len(results))]
    ptile75 = [np.nanpercentile(results.iloc[i]['max F (m)'],75) for i in range(0,len(results))]
    ax1.bar(np.arange(0,len(results)),medians,color=['b','grey','b','r','grey','b','r','grey','b','b'])
    for i in range(0,len(results)):
        ax1.plot((i,i),(medians[i],ptile25[i]),'k',linewidth=1)
        ax1.plot((i,i),(medians[i],ptile75[i]),'k',linewidth=1)
    ax[c].set_xticks(np.arange(0,len(results)))
    xl = ax[c].get_xlim()
    ax1.plot(xl,(0,0),'k',linewidth=.5)
    ax1.set_xlim(xl)
    ax1.set_ylabel('$F_{max}$ (m)')
    ax1.set_yticks([0,1])
    
    ax2 = plt.axes([0.1,0.1,0.25,0.15])
    storm='2019NE2'
    ax2.plot(np.arange(-1000,4000,5),results['max F (m)'][np.where(results['Name']==storm)[0][0]])
    ax2.set_ylabel('F (m)')
    ax2.set_xlabel('alongshore (m)')   
    # plt.text(.02, .95, storm,ha='left', va='top',transform=ax2.transAxes,fontsize=18)
    ax2.set_title(storm,fontsize=18)
    
    ax3 = plt.axes([0.4,0.1,0.25,0.15])
    storm='Teddy'
    ax3.plot(np.arange(-1000,4000,5),results['max F (m)'][np.where(results['Name']==storm)[0][0]])
    # plt.text(.02, .95, storm,ha='left', va='top',transform=ax3.transAxes,fontsize=18)
    ax3.set_title(storm,fontsize=18)      
     
    ax4 = plt.axes([0.7,0.1,0.25,0.15])
    storm='2021NE'
    ax4.plot(np.arange(-1000,4000,5),results['max F (m)'][np.where(results['Name']==storm)[0][0]])
    # plt.text(.02, .95, storm,ha='left', va='top',transform=ax4.transAxes,fontsize=18)
    ax4.set_title(storm,fontsize=18)
       
    ax5 = plt.axes([0.1,0.37,0.25,0.15])
    storm='2013NE'
    ax5.plot(np.arange(-1000,4000,5),results['max F (m)'][np.where(results['Name']==storm)[0][0]])
    # plt.text(.02, .95, storm,ha='left', va='top',transform=ax4.transAxes,fontsize=18)
    ax5.set_title(storm,fontsize=18)
    
    ax6 = plt.axes([0.4,0.37,0.25,0.15])
    storm='Riley'
    ax6.plot(np.arange(-1000,4000,5),results['max F (m)'][np.where(results['Name']==storm)[0][0]])
    # plt.text(.02, .95, storm,ha='left', va='top',transform=ax4.transAxes,fontsize=18)
    ax6.set_title(storm,fontsize=18)

    ax7 = plt.axes([0.7,0.37,0.25,0.15])
    storm='Dorian'
    ax7.plot(np.arange(-1000,4000,5),results['max F (m)'][np.where(results['Name']==storm)[0][0]])
    # plt.text(.02, .95, storm,ha='left', va='top',transform=ax4.transAxes,fontsize=18)
    ax7.set_title(storm,fontsize=18)
        
    
def plotHydroTimeseries():     
    
    import math
    
    # Get the wave and water level data #
    waves_all = []
    wl_all = []
    twl_all = []
    for yr in np.arange(2013,2022):
        if yr != 2013:
            ws = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/frf_17m_waves.xlsx'
            if yr!=2021:
                hydro = utils.hydroLab(int(str(yr)+'0101'),int(str(yr)+'1231'),station_wl=8651370,station_waves=ws,buoyDepth=17.8)  
            else:
                hydro = utils.hydroLab(int(str(yr)+'0101'),int(str(yr)+'0520'),station_wl=8651370,station_waves=ws,buoyDepth=17.8)              
            waves_all.append(hydro.waves)
            wl_all.append(hydro.wl)
            twl_all.append(hydro.calcTWL(beta=0.1,exceedance=2))
        else:  # Data from 17 m buoy missing for first few months of 2013. Use data from 26 m buoy and transform to 17 m with LWT #
            ws='/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/frf_26m_waves.xlsx'
            hydro = utils.hydroLab(int(str(yr)+'0101'),int(str(yr)+'1231'),station_wl=8651370,station_waves=ws,buoyDepth=26)  
            waves = hydro.waves
                       
            for i in range(0,len(waves)):
                dat = waves.iloc[i]
                d = 26
                H = dat['wvht (m)']
                T = dat['Tp (s)']
                theta = dat['MWD (degT)']
                
                if 0<theta<180:               
                    d_vec = np.arange(-d,-16,1) # Vector of water depths #
                    # Wave parameters at buoy #  
                    k0 = utils.newtRaph(T,d)
                    C0 = ((2*np.pi)/k0)/T
                    n0 = 0.5*(1+((2*k0*d)/np.sinh(2*k0*d)))
                    alpha0 = 90-theta
                    Hh = H
                    # Transform #
                    for h in d_vec[1:len(d_vec)]:
                        k = utils.newtRaph(T,-h)
                        C = ((2*np.pi)/k)/T
                        n = 0.5*(1+((2*k*-h)/np.sinh(2*k*-h))) 
                        alpha = np.degrees(math.asin((C/C0)*np.sin(np.radians(alpha0))))
                        Hs = math.sqrt((n0*C0)/(n*C))*math.sqrt(math.cos(math.radians(alpha0))/math.cos(math.radians(alpha)))*Hh
                        
                        k0 = k
                        C0 = C
                        n0 = n
                        alpha0 = alpha
                        Hh = Hs
                        
                    waves.replace({'wvht (m)':waves.iloc[i]['wvht (m)'],'MWD (degT)':waves.iloc[i]['MWD (degT)']},
                                          {'wvht (m)':Hh,'MWD (degT)':360-(270+alpha)},inplace=True)
                else:
                    waves.replace({'wvht (m)':waves.iloc[i]['wvht (m)'],'MWD (degT)':waves.iloc[i]['MWD (degT)']},
                                          {'wvht (m)':np.nan,'MWD (degT)':np.nan},inplace=True)     
                    
            waves_all.append(waves)
            wl_all.append(hydro.wl)
            hydro.waves = waves
            twl_all.append(hydro.calcTWL(beta=0.1,exceedance=2))
                
            
        
    waves = pd.concat(waves_all,ignore_index=True)
    wl = pd.concat(wl_all,ignore_index=True)
    twl = np.vstack(twl_all)
    
    dtime_waves = [datetime.datetime(int(waves.iloc[i]['yr']),int(waves.iloc[i]['mo']),
                                     int(waves.iloc[i]['day']),int(waves.iloc[i]['hr']),
                                     int(waves.iloc[i]['mm'])) for i in range(0,len(waves))]
    dtime_twl = twl[:,0]
        
    
    # Get the dates of the surveys
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '5km_T' in i])
    files_use = [i for i in files if '2011' not in i and '20170918' not in i] # Get rid of Irene and Jose #
    dates_all = []
    for file in files_use:       
        bdate = int(file[0:8])
        edate = int(file[9:17])
        dates = (bdate,edate)
        dates_all.append(dates)
    names = ['2013NE','Maria','Riley','Dorian','Humberto','2019NE1','2019NE2','Teddy','2021NE']   
        
    # Make the plot #
    fig,ax = plt.subplots(4,1,figsize=(6.5,5),sharex=True)
    ax[0].plot(dtime_waves,np.array(waves['wvht (m)']).astype(float),'.',markersize=1)
    ax[1].plot(dtime_waves,np.array(waves['Tp (s)']).astype(float),'.',markersize=1)
    ax[2].plot(dtime_waves,np.array(waves['MWD (degT)']).astype(float),'.',markersize=1)
    ax[3].plot(dtime_twl,twl[:,1],'.',markersize=2)
    ax[0].set_ylim(0,6)
    ax[1].set_ylim(0,25)
    ax[2].set_ylim(0,180);ax[2].set_yticks([0,90,180])
    ax[0].set_ylabel('$H_s$ (m)')
    ax[1].set_ylabel('$T_p$ (s)')
    ax[2].set_ylabel('D ($^o$)')
    ax[3].set_ylabel('TWL (m)')
      
    for ii in np.arange(0,4):
        for i in range(0,len(names)):
            yl = ax[ii].get_ylim()
            bdate = dates_all[i][0]
            edate = dates_all[i][1]
            x_start = datetime.datetime(int(str(bdate)[0:4]),int(str(bdate)[4:6]),int(str(bdate)[6:8]))
            x_end =  datetime.datetime(int(str(edate)[0:4]),int(str(edate)[4:6]),int(str(edate)[6:8]))
            dx = x_end-x_start
            rr = Rectangle((x_start,yl[0]),dx,np.diff(yl),edgecolor='r',facecolor='none')
            ax[ii].add_artist(rr)
            ax[ii].set_ylim(yl)
            
    # Make the separate subplots plot #
    fig,ax = plt.subplots(4,9,figsize=(6.5,4.5))
    twl_max = []
    for i in range(0,len(ax[0])):
        ax[0][i].plot(dtime_waves,np.array(waves['wvht (m)']).astype(float),'.',markersize=1)
        ax[1][i].plot(dtime_waves,np.array(waves['Tp (s)']).astype(float),'.',markersize=1)
        ax[2][i].plot(dtime_waves,np.array(waves['MWD (degT)']).astype(float),'.',markersize=1)
        ax[3][i].plot(dtime_twl,twl[:,1],'.',markersize=1)
        
        bdate = dates_all[i][0]
        edate = dates_all[i][1]
        x1 = datetime.datetime(int(str(bdate)[0:4]),int(str(bdate)[4:6]),int(str(bdate)[6:8]))
        x2 =  datetime.datetime(int(str(edate)[0:4]),int(str(edate)[4:6]),int(str(edate)[6:8]))               
        ax[0][i].set_xlim(x1,x2)
        ax[1][i].set_xlim(x1,x2)
        ax[2][i].set_xlim(x1,x2)
        ax[3][i].set_xlim(x1,x2)
        
        twl_max.append(np.nanmax(twl[:,1][np.logical_and(np.array(dtime_twl)>x1,np.array(dtime_twl)<x2)]))


        ax[0][i].set_ylim(0,6)
        ax[1][i].set_ylim(0,20)
        ax[2][i].set_ylim(0,180);ax[2][i].set_yticks([0,90,180])
        ax[3][i].set_ylim(0,4.5);ax[3][i].set_yticks([0,2,4])
       
        if i!=0:
            ax[0][i].set_yticklabels([])
            ax[1][i].set_yticklabels([])
            ax[2][i].set_yticklabels([])
            ax[3][i].set_yticklabels([])
        else:
            ax[0][i].text(x1-datetime.timedelta(15),ax[0][i].get_ylim()[0]+np.diff(ax[0][i].get_ylim())/2,'$H_s$ (m)',rotation=90,
                          ha='center',va='center')
            ax[1][i].text(x1-datetime.timedelta(15),ax[1][i].get_ylim()[0]+np.diff(ax[1][i].get_ylim())/2,'$T_p$ (s)',rotation=90,
                          ha='center',va='center')
            ax[2][i].text(x1-datetime.timedelta(15),ax[2][i].get_ylim()[0]+np.diff(ax[2][i].get_ylim())/2,'D ($^o$)',rotation=90,
                          ha='center',va='center')
            ax[3][i].text(x1-datetime.timedelta(15),ax[3][i].get_ylim()[0]+np.diff(ax[3][i].get_ylim())/2,'$TWL_{2\%}$ (m)',rotation=90,
                          ha='center',va='center')
            
        
        ax[0][i].set_xticks([x1+((x2-x1)/2)])
        ax[1][i].set_xticks([x1+((x2-x1)/2)])
        ax[2][i].set_xticks([x1+((x2-x1)/2)])
        ax[3][i].set_xticks([x1+((x2-x1)/2)])

        ax[0][i].set_xticklabels([])
        ax[1][i].set_xticklabels([])
        ax[2][i].set_xticklabels([])
        ax[3][i].set_xticklabels([str(int((x2-x1).days))+' days'],rotation=15)
        
        ax[0][i].set_title(names[i],fontweight='normal',rotation=90,fontsize=8)

    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/Hydrodynamics.png',dpi=400)
    plt.show()
            
    
 
    
 
    
 
    
 
    
 
    

def plotScarpPositionMap_Hovmoller(): 
    
    import datetime
    from matplotlib import cm 
    from matplotlib.colors import ListedColormap,LinearSegmentedColormap
    from matplotlib import patches
    
    def createRYBcmap():
        m = 256
        m1=m*0.5
        r = np.arange(0,m1)/(m1-1)
        g = np.arange(0,m1)/(m1-1)
        b = np.arange(0,m1)/(m1-1)
        r = np.vstack([r.reshape(-1,1),np.ones([len(r),1])])
        g = np.vstack([g.reshape(-1,1),np.flipud(g.reshape(-1,1))])
        b = np.vstack([b.reshape(-1,1),np.ones([len(b),1])])
        b = np.flipud(b)
        b = np.flipud(np.linspace(0,1,len(r))).reshape(-1,1)
        
        c = np.hstack([r,g,b,np.ones([len(r),1])])
        return np.flipud(c)
        
    c = createRYBcmap()
    ryb = ListedColormap(c)
    
    
    yy = np.array([datetime.datetime(2011,1,1)+datetime.timedelta(days=i) for i in range(0,(365*11)+20,100)])
    xx = np.arange(-1000,4100,100)
    
    fig = plt.figure(figsize=(4,6))
    plt.rcParams.update({'font.size': 8})
    ax = plt.axes([0.1,0.08,0.5,0.88])
    ax.set_xlim(-1000,4000)
    ax.set_ylim(datetime.datetime(2011,1,1),datetime.datetime(2021,12,13))
    ax.set_xlabel('FRF Y (m)')
    
    
            
    # Plot all storms for which pre+post storm data exist #
    dates_storms = [datetime.datetime(2011,8,27),datetime.datetime(2013,3,13),
                    datetime.datetime(2017,9,20),datetime.datetime(2017,9,25),
                    datetime.datetime(2018,3,6),datetime.datetime(2019,9,7),
                    datetime.datetime(2019,9,17),datetime.datetime(2019,10,13),
                    datetime.datetime(2019,11,16),datetime.datetime(2020,9,17,12),
                    datetime.datetime(2021,3,20,12)]
    names_storms = ['Irene','NorEaster','Jose','Maria','Riley','Dorian','Humberto','NorEaster','NorEaster','Teddy','NorEaster']
    for i in range(0,len(names_storms)):
        ax.plot((min(xx),max(xx)),(dates_storms[i],dates_storms[i]),'k--',linewidth=0.5)
    ax.text(4050,dates_storms[0],'Irene',va='center',ha='left')
    ax.text(4050,dates_storms[1],'NorEaster',va='center',ha='left')
    ax.text(4050,dates_storms[2]-datetime.timedelta(days=30),'Jose',va='center')
    ax.text(4050,dates_storms[3]+datetime.timedelta(days=30),'Maria',va='center')
    ax.text(4050,dates_storms[4],'Riley',va='center',ha='left')
    ax.text(4050,dates_storms[5]-datetime.timedelta(days=80),'Dorian',va='center',ha='left')
    ax.text(4050,dates_storms[6]-datetime.timedelta(days=20),'Humberto',va='center',ha='left')
    ax.text(4050,dates_storms[7]+datetime.timedelta(days=30),'NorEaster',va='center',ha='left')
    ax.text(4050,dates_storms[8]+datetime.timedelta(days=70),'NorEaster',va='center',ha='left')
    ax.text(4050,dates_storms[9],'Teddy',va='center',ha='left')
    ax.text(4050,dates_storms[10],'NorEaster',va='center',ha='left')
    
    for i in range(0,len(dates)):
        data = results.loc[results['Date']==dates[i]]
        # data = data.loc[data['Method']=='_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m']
        data = data.loc[data['Method']=='pd']

        # data = data.loc[data['Method']=='pd']
        if len(data)>0:
            numScarps = max(data['Scarp'])
            
            # file = [iii for iii in files if str(dates[i])+'-' in iii and '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in iii][0]
            file = [iii for iii in files if str(dates[i])+'-' in iii and 'pd' in iii][0]

            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
            scarpResults=pickle.load(f)
            scarpToes = scarpResults.scarpToes
            BT = scarpResults.BT
            T_scarps = scarpResults.T_scarps
            
            
            date_pre = datetime.datetime(int(str(dates[i])[0:4]),int(str(dates[i])[4:6]),int(str(dates[i])[6:8]))
            date_post = datetime.datetime(int(str(dates_post[i])[0:4]),int(str(dates_post[i])[4:6]),int(str(dates_post[i])[6:8]))
            date_plot = date_pre+((date_post-date_pre)/2)
            
            for scarp in range(0,numScarps):
                try:
                    loc = np.arange(min(scarpToes[scarp][0][:,1]),max(scarpToes[scarp][0][:,1]))
                    BT_region = np.mean(BT[scarp])
                    xx = loc
                    yy = np.array([date_plot+datetime.timedelta(i) for i in range(-20,21)])
                    a = ax.pcolor(xx,yy,np.tile(BT_region,[len(yy),len(xx)]),cmap=ryb,vmin=-.2,vmax=.2,alpha=1)
                    
                except:
                    pass
     
    cbax = plt.axes([0.26,0.45,0.2,0.01])
    cbax.set_xticks([])
    cbax.set_yticks([])
    plt.colorbar(a,cbax,orientation='horizontal',label='$B_T$')

 

           
def plotScarpPositionMap_Hovmoller_withExampleProfiles(): 
    
    import datetime
    from matplotlib import cm 
    from matplotlib.colors import ListedColormap,LinearSegmentedColormap
    from matplotlib import patches
    
    def createRYBcmap():
        m = 256
        m1=m*0.5
        r = np.arange(0,m1)/(m1-1)
        g = np.arange(0,m1)/(m1-1)
        b = np.arange(0,m1)/(m1-1)
        r = np.vstack([r.reshape(-1,1),np.ones([len(r),1])])
        g = np.vstack([g.reshape(-1,1),np.flipud(g.reshape(-1,1))])
        b = np.vstack([b.reshape(-1,1),np.ones([len(b),1])])
        b = np.flipud(b)
        b = np.flipud(np.linspace(0,1,len(r))).reshape(-1,1)
        
        c = np.hstack([r,g,b,np.ones([len(r),1])])
        return np.flipud(c)
        
    c = createRYBcmap()
    ryb = ListedColormap(c)
    
    # Create custom colormap #
    # top = cm.get_cmap('autumn', 128) # r means reversed version
    # bottom = cm.get_cmap('summer_r', 128)
    # newcolors = np.vstack((top(np.linspace(0, 1, 128)),
    #                    bottom(np.linspace(0, 1, 128))))
    # autumn_summer = ListedColormap(newcolors, name='AutumnSummer')
    
    yy = np.array([datetime.datetime(2011,1,1)+datetime.timedelta(days=i) for i in range(0,(365*11)+20,100)])
    xx = np.arange(-1000,4100,100)
    
    fig = plt.figure(figsize=(4,6))
    plt.rcParams.update({'font.size': 8})
    ax = plt.axes([0.1,0.08,0.5,0.88])
    ax.set_xlim(-1000,4000)
    ax.set_ylim(datetime.datetime(2011,1,1),datetime.datetime(2021,12,13))
    ax.set_xlabel('FRF Y (m)')
    
    
            
    # Plot all storms for which pre+post storm data exist #
    dates_storms = [datetime.datetime(2011,8,27),datetime.datetime(2013,3,13),
                    datetime.datetime(2017,9,20),datetime.datetime(2017,9,25),
                    datetime.datetime(2018,3,6),datetime.datetime(2019,9,7),
                    datetime.datetime(2019,9,17),datetime.datetime(2019,10,13),
                    datetime.datetime(2019,11,16),datetime.datetime(2020,9,17,12),
                    datetime.datetime(2021,3,20,12)]
    names_storms = ['Irene','NorEaster','Jose','Maria','Riley','Dorian','Humberto','NorEaster','NorEaster','Teddy','NorEaster']
    for i in range(0,len(names_storms)):
        ax.plot((min(xx),max(xx)),(dates_storms[i],dates_storms[i]),'k--',linewidth=0.5)
    ax.text(4050,dates_storms[0],'Irene',va='center',ha='left')
    ax.text(4050,dates_storms[1],'NorEaster',va='center',ha='left')
    ax.text(4050,dates_storms[2]-datetime.timedelta(days=30),'Jose',va='center')
    ax.text(4050,dates_storms[3]+datetime.timedelta(days=30),'Maria',va='center')
    ax.text(4050,dates_storms[4],'Riley',va='center',ha='left')
    ax.text(4050,dates_storms[5]-datetime.timedelta(days=80),'Dorian',va='center',ha='left')
    ax.text(4050,dates_storms[6]-datetime.timedelta(days=20),'Humberto',va='center',ha='left')
    ax.text(4050,dates_storms[7]+datetime.timedelta(days=30),'NorEaster',va='center',ha='left')
    ax.text(4050,dates_storms[8]+datetime.timedelta(days=70),'NorEaster',va='center',ha='left')
    ax.text(4050,dates_storms[9],'Teddy',va='center',ha='left')
    ax.text(4050,dates_storms[10],'NorEaster',va='center',ha='left')
    
    scarpUs = [99,4,1,0,0,2,1,0,0]
    for i in range(0,len(dates)):
        data = results.loc[results['Date']==dates[i]]
        # data = data.loc[data['Method']=='_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m']
        data = data.loc[data['Method']=='pd']

        # data = data.loc[data['Method']=='pd']
        if len(data)>0:
            numScarps = max(data['Scarp'])
            
            # file = [iii for iii in files if str(dates[i])+'-' in iii and '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in iii][0]
            file = [iii for iii in files if str(dates[i])+'-' in iii and 'pd' in iii][0]

            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
            scarpResults=pickle.load(f)
            scarpToes = scarpResults.scarpToes
            BT = scarpResults.BT
            T_scarps = scarpResults.T_scarps
            
            width = .15
            height = .1
            axsub = plt.axes([.8,(i*height+.05),width,height])
            axsub.set_ylim(0,8)
            
            lens = np.array([len(q) for q in BT])
            scarpU = scarpUs[i]
            tran_pre = T_scarps[scarpU][0][round(len(T_scarps[scarpU][0])/2)]
            tran_post = T_scarps[scarpU][1][round(len(T_scarps[scarpU][1])/2)]
            dToe_pre = scarpToes[scarpU][0][round(len(T_scarps[scarpU][0])/2)]
            dToe_post = scarpToes[scarpU][1][round(len(T_scarps[scarpU][0])/2)]
            hh1 = axsub.plot(tran_pre['X'],tran_pre['Z'],'k')
            hh2 = axsub.plot(tran_post['X'],tran_post['Z'],color='grey')
            axsub.plot(dToe_pre[0],dToe_pre[2],'k.')
            axsub.plot(dToe_post[0],dToe_post[2],'.',color='grey')
            axsub.set_xlim(min(axsub.get_xlim())-10,min(axsub.get_xlim())+70)
            axsub.set_xticks(axsub.get_xlim())
            if i==1:
                axsub.set_xticklabels([0,80])
                axsub.set_xlabel('x (m)',labelpad=0.1)
                axsub.set_ylabel('z (m)')
            else:
                axsub.set_xticklabels([])
                axsub.set_yticklabels([])
            
            date_pre = datetime.datetime(int(str(dates[i])[0:4]),int(str(dates[i])[4:6]),int(str(dates[i])[6:8]))
            date_post = datetime.datetime(int(str(dates_post[i])[0:4]),int(str(dates_post[i])[4:6]),int(str(dates_post[i])[6:8]))
            date_plot = date_pre+((date_post-date_pre)/2)
            
            for scarp in range(0,numScarps):
                try:
                    loc = np.arange(min(scarpToes[scarp][0][:,1]),max(scarpToes[scarp][0][:,1]))
                    BT_region = np.mean(BT[scarp])
                    xx = loc
                    yy = np.array([date_plot+datetime.timedelta(i) for i in range(-20,21)])
                    a = ax.pcolor(xx,yy,np.tile(BT_region,[len(yy),len(xx)]),cmap=ryb,vmin=-.2,vmax=.2,alpha=1)
                    
                except:
                    pass
     
    cbax = plt.axes([0.26,0.45,0.2,0.01])
    cbax.set_xticks([])
    cbax.set_yticks([])
    plt.colorbar(a,cbax,orientation='horizontal',label='$B_T$')
    fig.legend([hh1[0],hh2[0]],['Pre','Post'],ncol=2,loc=4,columnspacing=0.4,handletextpad=0.2)
            
        
            
            

def plotManualSummary():
    
    from datetime import datetime
    
    fig = plt.figure(figsize=(6.5,5))
    ax = plt.axes([0.1,0.6,0.85,0.35])
    plt.rcParams.update({'font.size': 8})
    ax.set_xlim(0,len(dates)+1)
    ax.set_ylim(-0.4,0.4)
    ax.plot((0,len(dates)+1),(0,0),'k')
    col = cm.get_cmap('cividis', 6).colors
    
    subLeftRegionWidth = 1/len(dates)
    subLeftRegionStarts = np.arange(0,1,subLeftRegionWidth)
    
    transectsUse = [[10,10],[10,15],[7,20,7],[10,10],[10]]
    labelV = ['a','b','c','d','e']
    
    x = 0
    for date in dates:
        x+=1
        data = results.loc[results['Date']==date]
        numScarps = max(data['Scarp'])
        if numScarps==6:
            xEdges = [-0.65,-0.45,-0.25,0,0.25,0.45]
        elif numScarps==3:
            xEdges=[-0.35,-0.05,0.25]
        elif numScarps==2:
            xEdges = [-0.3,0]
        elif numScarps==1:
            xEdges=[-0.15]
            
        for scarp in range(1,numScarps+1):
            data1 = data.loc[data['Scarp']==scarp]
            data2 = data1.loc[data1['Method']=='manual_transects_0.5mx0.5m']
            
            p = Rectangle((x+xEdges[scarp-1],0),0.3,float(data2['u']),facecolor=col[x-1,:],edgecolor='k')
            ax.add_patch(p)
            ax.plot( (x+xEdges[scarp-1]+0.15,x+xEdges[scarp-1]+0.15),(float(data2['u']),float(data2['u'])+float(data2['sigma'])),'k',linewidth=1 )
            ax.plot( (x+xEdges[scarp-1]+0.15,x+xEdges[scarp-1]+0.15),(float(data2['u']),float(data2['u'])-float(data2['sigma'])),'k',linewidth=1 )
            
            axs = plt.axes([subLeftRegionStarts[x-1]+.05,0.35-(0.15*(scarp-1)),subLeftRegionWidth-0.08,0.15])
            
            file = [iii for iii in files if str(date)+'-' in iii and 'manual_transects' in iii][0]
            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
            scarpResults=pickle.load(f)
            T_scarps = scarpResults.T_scarps
            scarpToes = scarpResults.scarpToes
            h1 = axs.plot(T_scarps[scarp-1][0][transectsUse[x-1][scarp-1]]['X'],T_scarps[scarp-1][0][transectsUse[x-1][scarp-1]]['Z'],'k')
            axs.plot(scarpToes[scarp-1][0][transectsUse[x-1][scarp-1]][0],scarpToes[scarp-1][0][transectsUse[x-1][scarp-1]][2],'k.')      
            h2 = axs.plot(T_scarps[scarp-1][1][transectsUse[x-1][scarp-1]]['X'],T_scarps[scarp-1][1][transectsUse[x-1][scarp-1]]['Z'],'grey')
            axs.plot(scarpToes[scarp-1][1][transectsUse[x-1][scarp-1]][0],scarpToes[scarp-1][1][transectsUse[x-1][scarp-1]][2],'.',color='grey')      
            axs.set_ylim(0,8)
            if scarp==1:
                axs.text(min(axs.get_xlim()),max(axs.get_ylim())+0.1,labelV[x-1],fontweight='bold')                
            if date==dates[-1]:
                fig.legend([h1[0],h2[0]],['Pre-storm','Post-storm'],loc='lower right')
                axs.set_xlabel('FRF X (m)')
                axs.set_ylabel('Z (m)')
        
        ax.set_xticks(np.arange(1,len(dates)+1))  
        xls = [datetime(int(str(dates[i])[0:4]),int(str(dates[i])[4:6]),int(str(dates[i])[6:8])).strftime("%B %Y") for i in range(0,len(dates))]
        ax.set_xticklabels(xls,rotation=10) 
        ax.set_ylabel('$B_T$')
        ax.text(1,0.25,'a',ha='center',fontweight='bold')
        ax.text(2,0.1,'b',ha='center',fontweight='bold')
        ax.text(3,0.1,'c',ha='center',fontweight='bold')
        ax.text(4,0.2,'d',ha='center',fontweight='bold')
        ax.text(5,0.15,'e',ha='center',fontweight='bold')
        # plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/manualMethodCompResults.png',dpi=350)
        


            
def plotMethodComp():
    totalScarps = sum([max(results.loc[results['Date']==i]['Scarp']) for i in dates])  
    
    fig = plt.figure(figsize=(8,3))
    ax = plt.axes([0.1,0.15,0.85,0.8])
    plt.rcParams.update({'font.size': 8})
    ax.set_xlim(0,totalScarps+1)
    ax.set_ylim(-1,1)
    ax.plot((0,totalScarps),(0,0),'k')
    col = cm.get_cmap('viridis', 6).colors
    
    c=0
    xls=[]
    for date in dates:
        scarp=1
        a=True
        while a is True:
            c+=1
            vals1 = results.loc[results['Date']==date]
            vals = vals1.loc[vals1['Scarp']==scarp]
            m = -0.3
            count=0
            h=[]
            for method in ['manual_pc','manual_transects','pd','mc','rr','ml']:
                dat = vals.loc[vals['Method']==method]
                if len(dat)>0:
                    # The bar #
                    p = Rectangle((c+m,0),width=0.1,height=float(dat['u']),color=col[count,:],ec='k')
                    h.append(ax.add_patch(p))
                    # The std bar #
                    ax.plot((c+(m+0.05),c+(m+0.05)),(float(dat['u']),float(dat['u'])+float(dat['sigma'])),'k',linewidth=1)
                    ax.plot((c+(m+0.05),c+(m+0.05)),(float(dat['u']),float(dat['u'])-float(dat['sigma'])),'k',linewidth=1)
                else:
                    pass
                    
                m+=0.1
                count+=1
            xls.append(str(date)+'\nscarp'+str(scarp))
            scarp+=1    
            a = max(vals1['Scarp'])>=scarp            
            
    ax.legend(h,['man. (pc)','man. (prof)','pd','mc','rr','ml'])
    ax.set_xticks(np.arange(1,totalScarps+1))
    ax.set_xticklabels(xls)
    ax.set_ylabel('$B_T')
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/methodCompResults.png',dpi=350)
    
        
    
def plotManualComp_20190904_20190910_Longshore():
    fig = plt.figure(figsize=(3,5))
    ax1 = plt.axes([0.2,0.4,0.25,0.5])
    ax2 = plt.axes([0.55,0.4,0.25,0.5])
    ax3 = plt.axes([0.2,0.1,0.25,0.2])
    ax4 = plt.axes([0.55,0.1,0.25,0.2])
    plt.rcParams.update({'font.size': 8})
    col = cm.get_cmap('viridis', 6).colors
    
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/20190904-20190910_5km_scarpResults_manual_pc.pkl','rb')
    pc = pickle.load(f)
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/20190904-20190910_5km_scarpResults_manual_transects.pkl','rb')
    tran = pickle.load(f)
    
    toes_pc_1_pre = pc.scarpToes[0][0]
    toes_pc_1_post = pc.scarpToes[0][1]
    toes_pc_2_pre = pc.scarpToes[1][0]
    toes_pc_2_post = pc.scarpToes[1][1]
    
    toes_tran_1_pre = tran.scarpToes[0][0]
    toes_tran_1_post = tran.scarpToes[0][1]
    toes_tran_2_pre = tran.scarpToes[1][0]
    toes_tran_2_post = tran.scarpToes[1][1]    
    
    ax1.plot(toes_pc_1_pre[:,0],toes_pc_1_pre[:,1],'b-')
    ax1.plot(toes_pc_1_post[:,0],toes_pc_1_post[:,1],'r-')  
    ax1.plot(toes_tran_1_pre[:,0],toes_tran_1_pre[:,1],'b--')
    ax1.plot(toes_tran_1_post[:,0],toes_tran_1_post[:,1],'r--') 
    ax1.text(45,1190,'Scarp 1')
    ax1.set_xlabel('FRF X (m)')
    ax1.set_ylabel('FRF Y (m)')
    
    ax2.plot(toes_pc_2_pre[:,0],toes_pc_2_pre[:,1],'b-')
    ax2.plot(toes_pc_2_post[:,0],toes_pc_2_post[:,1],'r-')  
    ax2.plot(toes_tran_2_pre[:,0],toes_tran_2_pre[:,1],'b--')
    ax2.plot(toes_tran_2_post[:,0],toes_tran_2_post[:,1],'r--')  
    ax2.text(-45,3800,'Scarp 2')
    ax2.set_xlabel('FRF X (m)')
    ax2.set_ylabel('FRF Y (m)')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    
    ax1.set_xlim(ax1.get_xlim())
    ax1.set_ylim(ax1.get_ylim())
    h1 = ax1.plot((-100,-99),(-100,-99),'b-')
    h2 = ax1.plot((-100,-99),(-100,-99),'b--')
    h3 = ax1.plot((-100,-99),(-100,-99),'r-')
    h4 = ax1.plot((-100,-99),(-100,-99),'r--')    
    fig.legend([h1[0],h2[0],h3[0],h4[0]],
               ['Pre-storm, PC','Pre-storm, prof','Post-storm, PC','Post-storm, prof'],
               loc='upper center' ,ncol=2,frameon=False,labelspacing=0.1,columnspacing=1)
    
    ax3.hist([pc.BTBf[0],tran.BTBf[0]],bins=np.linspace(-6,1,20),label=['PC','Prof'],density=True)
    ax3.legend(loc='upper left',handletextpad=0.01)
    ax3.set_xlabel('$B_T/B_f$')
    ax3.set_ylim(0,1)
    
    
    ax4.hist([pc.BTBf[1],tran.BTBf[1]],bins=np.linspace(-1,3,20),density=True)
    ax4.yaxis.tick_right()
    ax4.yaxis.set_label_position("right")
    ax4.set_xlabel('$B_T/B_f$')
    ax4.set_ylim(0,1)
        
    plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/20190904-20190910_5km_manualToePositionComp_Longshore.png',dpi=350)
    fig.show()

    
def examineAndPlot_TWLNov2019():
        
    file = [iii for iii in files if str(20191114)+'-' in iii and 'manual_transects' in iii][0]
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
    pc = pickle.load(f)    
    
    scarpToes = pc.scarpToes
    T_scarps = pc.T_scarps
    BT = pc.BT
    Bf = pc.Bf
    file = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/2019_17m.txt' 
    hydro = utils.hydroLab(20191114,20191119,station_wl=8651370,station_waves=file,buoyDepth=17.8)
    
    for scarp in range(0,len(Bf)):
        for tran in range(0,len(Bf[scarp])):
            twl = hydro.calcTWL(Bf[scarp][tran],exceedance=2)




def examineAndPlot_PreStormProfileVolumeAndBTBf():
    
    vals_all = []
    n=0
    for d in range(0,len(dates)):
        
        file = [iii for iii in files if str(dates[d])+'-' in iii and 'manual_transects' in iii][0]
    
    
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        pc = pickle.load(f)
        
        
        for ii in range(0,len(pc.T_scarps)):
            n+=1
            vs = []
            for i in range(0,len(pc.T_scarps[ii][0])):
                toe = pc.scarpToes[ii][0][i]
                z_high = toe[2]
                z_low=0.5
                
                x = pc.T_scarps[ii][0][i]['X']
                z = pc.T_scarps[ii][0][i]['Z']
                x_vol = x[np.logical_and(z<=z_high,z>=z_low)]
                z_vol = z[np.logical_and(z<=z_high,z>=z_low)]
                v = np.trapz(z_vol,x_vol)
                vs.append(v)
                
               
                fig,ax = plt.subplots(1)
                ax.plot(x,z)
                ax.plot(x_vol,z_vol,'r')
                ax.plot((min(x_vol),max(x_vol)),(z_low,z_low),'k')
                ax.plot((min(x_vol),min(x_vol)),(z_low,z_high),'k')
                ax.set_title('vol = '+str(round(v,3)))
                plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/BeachfaceVolVsBTBf_example'+str(n)+'.png',dpi=350)
                fig.show()
                plt.pause(2)
                plt.close('all')
                
            vals = np.stack([vs,pc.BTBf[ii]])
            vals_all.append(vals)
       
    fig,ax = plt.subplots(round(len(vals_all)/2),2)
    axx = []
    for i in ax:
        for ii in range(0,2):
            axx.append(i[ii])
    for s in range(0,len(vals_all)):
        axx[s].plot(vals_all[s][0],vals_all[s][1],'.')
        axx[s].axis('equal')
    axx[len(vals_all)-1].set_xlabel('V ($m^3/m$)')
    axx[len(vals_all)-1].set_ylabel('$B_T/B_f$')
 
       
def examineAndPlot_PreStormToeElevAndBTBf():
    
    vals_all = []
    n=0
    for d in range(0,len(dates)):
        
        file = [iii for iii in files if str(dates[d])+'-' in iii and 'manual_transects' in iii][0]
    
    
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        pc = pickle.load(f)
        
        
        for ii in range(0,len(pc.T_scarps)):
            n+=1
            elevs = []
            for i in range(0,len(pc.T_scarps[ii][0])):
                toe = pc.scarpToes[ii][0][i]
                toe_elev = toe[2]

                elevs.append(toe_elev)
                
                
            vals = np.stack([elevs,pc.BTBf[ii]])
            vals_all.append(vals)
       
    fig,ax = plt.subplots(round(len(vals_all)/2),2)
    axx = []
    for i in ax:
        for ii in range(0,2):
            axx.append(i[ii])
    for s in range(0,len(vals_all)):
        axx[s].plot(vals_all[s][0],vals_all[s][1],'.')
        # ax[s].axis('equal')
    axx[len(vals_all)-1].set_xlabel('Toe elev (m)')
    axx[len(vals_all)-1].set_ylabel('$B_T/B_f$')
    
    vals_all_avg = []
    for i in vals_all:
        vals_all_avg.append([np.nanmean(i[0]),np.nanmean(i[1])])
        
    col = cm.get_cmap('cividis', 6).colors
    fig,ax = plt.subplots(1)
    for i in range(0,len(vals_all_avg)):
        if i==0 or i==1:
            color = col[0,:]
        elif i==2 or i ==3:
            color = col[1,:]
        elif i==4 or i==5 or i==6:
            color=col[2,:]
        elif i==7 or i==8:
            color=col[3,:]
        elif i==9:
            color = col[4,:]
    
        ax.plot(vals_all_avg[i][0],vals_all_avg[i][1],'s',markerfacecolor=color,markeredgecolor=color)
    ax.set_xlabel('Avg. dune toe elev (m)')
    ax.set_ylabel('Avg. $B_T/B_f$')
        
                
 
def examineAndPlot_BTVsBf():  
    
    from sklearn.linear_model import LinearRegression
    
    vals_all = np.empty([0,2])
    for d in range(0,len(dates)):
        
        file = [iii for iii in files if str(dates[d])+'-' in iii and 'manual_transects' in iii][0]
    
    
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        pc = pickle.load(f)
        
        
        for ii in range(0,len(pc.T_scarps)):
            Bf = pc.Bf[ii]
            BT = pc.BT[ii]
            vals = np.transpose(np.vstack([Bf,BT]))
            vals_all = np.vstack([vals_all,vals])
            
            
    fig,ax = plt.subplots(1)
    ax.plot(vals_all[:,0],vals_all[:,1],'k.')
    ax.set_xlabel('Bf')
    ax.set_ylabel('BT')
    
    reg = LinearRegression().fit(vals_all[:,0][~np.isnan(vals_all[:,0])].reshape(-1,1),vals_all[:,1][~np.isnan(vals_all[:,0])].reshape(-1,1))
    yhat = reg.predict(vals_all[:,0][~np.isnan(vals_all[:,0])].reshape(-1,1))
    sst = np.sum((vals_all[:,1][~np.isnan(vals_all[:,0])]-np.mean(vals_all[:,1][~np.isnan(vals_all[:,0])]))**2)
    ssr = np.sum((yhat-np.mean(vals_all[:,1][~np.isnan(vals_all[:,0])]))**2)
    r2 = ssr/sst
    

def examineAndPlot_BTVsBf_ByOtherVars(): 
    '''
       Color scatter  by:
        - Mean water level
        - TWL relative to dune toe elev
        - Beach width
        - pre-storm volume
        - pre-storm dune slope
        - storm duration
        - IG energy
    '''
    
    from sklearn.linear_model import LinearRegression
    
    vals_all = np.empty([0,2])
    results = pd.DataFrame(columns=['BT','Bf','TWLRelativeToToe','Beach Width','Impact Duration','Beach Volume','Dune Volume'])
    for d in range(0,len(dates)):
        
        file = [iii for iii in files if str(dates[d])+'-' in iii and 'manual_transects' in iii][0]
    
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        pc = pickle.load(f)
        
        T_scarps = pc.T_scarps
        scarpToes = pc.scarpToes
        Bf = pc.Bf
        BT = pc.BT
        
        bdate = np.int64(file[0:8])
        edate = np.int64(file[9:17])
        hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=44056,buoyDepth=17.8)
        
        
        for scarp in range(0,len(T_scarps)):
            for ii in range(0,len(T_scarps[scarp][0])):
                
                BT_here = BT[scarp][ii]
                Bf_here = Bf[scarp][ii]
                
                # Calc mean water level during storm(s) #
                
                # Calc TWL relative to toe elev #
                twl = hydro.calcTWL(Bf[scarp][ii],exceedance=2)
                twl_relative_dtoe = max(twl[:,1])-scarpToes[scarp][0][ii][2]
                
                # Calc pre-storm beach width (toe to shoreline)
                x_toe = scarpToes[scarp][0][ii][0]
                try:
                    x_sl = utils.transectElevIntercept(0.36,T_scarps[scarp][0][ii]['X'],T_scarps[scarp][0][ii]['Z'])
                except IndexError: # Error when transects does not reach MHW #
                    x_sl = np.nan
                width = x_sl-x_toe
                
                # Calc impact duration #
                impact_dur_i = np.where(twl[:,1]>scarpToes[scarp][0][ii][2])[0]
                impact_dur = abs((len(impact_dur_i)-1)/2) # Assuming sampling is half hour #
                
                # Calc pre-storm volume (beach) #
                z_high = scarpToes[scarp][0][ii][2]
                z_low=0.75
                x = T_scarps[scarp][0][ii]['X']
                z = T_scarps[scarp][0][ii]['Z']
                x_vol = x[np.logical_and(z<=z_high,z>=z_low)]
                z_vol = z[np.logical_and(z<=z_high,z>=z_low)]
                vol_beach = np.trapz(z_vol,x_vol)
                
                # Calc pre-storm volume (dune) #
                z_high = max(T_scarps[scarp][0][ii]['Z'])
                z_low=scarpToes[scarp][0][ii][2]
                x = T_scarps[scarp][0][ii]['X']
                z = T_scarps[scarp][0][ii]['Z']
                x_vol = x[np.logical_and(z<=z_high,z>=z_low)]
                z_vol = z[np.logical_and(z<=z_high,z>=z_low)]
                vol_dune = np.trapz(z_vol,x_vol)
                
                results = results.append({'BT':BT_here,
                                              'Bf':Bf_here,
                                              'TWLRelativeToToe':twl_relative_dtoe,
                                              'Beach Width':width,
                                              'Impact Duration':impact_dur,
                                              'Beach Volume':vol_beach,
                                              'Dune Volume':vol_dune},
                                             ignore_index=True)
                
        for colorVal in ['TWLRelativeToToe','Beach Width','Impact Duration','Beach Volume','Dune Volume'] :
            fig = plt.figure()
            ax = plt.axes([0.15,0.15,0.8,0.8])
            h = ax.scatter(results['Bf'],results['BT'],5,results[colorVal])
            fig.colorbar(h,label=colorVal)
            ax.set_xlabel('Bf')
            ax.set_ylabel('BT')
            ax.plot(ax.get_xlim(),(0,0),'k',linewidth=1)
            plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/BfVsBT_by'+str(colorVal.replace(' ',''))+'.png',dpi=350)

        
        
        
    
def examineAndPlot_WLAndPostStormScarpPos():
    
    from sklearn.linear_model import LinearRegression
    
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if 'scarpResults' in i])
    dates_b = np.unique([int(i[0:8]) for i in files])
    dates_e = np.unique([int(i[9:17]) for i in files])

    res_exc = []
    for exc in [2,7,16,23,30,50]:
        res_date = []
        for i in range(0,len(dates_b)):
            bdate = dates_b[i]
            edate = dates_e[i]
            hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=44056,buoyDepth=17.8)
            
            file = [iii for iii in files if str(dates_b[i])+'-' in iii and 'manual_transects' in iii][0]
            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
            scarpResults=pickle.load(f)
            Bf = scarpResults.Bf
            scarpToes = scarpResults.scarpToes
            
            for scarp in range(0,len(Bf)):
                res_scarp = np.empty([0,2])
                Bf_scarp = Bf[scarp]
                scarpToes_scarp = scarpToes[scarp][1]
                for s,t in zip(Bf_scarp,scarpToes_scarp):
                    twl = hydro.calcTWL(s,exceedance=exc)
                    twl_max = max(twl[:,1])
                    
                    res_scarp = np.vstack([res_scarp,np.hstack([twl_max,t[2]])])
                res_date.append(res_scarp)
        res_exc.append(res_date)
    
    fig,ax = plt.subplots(3,2)
    axx = []
    for i in ax:
        for ii in range(0,2):
            axx.append(i[ii])
    for s in range(0,6):
        allv = np.empty([0,2])
        meanv = np.empty([0,2])
        for scarp in range(0,10):
            axx[s].plot(res_exc[s][scarp][:,0],res_exc[s][scarp][:,1],'.')
            allv = np.vstack([allv,res_exc[s][scarp]])
            meanv = np.vstack([allv,np.mean(res_exc[s][scarp],axis=0)])
        
        axx[s].set_title('R'+str([2,7,16,23,30,50][s])+'%')
        axx[s].set_xlabel('Max TWL (m)')
        axx[s].set_ylabel('Post-storm toe')
        
        reg = LinearRegression().fit(allv[:,0][~np.isnan(allv[:,0])].reshape(-1,1),allv[:,1][~np.isnan(allv[:,0])].reshape(-1,1))
        yhat = reg.predict(allv[:,0][~np.isnan(allv[:,0])].reshape(-1,1))
        sst = np.sum((allv[:,1][~np.isnan(allv[:,0])]-np.mean(allv[:,1][~np.isnan(allv[:,0])]))**2)
        ssr = np.sum((yhat-np.mean(allv[:,1][~np.isnan(allv[:,0])]))**2)
        r2 = ssr/sst
        axx[s].text(min(axx[s].get_xlim())+.2,min(axx[s].get_ylim()),'$R^2=$'+str(round(r2,2)))
        
    fig,ax = plt.subplots(3,2)
    axx = []
    for i in ax:
        for ii in range(0,2):
            axx.append(i[ii])
    for s in range(0,6):
        meanv = np.empty([0,2])
        for scarp in range(0,10):
            meanv = np.vstack([meanv,np.nanmean(res_exc[s][scarp],axis=0)])
        axx[s].plot(meanv[:,0],meanv[:,1],'s')
        axx[s].set_title('R'+str([2,7,16,23,30,50][s])+'%')
        axx[s].set_xlabel('avg Max TWL (m)')
        axx[s].set_ylabel('avg Post-storm toe')
        
        reg = LinearRegression().fit(meanv[:,0][~np.isnan(meanv[:,0])].reshape(-1,1),meanv[:,1][~np.isnan(meanv[:,0])].reshape(-1,1))
        yhat = reg.predict(meanv[:,0][~np.isnan(meanv[:,0])].reshape(-1,1))
        sst = np.sum((meanv[:,1][~np.isnan(meanv[:,0])]-np.mean(meanv[:,1][~np.isnan(meanv[:,0])]))**2)
        ssr = np.sum((yhat-np.mean(meanv[:,1][~np.isnan(meanv[:,0])]))**2)
        r2 = ssr/sst
        axx[s].text(min(axx[s].get_xlim())+.2,min(axx[s].get_ylim()),'$R^2=$'+str(round(r2,2)))
         
    

def examineAndPlot_WLAndBTQuadrants():
    
    from matplotlib import cm 
    from matplotlib.colors import ListedColormap,LinearSegmentedColormap
    from sklearn.linear_model import LinearRegression
    from scipy import stats

        
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if 'scarpResults' in i])
    dates_b = np.unique([int(i[0:8]) for i in files])
    dates_e = np.unique([int(i[9:17]) for i in files])

    res_exc = []
    for exc in [2,7,16,23,30,50]:
        res_date = []
        for i in range(0,len(dates_b)):
            bdate = dates_b[i]
            edate = dates_e[i]
            if int(str(bdate)[0:4])<2020:
                ws = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/'+str(bdate)[0:4]+'_17m.txt'
            else:
                ws = 44056
            hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=17.8)

            
            file = [iii for iii in files if str(dates_b[i])+'-' in iii and '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in iii]
            if len(file)>0:
                file = file[0]
                f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
                scarpResults=pickle.load(f)
                Bf = scarpResults.Bf
                BT = scarpResults.BT
                scarpToes = scarpResults.scarpToes
                T_scarps = scarpResults.T_scarps
                
                for scarp in range(0,len(Bf)):
                    res_scarp = np.empty([0,3])
                    BT_scarp = BT[scarp]
                    Bf_scarp = Bf[scarp]
                    scarpToes_scarp = scarpToes[scarp][0]
                    scarpToes_scarp_aft = scarpToes[scarp][1]
                    T_scarps_scarp_bef = T_scarps[scarp][0]
                    T_scarps_scarp_aft = T_scarps[scarp][1]
                    for s,t,toe,toe_aft,tran_bef,tran_aft in zip(Bf_scarp,BT_scarp,scarpToes_scarp,scarpToes_scarp_aft,
                                                                 T_scarps_scarp_bef,T_scarps_scarp_aft):
                        twl = hydro.calcTWL(s,exceedance=exc)
                        twl_max = max(twl[:,1])
                        
                        # fig,ax = plt.subplots(1)
                        # ax.plot(tran_bef['X'],tran_bef['Z'],'k')
                        # ax.plot(tran_aft['X'],tran_aft['Z'],'grey')
                        # ax.plot(toe[0],toe[2],'k')
                        # ax.plot(toe_aft[0],toe_aft[2])
                        # ax.plot(ax.get_xlim(),(twl_max,twl_max),'c')
                        # fig.show()
                        # plt.pause(2)
                        # plt.close('all')
                            
                        res_scarp = np.vstack([res_scarp,np.hstack([twl_max-toe[2],t,s])])
                    res_date.append(res_scarp)
                else:
                    pass
        if len(res_date)>0:
            res_exc.append(res_date)
    
    brg = cm.get_cmap('jet', len(res_exc[0]))
    cols = brg(range(len(res_exc[0])))
    
    fig,ax = plt.subplots(3,2)
    axx = []
    for i in ax:
        for ii in range(0,2):
            axx.append(i[ii])
    for s in range(0,6):
        allv = np.empty([0,3])
        meanv = np.empty([0,3])
        for scarp in range(0,len(res_exc[s])):
            axx[s].plot(res_exc[s][scarp][:,0],res_exc[s][scarp][:,1],'.',markersize=1.5,color=cols[scarp,:])
            allv = np.vstack([allv,res_exc[s][scarp]])
            meanv = np.vstack([allv,np.mean(res_exc[s][scarp],axis=0)])
        
        axx[s].set_title('R'+str([2,7,16,23,30,50][s])+'%')
        axx[s].set_xlabel('Max TWL - toe elev (m)')
        axx[s].set_xticks([-1,0,1,2])
        axx[s].set_xticklabels([-1,0,1,2])
        axx[s].set_ylabel('BT')
        axx[s].set_xlim(-2,2)
        axx[s].set_ylim(-3,1)
        # axx[s].set_ylim(-2,2)
        # axx[s].set_xlim(-2,2)
        # axx[s].axis('equal')
        
        # coefs,yhat,r2,p_values = utils.linearRegression(allv[:,0][~np.isnan(allv[:,0])],allv[:,1][~np.isnan(allv[:,0])])
        
        # axx[s].text(max(axx[s].get_xlim())-(np.diff(axx[s].get_xlim())/10),-1.7,'$R^2=$'+str(round(r2,2)))
        # axx[s].plot(allv[:,0][~np.isnan(allv[:,0])],yhat,'k--',linewidth=1)
        axx[s].plot((0,0),axx[s].get_ylim(),'k',linewidth=1)
        axx[s].plot(axx[s].get_xlim(),(0,0),'k',linewidth=1) 
        
        
        
        
        
        
       
        
 

        
def examineAndPlot_WLAndBTQuadrants_ByOtherVars():
    
    import matplotlib
    import matplotlib.cm as cm
        
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in i])
    dates_b = np.unique([int(i[0:8]) for i in files])
    dates_e = np.unique([int(i[9:17]) for i in files])


    results = pd.DataFrame(columns=['BT','Bf','Pre-storm toe elev','Post-storm toe elev','Toe dx','Toe dz',
                                    'Pre-storm scarp slope',
                                    'F','Beach Width','Impact Duration',
                                    '$\Delta V_B$','$\Delta V_D$','Dune Height'])
    for i in range(0,len(dates_b)):
        
        bdate = dates_b[i]
        edate = dates_e[i]
        ws = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/'+str(bdate)[0:4]+'_17m.txt'
        hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=17.8)
        
        file = [iii for iii in files if str(dates_b[i])+'-' in iii and '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in iii]
        if len(file)>0:
            file = file[0]
            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
            scarpResults=pickle.load(f)
            Bf = scarpResults.Bf
            BT = scarpResults.BT
            scarpToes = scarpResults.scarpToes
            T_scarps = scarpResults.T_scarps
            
            for scarp in range(0,len(T_scarps)):
                for ii in range(0,len(T_scarps[scarp][0])):
                    
                    BT_here = BT[scarp][ii]
                    Bf_here = Bf[scarp][ii]
                    
                    # Calc toe elev pre and post #
                    toe_pre = scarpToes[scarp][0][ii][2]
                    toe_post = scarpToes[scarp][1][ii][2]
                    
                    # Calc TWL freeboard #
                    twl = hydro.calcTWL(Bf_here,exceedance=2)
                    F = max(twl[:,1])-toe_pre
                    
                    # Calc toe dx #
                    dx = scarpToes[scarp][1][ii][0] - scarpToes[scarp][0][ii][0]
                    
                    # Calc toe dz #
                    dz = scarpToes[scarp][1][ii][2] - scarpToes[scarp][0][ii][2]
                    
                    # Calc pre-storm slope from toe to +0.5 m z #
                    t = T_scarps[scarp][0][ii]
                    # fig,ax = plt.subplots(1)
                    # ax.plot(t['X'],t['Z'])
                    t = t[np.logical_and(t['X']<=scarpToes[scarp][0][ii][0],t['Z']-scarpToes[scarp][0][ii][2]<=0.5,t['Z']-scarpToes[scarp][0][ii][2]>=0)]
                    m = np.mean(np.diff(t['Z'])/np.diff(t['X']))
                    # inter,m = utils.linearRegression(t['X'],t['Z'])
                    slope = abs(m)
                    # ax.plot(t['X'],t['Z'],'c')
                    # ax.plot(t['X'],(t['X']*m)+inter,'r')
                    # ax.set_title(str(m))
                    # fig.show()
                    # plt.pause(2)
                    # plt.close('all')
                    
                    # Calc pre-storm beach width (toe to shoreline)
                    x_toe = scarpToes[scarp][0][ii][0]
                    try:
                        x_sl = utils.transectElevIntercept(0.36,T_scarps[scarp][0][ii]['X'],T_scarps[scarp][0][ii]['Z'])
                    except IndexError: # Error when transects does not reach MHW #
                        x_sl = np.nan
                    width = x_sl-x_toe
                    
                    # Calc impact duration #
                    impact_dur_i = np.where(twl[:,1]>scarpToes[scarp][0][ii][2])[0]
                    if len(impact_dur_i)==0:
                        impact_dur=0
                    elif len(impact_dur_i)==1:
                        impact_dur=0.5
                    else:
                        impact_dur = (len(impact_dur_i)-1)/2 # Assuming sampling is half hour #
                    
                    # Calc pre- and post-storm volume (beach) #
                    vol_beach = []
                    for d in range(0,2):
                        z_high = scarpToes[scarp][d][ii][2]
                        z_low=0.5
                        x = T_scarps[scarp][d][ii]['X']
                        z = T_scarps[scarp][d][ii]['Z']
                        x_vol = x[np.logical_and(z<=z_high,z>=z_low)]
                        z_vol = z[np.logical_and(z<=z_high,z>=z_low)]
                        vol_beach.append(np.trapz(z_vol,x_vol))
                    vol_beach = vol_beach[1]-vol_beach[0]
                    
                    # Calc pre- and post-storm volume (dune) #
                    vol_dune1 = []
                    for d in range(0,2):
                        z_high = 6#max(T_scarps[scarp][d][ii]['Z'])
                        z_low=scarpToes[scarp][d][ii][2]
                        x_vol = x[np.logical_and(z<=z_high,z>=z_low)]
                        z_vol = z[np.logical_and(z<=z_high,z>=z_low)]
                        vol_dune1.append(np.trapz(z_vol,x_vol,dx=0.25))
                    vol_dune = vol_dune1[1]-vol_dune1[0]
                    
                    # Calc pre-storm dune height #
                    dHeight = max(T_scarps[scarp][0][ii]['Z'])-toe_pre
                    
             
                    results = results.append({'BT':BT_here,
                                              'Bf':Bf_here,
                                              'Pre-storm toe elev':toe_pre,
                                              'Post-storm toe elev':toe_post,
                                              'Toe dx':dx,
                                              'Toe dz':dz,
                                              'Pre-storm scarp slope':slope,
                                              'F':F,
                                              'Beach Width':width,
                                              'Impact Duration':impact_dur,
                                              '$\Delta V_B$':vol_beach,
                                              '$\Delta V_D$':vol_dune,
                                              'Dune Height':dHeight},
                                             ignore_index=True)
        else:
            pass

    # for colorVal,clims in zip(['Bf','Pre-storm scarp slope','Beach Volume','Toe dx','F','Impact Duration'],
    #                           [(0.05,0.12),(2,4),(2,4),(-10,0),(-1,1),(0,1),(10,50),(0,5),(0,60),(75,200)]):
    #     fig,ax = plt.subplots(3,2)
    #     axx = []
    #     for i in ax:
    #         for ii in range(0,2):
    #             axx.append(i[ii])
    #     for s in range(0,6):
    #         allv = np.empty([0,2])
    #         h = axx[s].scatter(res_exc[s]['TWLRelativeToToe'],res_exc[s]['BT'],2.5,res_exc[s][colorVal],vmin=clims[0],vmax=clims[1])
    #         if s==5:
    #             fig.colorbar(h,label=colorVal)
            
    #         allv = np.vstack([allv,np.transpose(np.vstack([res_exc[s]['TWLRelativeToToe'],res_exc[s]['BT']]))])
            
    #         axx[s].set_title('R'+str([2,7,16,23,30,50][s])+'%')
    #         axx[s].set_xlabel('Max TWL - toe elev (m)')
    #         axx[s].set_ylabel('BT')
    #         # axx[s].set_xlim(-2,1)
    #         axx[s].set_ylim(-2,1)
    #         axx[s].axis('equal')
            
    #         coefs,yhat,r2,p_values = utils.linearRegression(allv[:,0][~np.isnan(allv[:,0])],allv[:,1][~np.isnan(allv[:,0])])
            
    #         axx[s].text(max(axx[s].get_xlim())-1,min(axx[s].get_ylim())+0.1,'$R^2=$'+str(round(r2,2)))
    #         axx[s].plot(allv[:,0][~np.isnan(allv[:,0])],yhat,'k--',linewidth=1)
    #         axx[s].plot((0,0),axx[s].get_ylim(),'k',linewidth=1)
    #         axx[s].plot(axx[s].get_xlim(),(0,0),'k',linewidth=1)      
    
    # # Colored scatter plot with x=freeboard,y=BT, and color=var
    # focusExc = 2
    # fig,ax = plt.subplots(3,4)
    # axx = []
    # for i in ax:
    #     for ii in range(0,4):
    #         axx.append(i[ii])
    # axx[-1].axis('off')
        
    # for colorVal,clims,s in zip(['Bf','Pre-storm toe elev','Post-storm toe elev','Toe dx','Toe dz','Pre-storm scarp slope','Beach Width','Impact Duration','Beach Volume','Dune Volume'],
    #                           [(0.05,0.12),(2,4),(2,4),(-10,0),(-1,1),(0,1),(10,50),(0,5),(0,60),(75,200)],
    #                           range(0,12)):
        
    #     h = axx[s].scatter(res_exc[focusExc]['TWLRelativeToToe'],res_exc[focusExc]['BT'],2.5,res_exc[focusExc][colorVal],vmin=clims[0],vmax=clims[1])
    #     fig.colorbar(h,ax=axx[s])
            
            
    #     axx[s].set_title(colorVal)
    #     axx[s].set_xlabel('Max TWL - toe elev (m)')
    #     axx[s].set_ylabel('BT')
    #     # axx[s].set_xlim(-2,1)
    #     axx[s].set_ylim(-2,1)
    #     axx[s].axis('equal')
    #     axx[s].plot((0,0),axx[s].get_ylim(),'k',linewidth=1)
    #     axx[s].plot(axx[s].get_xlim(),(0,0),'k',linewidth=1) 
        
    # Binned bar plot of BT vs var #
    from matplotlib.patches import Rectangle
    binsize = 0.1
    binedges = np.arange(-.6,.5,binsize)
    iBins = np.digitize(np.array(results['BT']),binedges)
    iPos = np.where(np.array(results['BT'])>0)[0]
    iNeg = np.where(np.array(results['BT'])<0)[0]
    
    fig,ax = plt.subplots(2,3,figsize=(8,4))
    axx = []
    for i in ax:
        for ii in range(0,3):
            axx.append(i[ii])
        
    for colorVal,s in zip(['Bf','Pre-storm scarp slope','Dune Height','Toe dx','F','Impact Duration'],range(0,6)):
        
        val_pos_mean = np.mean(results[colorVal][iPos])
        val_neg_mean = np.mean(results[colorVal][iNeg])
               
        for binn in range(1,len(binedges)):
            val_mean = np.mean(results[colorVal][np.where(iBins==binn)[0]])
            val_std = np.std(results[colorVal][np.where(iBins==binn)[0]])
            r = Rectangle((0,binedges[binn-1]),val_mean,binedges[binn]-binedges[binn-1],
                          edgecolor='k',facecolor='grey')
            axx[s].add_patch(r)
            axx[s].plot((val_mean,val_mean+val_std),(binedges[binn-1]+((binedges[binn]-binedges[binn-1])/2),binedges[binn-1]+((binedges[binn]-binedges[binn-1])/2)),'k')
            axx[s].plot((val_mean,val_mean-val_std),(binedges[binn-1]+((binedges[binn]-binedges[binn-1])/2),binedges[binn-1]+((binedges[binn]-binedges[binn-1])/2)),'k')
            axx[s].set_ylim(min(binedges),max(binedges))
            axx[s].set_xlabel(colorVal)
            if colorVal=='Bf':
                axx[s].text(0.033,((binedges[binn-1]-binedges[binn])/2)+binedges[binn],'n='+str(len(np.where(iBins==binn)[0])),va='center',ha='center')
            if s==0 or s==3:
                axx[s].set_ylabel('BT')
                
        m1 = axx[s].plot(val_pos_mean,min(binedges),'b*',alpha=1)
        m2 = axx[s].plot(val_neg_mean,min(binedges),'r*',alpha=1)
        m1[0].set_clip_on(False)
        m2[0].set_clip_on(False)
        fig.legend([m1[0],m2[0]],['Mean for (+)BT','Mean for (-)BT'],loc='upper center',ncol=2,handletextpad=0.1)            
    
    
    
    # # Uncolored scatter plot of variable vs. BT #
    fig,ax = plt.subplots(3,4)
    axx = []
    for i in ax:
        for ii in range(0,4):
            axx.append(i[ii])
    axx[-1].axis('off')
        
    for colorVal,s in zip(['Bf','Pre-storm scarp slope','Dune Height','Toe dx','F','Impact Duration'],range(0,6)):
                            
        
        axx[s].plot(results[colorVal],results['BT'],'.')
      
            
            
        axx[s].set_xlabel(colorVal)
        axx[s].set_ylabel('BT')
        # axx[s].set_xlim(-2,1)
        axx[s].set_ylim(-2,1)
        # axx[s].axis('equal')
            
   
            
def examineAndPlot_WLAndBTQuadrants_ByOtherVars_ForEachStorm():   
    
    import matplotlib
    import matplotlib.cm as cm
        
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in i])
    dates_b = [20130306, 20170918, 20170922, 20180303, 20190904, 20190910, 20191011, 20191114, 20200910, 20210318]
    dates_e = [20130320, 20170922, 20170929, 20180309, 20190910, 20190924, 20191015, 20191119, 20200925, 20210323]
    names = ['2013NE','Jose','Maria','Riley','Dorian','Humberto','2019NE1','2019NE2','Teddy','2021NE']

    results_byStorm = pd.DataFrame(columns=['Name','results'])
    for i in range(0,len(dates_b)):
        
        results = pd.DataFrame(columns=['BT','Bf','Pre-storm toe elev','Post-storm toe elev','Toe dx','Toe dz',
                                'Pre-storm scarp slope',
                                'F','Beach Width','Impact Duration',
                                '$\Delta V_B$','$\Delta V_D$','Dune Height'])
        bdate = dates_b[i]
        edate = dates_e[i]
        ws = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/'+str(bdate)[0:4]+'_17m.txt'
        hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=17.8)

        
        file = [iii for iii in files if str(dates_b[i])+'-' in iii and '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in iii]
        if len(file)>0:
            file = file[0]
            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
            scarpResults=pickle.load(f)
            Bf = scarpResults.Bf
            BT = scarpResults.BT
            scarpToes = scarpResults.scarpToes
            T_scarps = scarpResults.T_scarps
            
            for scarp in range(0,len(T_scarps)):
                for ii in range(0,len(T_scarps[scarp][0])):
                    
                    BT_here = BT[scarp][ii]
                    Bf_here = Bf[scarp][ii]
                    
                    # Calc toe elev pre and post #
                    toe_pre = scarpToes[scarp][0][ii][2]
                    toe_post = scarpToes[scarp][1][ii][2]
                    
                    # Calc TWL freeboard #
                    twl = hydro.calcTWL(Bf_here,exceedance=2)
                    F = max(twl[:,1])-toe_pre
                    
                    # Calc toe dx #
                    dx = scarpToes[scarp][1][ii][0] - scarpToes[scarp][0][ii][0]
                    
                    # Calc toe dz #
                    dz = scarpToes[scarp][1][ii][2] - scarpToes[scarp][0][ii][2]
                    
                    # Calc pre-storm slope from toe to +0.5 m z #
                    t = T_scarps[scarp][0][ii]
                    # fig,ax = plt.subplots(1)
                    # ax.plot(t['X'],t['Z'])
                    t = t[np.logical_and(t['X']<=scarpToes[scarp][0][ii][0],t['Z']-scarpToes[scarp][0][ii][2]<=0.5,t['Z']-scarpToes[scarp][0][ii][2]>=0)]
                    m = np.mean(np.diff(t['Z'])/np.diff(t['X']))
                    # inter,m = utils.linearRegression(t['X'],t['Z'])
                    slope = abs(m)
                    # ax.plot(t['X'],t['Z'],'c')
                    # ax.plot(t['X'],(t['X']*m)+inter,'r')
                    # ax.set_title(str(m))
                    # fig.show()
                    # plt.pause(2)
                    # plt.close('all')
                    
                    # Calc pre-storm beach width (toe to shoreline)
                    x_toe = scarpToes[scarp][0][ii][0]
                    try:
                        x_sl = utils.transectElevIntercept(0.36,T_scarps[scarp][0][ii]['X'],T_scarps[scarp][0][ii]['Z'])
                    except IndexError: # Error when transects does not reach MHW #
                        x_sl = np.nan
                    width = x_sl-x_toe
                    
                    # Calc impact duration #
                    impact_dur_i = np.where(twl[:,1]>scarpToes[scarp][0][ii][2])[0]
                    if len(impact_dur_i)==0:
                        impact_dur=0
                    elif len(impact_dur_i)==1:
                        impact_dur=0.5
                    else:
                        impact_dur = (len(impact_dur_i)-1)/2 # Assuming sampling is half hour #
                    
                    # Calc pre- and post-storm volume (beach) #
                    vol_beach = []
                    for d in range(0,2):
                        z_high = scarpToes[scarp][d][ii][2]
                        z_low=0.5
                        x = T_scarps[scarp][d][ii]['X']
                        z = T_scarps[scarp][d][ii]['Z']
                        x_vol = x[np.logical_and(z<=z_high,z>=z_low)]
                        z_vol = z[np.logical_and(z<=z_high,z>=z_low)]
                        vol_beach.append(np.trapz(z_vol,x_vol))
                    vol_beach = vol_beach[1]-vol_beach[0]
                    
                    # Calc pre- and post-storm volume (dune) #
                    vol_dune1 = []
                    for d in range(0,2):
                        z_high = 6#max(T_scarps[scarp][d][ii]['Z'])
                        z_low=scarpToes[scarp][d][ii][2]
                        x_vol = x[np.logical_and(z<=z_high,z>=z_low)]
                        z_vol = z[np.logical_and(z<=z_high,z>=z_low)]
                        vol_dune1.append(np.trapz(z_vol,x_vol,dx=0.25))
                    vol_dune = vol_dune1[1]-vol_dune1[0]
                    
                    # fig,ax = plt.subplots(1)
                    # ax.plot(T_scarps[scarp][0][ii]['X'],T_scarps[scarp][0][ii]['Z'],'k')
                    # ax.plot(scarpToes[scarp][0][ii][0],scarpToes[scarp][0][ii][2],'k.')
                    # ax.plot(T_scarps[scarp][0][ii]['X'][np.where(abs(T_scarps[scarp][0][ii]['Z']-6) == min(abs(T_scarps[scarp][0][ii]['Z']-6)))],
                    #            T_scarps[scarp][0][ii]['Z'][np.where(abs(T_scarps[scarp][0][ii]['Z']-6) == min(abs(T_scarps[scarp][0][ii]['Z']-6)))],'r.')
                    # ax.plot(T_scarps[scarp][1][ii]['X'],T_scarps[scarp][1][ii]['Z'],'grey')
                    # ax.plot(scarpToes[scarp][1][ii][0],scarpToes[scarp][1][ii][2],'.',color='grey')
                    # ax.plot(T_scarps[scarp][1][ii]['X'][np.where(abs(T_scarps[scarp][1][ii]['Z']-6) == min(abs(T_scarps[scarp][1][ii]['Z']-6)))],
                    #            T_scarps[scarp][1][ii]['Z'][np.where(abs(T_scarps[scarp][1][ii]['Z']-6) == min(abs(T_scarps[scarp][1][ii]['Z']-6)))],'b.')
                    # fig.suptitle(vol_dune)
                    # plt.pause(5)
                    # plt.close('all')
                    
                    # Calc pre-storm dune height #
                    dHeight = max(T_scarps[scarp][0][ii]['Z'])-toe_pre
                    
             
                    results = results.append({'BT':BT_here,
                                              'Bf':Bf_here,
                                              'Pre-storm toe elev':toe_pre,
                                              'Post-storm toe elev':toe_post,
                                              'Toe dx':dx,
                                              'Toe dz':dz,
                                              'Pre-storm scarp slope':slope,
                                              'F':F,
                                              'Beach Width':width,
                                              'Impact Duration':impact_dur,
                                              '$\Delta V_B$':vol_beach,
                                              '$\Delta V_D$':vol_dune,
                                              'Dune Height':dHeight},
                                             ignore_index=True)
                    
            name = names[i]
            results_byStorm = results_byStorm.append({'Name':name,'results':results},ignore_index=True)
        else:
            pass



    fig,ax = plt.subplots(2,3,figsize=(8,4))
    axx = []
    for i in ax:
        for ii in range(0,3):
            axx.append(i[ii])
          
    for colorVal,s in zip(['Bf','Pre-storm scarp slope','Dune Height','Toe dx','F','Impact Duration'],range(0,6)):
        h = []
        for ii in range(0,len(results_byStorm)):
            name = results_byStorm.loc[ii]['Name']
            results = results_byStorm.loc[ii]['results']
            BT_mean = np.nanmean(results['BT'])
            BT_std = np.std(results['BT'])
            BT_25 = np.percentile(results['BT'],25)
            BT_75 = np.percentile(results['BT'],75)
            val_mean = np.nanmean(results[colorVal])
            val_std = np.std(results[colorVal])
            val_25 = np.percentile(results[colorVal],25)
            val_75 = np.percentile(results[colorVal],75)
            
            n = len(results[colorVal])
            alpha_scaled = ((n-20)/(345-20))
            if alpha_scaled<0.1:
                alpha_scaled=0.1
            # print(alpha_scaled)
            
            
            
            axx[s].plot( (val_mean,val_25),(BT_mean,BT_mean),'k',linewidth=1 )
            axx[s].plot( (val_mean,val_75),(BT_mean,BT_mean),'k',linewidth=1 )
            axx[s].plot( (val_mean,val_mean),(BT_mean,BT_25),'k',linewidth=1 )
            axx[s].plot( (val_mean,val_mean),(BT_mean,BT_75),'k',linewidth=1 )
            hh = axx[s].plot(val_mean,BT_mean,'s',alpha=alpha_scaled,markersize=14)
            h.append(hh[0])
            axx[s].set_xlabel(colorVal)
            if s==0 or s==3:
                axx[s].set_ylabel('BT')
    fig.legend(h,['2013NE','Jose','Riley','Dorian','2019NE1','2019NE2','Teddy','2021NE'],loc=7,frameon=False,handletextpad=0.1)
            
        
     

           


        
     
        
     
        

def examineAndPlot_EnvParamsInScarpedVsNonScarpedRegions():  
    
    from scipy.stats import ks_2samp
    
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if 'scarpResults' in i])
    dates_b = np.unique([int(i[0:8]) for i in files])
    dates_e = np.unique([int(i[9:17]) for i in files]) 
    
    for i in range(0,len(dates_b)):
        bdate = dates_b[i]
        edate = dates_e[i]
        if int(str(bdate)[0:4])==2018:
            ws = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/2018_17m.txt'
        elif int(str(bdate)[0:4])==2019:
            ws = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/2019_17m.txt'
        else:
            ws = 44056
        hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=17.8)
        
        # Load all the transects #
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+str(bdate)+'-'+str(edate)+'_5km_T.pkl','rb')
        T = pickle.load(f)
       
        # Load the identified scarp transects and other scarp vars #
        file = [iii for iii in files if str(dates_b[i])+'-' in iii and 'manual_transects' in iii][0]
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        scarpResults=pickle.load(f)
        Bf = scarpResults.Bf
        BT = scarpResults.BT
        scarpToes = scarpResults.scarpToes
        T_scarps = scarpResults.T_scarps
        
        # Create a boolean array where locations in T that are also in T_scarps are False #
        y_scarp = []
        for scarp in T_scarps:
            y_scarp_i = [scarp[0][i]['Y'][0] for i in range(0,len(scarp[0]))]
            y_scarp.append(y_scarp_i)
        # y_scarp = [item for sublist in y_scarp for item in sublist]
        
        y_all = [T[0][i]['Y'][0] for i in range(0,len(T[0]))]
        # booa = [i not in y_scarp for i in y_all] # True everywhere the transect isn't a scarp transect #
        booa11 = []
        for yy in y_scarp:
            booa1 = [i<np.nanmin(yy) or i>np.nanmax(yy) for i in y_all]
            booa11.append(booa1)
        booa = []
        for ii in range(0,len(y_all)):
            if not booa11[0][ii] or not booa11[1][ii]:
                booa.append(False)
            else:
                booa.append(True)
        
        
        # Calculate the dune toe locations #
        def indexListWithBool(nList,nBool):
            keepList = []
            for i in range(0,len(nList)):
                if nBool[i]:
                    keepList.append(nList[i])
            return keepList
        T_dune = [indexListWithBool(T[0],booa),indexListWithBool(T[1],booa)]
        
        lenList = [len(t) for t in T_dune[0]]
        T_dune = [indexListWithBool(T_dune[0],np.array(lenList)>1),indexListWithBool(T_dune[1],np.array(lenList)>1)]
        
        scarpLab = claris.scarpManager()
        
        cats = ['Bf','Pre-storm toe elev','Pre-storm scarp slope','TWLRelativeToToe','Beach Width','Impact Duration',
                    'Beach Volume','Dune Volume']
        labels=['$B_f$','$Z_{toe}$','$B_d$','$\Delta_{twl}$','$X_b$','$t_i$','$V_b$','$V_d$']
        fig,ax = plt.subplots(len(cats),1,sharex=True,figsize=(3,6.5))
        fig.subplots_adjust(left=0.245,bottom=0.058,right=0.9,top=0.898,wspace=0.2,hspace=0.71)
        ax[-1].set_xlim(.5,1.5)
        ax[-1].set_xticks([1])
        ax[-1].set_xticklabels(['contour'])
        fig.suptitle(str(bdate)+'-'+str(edate))
        c=0
        p_all = []
        for m in ['contour']:
            c+=1
            
            duneToes = scarpLab.idDuneToe(T_dune,method=m)
            
            # Plot the toe locations for non-scarped vs scarped regions #
            pre = []
            post = []
            for d in range(0,len(scarpToes)):
                pre1 = scarpToes[d][0]
                post1 = scarpToes[d][1]
                pre.append(pre1)
                post.append(post1)
            pre = [item for sublist in pre for item in sublist]
            post = [item for sublist in post for item in sublist]
            
            # fig,ax = plt.subplots(1,2,sharey=True,sharex=True)
            # ax[0].plot(np.array(pre)[:,0],np.array(pre)[:,1],'r.')
            # ax[0].plot(duneToes[0][:,0],duneToes[0][:,1],'b.')
            # ax[0].set_xlabel('FRF X (m)')
            # ax[0].set_ylabel('FRF Y (m)')
            # h1 = ax[1].plot(np.array(post)[:,0],np.array(post)[:,1],'r.')
            # h2 = ax[1].plot(duneToes[1][:,0],duneToes[1][:,1],'b.')
            # ax[1].set_xlabel('FRF X (m)')
            # ax[1].legend([h1[0],h2[0]],['Scarp','Not scarp'])
            
            
            # Calculate the parameters for scarped regions #
            results_scarp = pd.DataFrame(columns=['Bf','Pre-storm toe elev',
                                        'Pre-storm scarp slope',
                                        'TWLRelativeToToe','Beach Width','Impact Duration',
                                        'Beach Volume','Dune Volume'])
            for scarp in range(0,len(T_scarps)):
                for ii in range(0,len(T_scarps[scarp][0])):            
                    
                    Bf_here = Bf[scarp][ii]
                    
                    # Calc toe elev pre and post #
                    toe_pre = scarpToes[scarp][0][ii][2]
                    
                    # Calc pre-storm slope from toe to +0.5 m z #
                    t = T_scarps[scarp][0][ii]
                    t = t[np.logical_and(t['Z']-scarpToes[scarp][0][ii][2]<=0.5,t['Z']-scarpToes[scarp][0][ii][2]>=0)]
                    coefs,yhat,r2,p_values = utils.linearRegression(t['X'],t['Z'])
                    m = coefs[1]
                    slope = abs(m)
                    
                    # Calc TWL relative to toe elev #
                    twl = hydro.calcTWL(Bf[scarp][ii],exceedance=2)
                    twl_relative_dtoe = max(twl[:,1])-scarpToes[scarp][0][ii][2]
                    
                    # Calc pre-storm beach width (toe to shoreline)
                    x_toe = scarpToes[scarp][0][ii][0]
                    try:
                        x_sl = utils.transectElevIntercept(0.36,T_scarps[scarp][0][ii]['X'],T_scarps[scarp][0][ii]['Z'])
                    except IndexError: # Error when transects does not reach MHW #
                        x_sl = np.nan
                    width = x_sl-x_toe
                    
                    # Calc impact duration #
                    impact_dur_i = np.where(twl[:,1]>scarpToes[scarp][0][ii][2])[0]
                    if len(impact_dur_i)==0:
                        impact_dur=0
                    elif len(impact_dur_i)==1:
                        impact_dur=0.5
                    else:
                        impact_dur = (len(impact_dur_i)-1)/2 # Assuming sampling is half hour #
                    
                    # Calc pre-storm volume (beach) #
                    z_high = scarpToes[scarp][0][ii][2]
                    z_low=0.75
                    x = T_scarps[scarp][0][ii]['X']
                    z = T_scarps[scarp][0][ii]['Z']
                    x_vol = x[np.logical_and(z<=z_high,z>=z_low)]
                    z_vol = z[np.logical_and(z<=z_high,z>=z_low)]
                    vol_beach = np.trapz(z_vol,x_vol)
                    
                    # Calc pre-storm volume (dune) #
                    z_high = max(T_scarps[scarp][0][ii]['Z'])
                    z_low=scarpToes[scarp][0][ii][2]
                    x = T_scarps[scarp][0][ii]['X']
                    z = T_scarps[scarp][0][ii]['Z']
                    x_vol = x[np.logical_and(z<=z_high,z>=z_low)]
                    z_vol = z[np.logical_and(z<=z_high,z>=z_low)]
                    vol_dune = np.trapz(z_vol,x_vol)
                    
                    results_scarp = results_scarp.append({'Bf':Bf_here,
                                                  'Pre-storm toe elev':toe_pre,
                                                  'Pre-storm scarp slope':slope,
                                                  'TWLRelativeToToe':twl_relative_dtoe,
                                                  'Beach Width':width,
                                                  'Impact Duration':impact_dur,
                                                  'Beach Volume':vol_beach,
                                                  'Dune Volume':vol_dune},
                                                 ignore_index=True)
                
    
            # Calculate the parameters for non-scarped regions #
            results_dune = pd.DataFrame(columns=['Bf','Pre-storm toe elev',
                                        'Pre-storm scarp slope',
                                        'TWLRelativeToToe','Beach Width','Impact Duration',
                                        'Beach Volume','Dune Volume'])
            for ii in range(0,len(T_dune[0])):            
                
                x_beta = T_dune[0][ii]['X'][T_dune[0][ii]['X']>=duneToes[0][ii][0]]
                z_beta = T_dune[0][ii]['Z'][T_dune[0][ii]['X']>=duneToes[0][ii][0]]
                try:
                    m = abs((z_beta[-1]-z_beta[0])/(x_beta[-1]-x_beta[0]))
                except:
                    m = np.nan
                Bf_here = m
                
                # Calc toe elev pre and post #
                toe_pre = duneToes[0][ii][2]
                
                # Calc pre-storm slope from toe to +0.5 m z #
                t = T_dune[0][ii]
                t = t[np.logical_and(t['Z']-toe_pre<=0.5,t['Z']-toe_pre>=0)]
                coefs,yhat,r2,p_values = utils.linearRegression(t['X'],t['Z'])
                m = coefs[1]
                slope = abs(m)
                
                # Calc TWL relative to toe elev #
                twl = hydro.calcTWL(Bf_here,exceedance=2)
                twl_relative_dtoe = max(twl[:,1])-toe_pre
                
                # Calc pre-storm beach width (toe to shoreline)
                x_toe = duneToes[0][ii][0]
                try:
                    x_sl = utils.transectElevIntercept(0.36,T_dune[0][ii]['X'],T_dune[0][ii]['Z'])
                except IndexError: # Error when transects does not reach MHW #
                    x_sl = np.nan
                width = x_sl-x_toe
                
                # Calc impact duration #
                impact_dur_i = np.where(twl[:,1]>toe_pre)[0]
                if len(impact_dur_i)==0:
                    impact_dur=0
                elif len(impact_dur_i)==1:
                    impact_dur=0.5
                else:
                    impact_dur = (len(impact_dur_i)-1)/2 # Assuming sampling is half hour #
                
                # Calc pre-storm volume (beach) #
                z_high = toe_pre
                z_low=0.75
                x = T_dune[0][ii]['X']
                z = T_dune[0][ii]['Z']
                x_vol = x[np.logical_and(z<=z_high,z>=z_low)]
                z_vol = z[np.logical_and(z<=z_high,z>=z_low)]
                vol_beach = np.trapz(z_vol,x_vol)
                
                # Calc pre-storm volume (dune) #
                z_high = max(T_dune[0][ii]['Z'])
                z_low=toe_pre
                x = T_dune[0][ii]['X']
                z = T_dune[0][ii]['Z']
                x_vol = x[np.logical_and(z<=z_high,z>=z_low)]
                z_vol = z[np.logical_and(z<=z_high,z>=z_low)]
                vol_dune = np.trapz(z_vol,x_vol)
                
                results_dune = results_dune.append({'Bf':Bf_here,
                                              'Pre-storm toe elev':toe_pre,
                                              'Pre-storm scarp slope':slope,
                                              'TWLRelativeToToe':twl_relative_dtoe,
                                              'Beach Width':width,
                                              'Impact Duration':impact_dur,
                                              'Beach Volume':vol_beach,
                                              'Dune Volume':vol_dune},
                                             ignore_index=True)   
                
            
            p_cat = []
            for ii in range(0,len(cats)):
                
                # ks_stat,p = ks_2samp(results_dune[cats[ii]][~np.isnan(results_dune[cats[ii]])],
                #                      results_scarp[cats[ii]][~np.isnan(results_scarp[cats[ii]])],
                #                      alternative='two-sided',mode='exact')
                
                
                H,X1 = np.histogram( results_dune[cats[ii]][~np.isnan(results_dune[cats[ii]])], bins = 100, density = True )
                dx = X1[1] - X1[0]
                ecdf1_data = np.cumsum(H)*dx
                H,X1 = np.histogram( results_scarp[cats[ii]][~np.isnan(results_scarp[cats[ii]])], bins = 100, density = True )
                dx = X1[1] - X1[0]
                ecdf2_data = np.cumsum(H)*dx
                stat_data = np.max(np.abs(ecdf1_data-ecdf2_data))
                
                import random
                N = len(results_dune[cats[ii]][~np.isnan(results_dune[cats[ii]])])+len(results_scarp[cats[ii]][~np.isnan(results_scarp[cats[ii]])])
                allVals1 = [results_dune[cats[ii]][~np.isnan(results_dune[cats[ii]])],results_scarp[cats[ii]][~np.isnan(results_scarp[cats[ii]])]]
                allVals = pd.concat(allVals1,ignore_index=True)
                iTotal = np.arange(0,N)
                stat_all = []
                for rep in range(0,1000):
                    i1 = random.sample(list(iTotal), len(results_dune[cats[ii]][~np.isnan(results_dune[cats[ii]])]))
                    i2 = [i for i in iTotal if i not in i1]
                    
                    vals1 = allVals[i1]
                    vals2 = allVals[i2]
                    
                    H,X1 = np.histogram( vals1, bins = 100, density = True )
                    dx = X1[1] - X1[0]
                    ecdf1 = np.cumsum(H)*dx
                    H,X1 = np.histogram( vals2, bins = 100, density = True )
                    dx = X1[1] - X1[0]
                    ecdf2= np.cumsum(H)*dx
                    
                   
                    stat = np.max(np.abs(ecdf1-ecdf2))
                    stat_all.append(stat)
                    
                p = len(np.where(stat_all>stat_data)[0])/len(stat_all)
                p_cat.append(p)
                
                
                flierprops = dict(marker='x', markerfacecolor='k', markersize=3)
                ax[ii].boxplot([results_dune[cats[ii]][~np.isnan(results_dune[cats[ii]])],
                               results_scarp[cats[ii]][~np.isnan(results_scarp[cats[ii]])]],
                              positions=[c-0.1,c+0.1],manage_ticks=False,flierprops=flierprops)
                if c==1:
                    ax[ii].text(0.2,np.mean(ax[ii].get_ylim()),labels[ii],rotation=90)
                # ax[ii].text(c,max(ax[ii].get_ylim())+((max(ax[ii].get_ylim())-min(ax[ii].get_ylim()))/6),str(round(p,3)))
            p_all.append(p_cat)
        for meth in range(0,1):
            for cat in range(0,8):
                ax[cat].text(meth+1,max(ax[cat].get_ylim())+((max(ax[cat].get_ylim())-min(ax[cat].get_ylim()))/6),str(round(p_all[meth][cat],2)))

    
    
def examineAndPlot_ParamsForEachStorm():
    
    import datetime
    import matplotlib
    import matplotlib.cm as cm
        
    # files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in i])
    # dates_b = np.unique([int(i[0:8]) for i in files])
    # dates_e = np.unique([int(i[9:17]) for i in files])
    
    dates_b = [20130306, 20170918, 20170922, 20180303, 20190904, 20190910, 20191011, 20191114, 20200910, 20210318]
    dates_e = [20130320, 20170922, 20170929, 20180309, 20190910, 20190924, 20191015, 20191119, 20200925, 20210323]
    
    # fig,ax = plt.subplots(4,1)
    results = pd.DataFrame(columns=['date','BT','$\sum E /s$','maxE','$I_1$','$I_2$','$I_3$'])
    for i in range(0,len(dates_b)):
        bdate = dates_b[i]
        edate = dates_e[i]
        ws = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/'+str(bdate)[0:4]+'_17m.txt'
        hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=17.8)
        
        # dn = []
        # for ii in range(0,len(hydro.waves)):
        #     d = hydro.waves.iloc[ii]
        #     dn.append(datetime(int(d['yr']),int(d['mo']),int(d['day']),int(d['hr']),int(d['mm'])))
            
        # dn_wl = []
        # for ii in range(0,len(hydro.wl)):
        #     d = hydro.wl['Time'][ii]
        #     dn_wl.append(datetime(int(d[0:4]),int(d[5:7]),int(d[8:10]),int(d[11:13]),int(d[14:16])))
            
        # ax[0].plot(dn,np.array(hydro.waves['wvht (m)']).astype(float))
        # ax[0].set_xticklabels([])
        # ax[0].set_ylabel('$H_s$')
        
        # ax[1].plot(dn,np.array(hydro.waves['DPD (sec)']).astype(float))
        # ax[1].set_xticklabels([])
        # ax[1].set_ylabel('$T$')        
        
        # ax[2].plot(dn_wl,hydro.wl['wl_obs'])
        # ax[2].set_xticklabels([])
        # ax[2].set_ylabel('$\eta$')
        
        file = [iii for iii in files if str(dates_b[i])+'-' in iii and '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in iii]
        if len(file)>0:
            file = file[0]
            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
            scarpResults=pickle.load(f)
            BT = scarpResults.BT
            BT_all = [item for sublist in BT for item in sublist]
            BT_val = np.mean(BT_all)
            Bf = scarpResults.Bf
            Bf_all = [item for sublist in Bf for item in sublist] 
            Bf_mean = np.nanmean(Bf_all)
            toes = scarpResults.scarpToes

        else:
            Bf_mean = 0.07
            BT_val = 0
            toes = None
            
            # ax[3].plot(dn_wl[int(len(dn_wl)/2)],Bf_mean,'s')
            # ax[3].set_ylabel('B')
            
        # Calculate sumE #
        E = (1/16)*1025*9.81*(np.array(hydro.waves['wvht (m)']).astype(float)[np.array(hydro.waves['wvht (m)']).astype(float)>2]**2)
        sumE = np.trapz(E,x=None,dx=30*60)
        sumE = sumE/(datetime.datetime(int(str(edate)[0:4]),int(str(edate)[4:6]),int(str(edate)[6:8]))-
                     datetime.datetime(int(str(bdate)[0:4]),int(str(bdate)[4:6]),int(str(bdate)[6:8]))).total_seconds()
        
        # Calculate max E #
        maxE = max(E)
        
        # Calculate TWL and impact durations #
        twl = hydro.calcTWL(beta=Bf_mean,exceedance=2)
        durations = []
        for level in [1,2,3]:
            iAbove = np.where(twl[:,1]>level)[0]                
            if len(iAbove)==0:
                impact_dur=0
            elif len(iAbove)==1:
                impact_dur=0.5
            else:
                impact_dur = (len(iAbove)-1)/2 # Assuming sampling is half hour #
            durations.append(impact_dur)
        dur1 = durations[0]
        dur2 = durations[1]
        dur3 = durations[2]
        
        # Calculate SWL exceedance above dToe everywhere #
        if toes:
            swl_dur = []
            BTBf = []
            for scarp in range(0,len(toes)):
                swl_dur_transect = []
                BTBf_transect = []
                for iii in range(0,len(toes[scarp][0])):
                    toe = toes[scarp][0][iii]
                    Bf_here = Bf[scarp][iii]
                    BT_here = BT[scarp][iii]
                    BTBf_here = BT_here/Bf_here
                    swl = hydro.calcTWL(beta=Bf_here,exceedance=50)
                    iAbove = np.where(swl[:,1]>toe[2])[0] 
                    if not np.isnan(Bf_here):
                        if len(iAbove)==0:
                            impact_dur=0
                        elif len(iAbove)==1:
                            impact_dur=0.5
                        else:
                            impact_dur = (len(iAbove)-1)/2 # Assuming sampling is half hour #
                    else:
                        impact_dur = np.nan
                    swl_dur_transect.append(impact_dur)
                    BTBf_transect.append(BTBf_here)
                [swl_dur.append(p) for p in swl_dur_transect]
                [BTBf.append(o) for o in BTBf_transect]
        else:
            swl_dur = []
            BTBf = []
                

        
        
        results = results.append({'date':bdate,
                                  'BT':BT_val,
                                  '$\sum E /s$':sumE,
                                  'maxE':maxE,
                                  '$I_1$':dur1,
                                  '$I_2$':dur2,
                                  '$I_3$':dur3,
                                  'MWL>dToe':swl_dur,
                                  'BTBf':BTBf},
                                   ignore_index=True)
       
            
    fig,axx = plt.subplots(1,5,figsize=(8,2))  
    names = ['2013NE','Jose','Maria','Riley','Dorian','Humberto','2019NE1','2019NE2','Teddy','2021NE']
    for colorVal,s in zip(['$\sum E /s$','maxE','$I_1$','$I_2$','$I_3$'],range(0,len(names))):
        for ii in range(0,len(names)):
            axx[s].plot(results[colorVal][ii],results['BT'][ii],'ws')
            if ii==2 or ii==5:
                axx[s].text(results[colorVal][ii],results['BT'][ii],names[ii],va='center',ha='center',fontsize=8,color='r')
            else:
                axx[s].text(results[colorVal][ii],results['BT'][ii],names[ii],va='center',ha='center',fontsize=8)
               
        # axx[s].plot(results[colorVal],results['BT'],'s')
        axx[s].set_xlabel(colorVal)
        if s==0:
            axx[s].set_ylabel('BT')
        else:
            axx[s].set_yticklabels([])

    
    
    
def examineAndPlot_AnteceedantMorphology():
    
    def cropNormalizeGrid(transects1,xi):
        
        import copy       
        
        # Crop, normalize, and grid transects #
        transects = []
        transects_c = []
        transects_c_n = []
        transects_c_n_i = []
        for t in transects1:
            if max(t['Z'])>3.5 and min(t['Z'])<1: # Make sure looking at full profiles #
                t_c = t[np.logical_and(t['Z']>=1,t['Z']<=3)] # Crop the transect #
                
                dx = np.diff(t_c['X'][~np.isnan(t_c['Z'])])
                if len(dx)>0:
                    if max(dx)<0.3: # Make sure there is not a bunch of missing data in the profile #
                        t_c_n = copy.copy(t_c) # Normalize the transect #
                        t_c_n['X'] = (t_c_n['X']-min(t_c_n['X']))/(max(t_c_n['X'])-min(t_c_n['X']))
                        
                        zi = np.interp(xi,t_c_n['X'],t_c_n['Z'])
                        t_c_n_i = np.zeros([len(xi)],dtype=[('X','<f8'),('Z','<f8')])
                        t_c_n_i['X'] = xi
                        t_c_n_i['Z'] = zi
                        
                        
                        transects.append(t)
                        transects_c.append(t_c)
                        transects_c_n.append(t_c_n)
                        transects_c_n_i.append(t_c_n_i)
 
            else:
                pass
            
            
            # for i in range(0,len(transects)):
            #     fig,ax = plt.subplots(1,2,sharey=True)
            #     ax[0].plot(transects[i]['X'],transects[i]['Z'],'k.')
            #     ax[0].plot(transects_c[i]['X'],transects_c[i]['Z'],'r.')
            #     ax[1].plot(transects_c_n[i]['X'],transects_c_n[i]['Z'],'k.')
            #     plt.pause(1)
            #     plt.close('all')
            
            
        return transects,transects_c,transects_c_n,transects_c_n_i
        
        
    # First, calc the average profile. This only needs to be done once. #       
    direcs = sorted(os.listdir('/Users/frfuser/Documents/pyCLARIS_project/data')) 
    direcs = direcs[1:-1]
    direcs = ['/Users/frfuser/Documents/pyCLARIS_project/data/'+i for i in direcs]
    analysisLen=5
    xi = np.arange(0,1.01,0.01)
    
    # transects_all = []
    # for direc in direcs:
    #     # Create the cropped point cloud if i doesn't exist #
    #     if not os.path.exists(direc+'/FRF_'+str(analysisLen)+'km.las'):
    #         if analysisLen == 1:
    #             croperU = 'frf'
    #         else:
    #             croperU = '5km'
    #         claris.createFRFLas(direc,croper=croperU)
    #         os.rename(direc+'/FRF.las',direc+'/FRF_'+str(analysisLen)+'km.las')
            
    #     # Pull transects #
    #     file = direc+'/FRF_'+str(analysisLen)+'km.las'
    #     pc = claris.pcManager(file)
    #     transects1 = pc.createTransects(dy=5)
        
    #     transects,transects_c,transects_c_n,transects_c_n_i = cropNormalizeGrid(transects1,xi)
        
    #     # Average all profiles for this survey #
    #     transect_mean = np.zeros([len(xi)],dtype=[('X','<f8'),('Z','<f8')])
    #     transect_mean['X'] = xi
    #     for x in range(0,len(xi)):
    #         vals = [i['Z'][x] for i in transects_c_n_i]
    #         transect_mean['Z'][x] = np.mean(vals)
            
    #     fig,ax = plt.subplots(1)
    #     for t in transects_c_n_i:
    #         ax.plot(t['X'],t['Z'],'k',linewidth=0.5)
    #     ax.plot(transect_mean['X'],transect_mean['Z'],'r',linewidth=2)
    #     fig.show()
                
    #     transects_all.append(transect_mean)  
        
    # transect_mean = np.zeros([len(xi)],dtype=[('X','<f8'),('Z','<f8')])
    # transect_mean['X'] = xi
    # for x in range(0,len(xi)):
    #     vals = [i['Z'][x] for i in transects_all]
    #     transect_mean['Z'][x] = np.mean(vals)
    # meanProf = transect_mean
   
    ########################################################################
    
    
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/meanProf.pkl','rb')
    meanProf = pickle.load(f)
    
    
    # # Calc deviation of profiles from every survey from mean profile #
    dates_b = [20130306, 20170918, 20170922, 20180303, 20190904, 20190910, 20191011, 20191114, 20200910, 20210318]
    # direcs_preStorm = []
    # for d in dates_b:
    #     direc = [i for i in direcs if str(d) in i][0]
    #     direcs_preStorm.append(direc)
    
    # results = pd.DataFrame(columns=['date','deltaV','deltaV_high','deltaV_mid','deltaV_low'])
    # s=-1
    # for direc in direcs_preStorm:
        
    #     s+=1
            
    #     # Pull transects #
    #     analysisLen = 5
    #     file = direc+'/FRF_'+str(analysisLen)+'km.las'
    #     pc = claris.pcManager(file)
    #     transects1 = pc.createTransects(dy=5)
        
    #     # Crop, normalize, and grid transects #
    #     t = cropNormalizeGrid(transects1,xi)
    #     transects,transects_c,transects_c_n,transects_c_n_i = cropNormalizeGrid(transects1,xi)
        
    #     if len(transects_c_n)>0:
    #         # Calc deviations of each profile from the mean #
    #         vol_mean = np.trapz(meanProf['Z'],meanProf['X'])
    #         iLow = meanProf['X']>=0.67
    #         iMid = np.logical_and(meanProf['X']>=0.33,meanProf['X']<0.67)
    #         iHigh = meanProf['X']<0.33
            
    #         deltaV = []
    #         deltaV_low = []
    #         deltaV_mid = []
    #         deltaV_high = []
    #         for t in transects_c_n_i:
    #             # Volume deviation #
    #             vol_prof = np.trapz(t['Z'],t['X'])
    #             deltaV.append(vol_prof-vol_mean)
                
    #             # Volume deviation in each third in the x-direction #
    #             vols_low = (np.trapz(meanProf['Z'][iLow],meanProf['X'][iLow]),np.trapz(t['Z'][iLow],t['X'][iLow]))
    #             vols_mid = (np.trapz(meanProf['Z'][iMid],meanProf['X'][iMid]),np.trapz(t['Z'][iMid],t['X'][iMid]))
    #             vols_high = (np.trapz(meanProf['Z'][iHigh],meanProf['X'][iHigh]),np.trapz(t['Z'][iHigh],t['X'][iHigh]))
    #             deltaV_low.append(np.diff(vols_low))
    #             deltaV_mid.append(np.diff(vols_mid))
    #             deltaV_high.append(np.diff(vols_high))
              
                
    #         fig = plt.figure(figsize=(8,4))
    #         ax = plt.axes([0.35,0.6,0.3,0.3])
    #         axs1 = plt.axes([0.05,0.1,0.3,0.3],sharex=ax,sharey=ax)
    #         axs2 = plt.axes([0.35,0.1,0.3,0.3],sharex=ax,sharey=ax)
    #         # axs2.set_yticklabels([])
    #         axs3 = plt.axes([0.65,0.1,0.3,0.3],sharex=ax,sharey=ax)
    #         # axs3.set_yticklabels([])
            
    #         ax.plot(deltaV)
    #         axs1.plot(deltaV_low)
    #         axs2.plot(deltaV_mid)
    #         axs3.plot(deltaV_high)
            
    #         results = results.append({'date':dates_b[s],'deltaV':deltaV,'deltaV_high':deltaV_high,'deltaV_mid':deltaV_mid,'deltaV_low':deltaV_low},
    #                             ignore_index=True)
    #     else:
    #         results = results.append({'date':dates_b[s],'deltaV':np.nan,'deltaV_high':np.nan,'deltaV_mid':np.nan,'deltaV_low':np.nan},
    #                             ignore_index=True)
            
        
    
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/meanProf_compResults.pkl','rb')
    results = pickle.load(f)
    
    names = ['2013NE','Jose','Maria','Riley','Dorian','Humberto','2019NE1','2019NE2','Teddy','2021NE']
    updown = ['up','up','na','up','down','na','up','down','up','up']
    
    
    fig = plt.figure(figsize=(8,4))
    ax = plt.axes([0.1,0.1,0.85,0.85])
    for i in range(0,len(dates_b)):
        
        low = np.array(results['deltaV_low'][results['date']==dates_b[i]])[0]
        mid = np.array(results['deltaV_mid'][results['date']==dates_b[i]])[0]
        high = np.array(results['deltaV_high'][results['date']==dates_b[i]])[0]
        
        if updown[i] == 'up':
            fc = 'b'
        elif updown[i] == 'down':
            fc = 'r'
        else:
            fc = 'k'
        rH = Rectangle((i-0.3,0),0.2,np.nanmean(high),edgecolor='k',facecolor=fc)
        rM = Rectangle((i-0.1,0),0.2,np.nanmean(mid),edgecolor='k',facecolor=fc)
        rL = Rectangle((i+0.1,0),0.2,np.nanmean(low),edgecolor='k',facecolor=fc)
        
        ax.add_artist(rH)
        ax.add_artist(rM)
        ax.add_artist(rL)
        
    ax.set_xlim(-1,len(dates_b))
    ax.set_ylim(-0.15,0.15)
    ax.plot((-1,len(dates_b)),(0,0),'k')
    ax.set_xticks(np.arange(0,len(dates_b)))
    ax.set_xticklabels(names)
    ax.set_ylabel('$\delta V$')
    
    # Look at example profiles compared to meanProf #
    # Teddy #
    direc = [i for i in direcs if str(20200910) in i][0] 
    file = direc+'/FRF_'+str(analysisLen)+'km.las'
    pc = claris.pcManager(file)
    transects1 = pc.createTransects(dy=5)
    t = cropNormalizeGrid(transects1,xi)
    transects,transects_c,transects_c_n,transects_c_n_i = cropNormalizeGrid(transects1,xi)
    for t in np.flip(transects_c_n):
        fig,ax = plt.subplots(1)
        ax.set_title('Teddy')        
        ax.plot(meanProf['X'],meanProf['Z'],'k',linewidth=2,zorder=5)
        ax.plot(t['X'],t['Z'],'b',zorder=5)
        ax.set_ylim(ax.get_ylim())
        ax.plot((0.67,0.67),ax.get_ylim(),'k',zorder=1)
        ax.plot((0.33,0.33),ax.get_ylim(),'k',zorder=1)
        plt.show()
        plt.pause(2)
        plt.close('all')
        
    # 2019NE1 #
    direc = [i for i in direcs if str(20191011) in i][0] 
    file = direc+'/FRF_'+str(analysisLen)+'km.las'
    pc = claris.pcManager(file)
    transects1 = pc.createTransects(dy=5)
    t = cropNormalizeGrid(transects1,xi)
    transects,transects_c,transects_c_n,transects_c_n_i = cropNormalizeGrid(transects1,xi)
    for t in transects_c_n:
        fig,ax = plt.subplots(1)
        ax.set_title('Teddy')        
        ax.plot(meanProf['X'],meanProf['Z'],'k',linewidth=2,zorder=5)
        ax.plot(t['X'],t['Z'],'b',zorder=5)
        ax.set_ylim(ax.get_ylim())
        ax.plot((0.67,0.67),ax.get_ylim(),'k',zorder=1)
        ax.plot((0.33,0.33),ax.get_ylim(),'k',zorder=1)
        plt.show()
        plt.pause(2)
        plt.close('all')
        
    # Jose #
    direc = [i for i in direcs if str(20170918) in i][0] 
    file = direc+'/FRF_'+str(analysisLen)+'km.las'
    pc = claris.pcManager(file)
    transects1 = pc.createTransects(dy=5)
    for t in transects1:
        fig,ax = plt.subplots(1)
        ax.set_title('Jose')        
        ax.plot(t['X'],t['Z'],'b',zorder=5) 
        plt.show()
        plt.pause(0.2)
        plt.close('all')
    t = cropNormalizeGrid(transects1,xi)
    transects,transects_c,transects_c_n,transects_c_n_i = cropNormalizeGrid(transects1,xi)
    for t in transects_c_n:
        fig,ax = plt.subplots(1)
        ax.set_title('Teddy')        
        ax.plot(meanProf['X'],meanProf['Z'],'k',linewidth=2,zorder=5)
        ax.plot(t['X'],t['Z'],'b',zorder=5)
        ax.set_ylim(ax.get_ylim())
        ax.plot((0.67,0.67),ax.get_ylim(),'k',zorder=1)
        ax.plot((0.33,0.33),ax.get_ylim(),'k',zorder=1)
        plt.show()
        plt.pause(2)
        plt.close('all')


    
    
    
   
    
    
