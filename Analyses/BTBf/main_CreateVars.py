# Standard library imports #
import os

# 3rd party inputs #
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.interpolate import griddata

# Project imports
import pyCLARIS.pyCLARISAnalysis as claris


for date in [ ['20210318','20210323'],['20200910','20200925'],['20200305','20200310'],['20191114','20191119'],
                ['20191011','20191015'],['20190910','20190924'],['20190904','20190910'],
                ['20180303','20180309'],['20170922','20170929'],['20130306','20130320'] ]:
    

        
    date_pre = date[0]
    date_post = date[1]
    savedVars = False
    analysisLen = 5
    direcs = ['/Users/frfuser/Documents/pyCLARIS_project/data/'+i for i in sorted(os.listdir('/Users/frfuser/Documents/pyCLARIS_project/data')) if date_pre in i or date_post in i]

            
    if analysisLen == 1:
        xx = np.arange(50,150,0.5)
        yy = np.arange(0,1000,0.5)
    elif analysisLen == 5:
        xx = np.arange(-200,150,0.5)
        yy = np.arange(0,4000,0.5)
        

    # Create or load in the DSMs and transects #
    if not savedVars:

        dsms = list()
        T = list()
        for direc in direcs:
            
            # Create FRF point cloud if it hasn't already been made
            if not os.path.exists(direc+'/FRF_'+str(analysisLen)+'km.las'):
                if analysisLen == 1:
                    croperU = 'frf'
                else:
                    croperU = '5km'
                claris.createFRFLas(direc,croper=croperU)
                os.rename(direc+'/FRF.las',direc+'/FRF_'+str(analysisLen)+'km.las')

            # Create the PC object #
            pc = claris.pcManager(direc+'/FRF_'+str(analysisLen)+'km.las') 

            # Grid the PC into a DSM and a colored image #
            print('Making DSM')
            dsm = pc.gridPC_Slocum(xx,yy,z_val='z',function='min')

            # Pull xshore transects from the PC #
            print('Making transects')
            transects = pc.createTransects(dy=5,y_min=min(yy),y_max=max(yy))
            # transects = pc.createTransects_fromDSM(xx,yy,dsm,dy=5)

            # Save the DSM and transects #
            dsms.append(dsm)
            T.append(transects)

        # Extract areas of change and plot on pre-storm DEM and RGB image #
        print('Extracting regions')
        regions_agg,regions_deg = claris.extractChangeAreas(xx,yy,dsms[0],dsms[1],thresh=0.5)
        
        for val in ['regions_agg','regions_deg']:
            with open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_'+val+'_0.5m.pkl','wb') as f:
                pickle.dump(eval(val),f)

        for val in ['dsms','T']:
            with open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_'+val+'.pkl','wb') as f:
                pickle.dump(eval(val),f)




    else:
        
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_dsms.pkl','rb'); dsms = pickle.load(f)
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_T.pkl','rb'); T = pickle.load(f)
        

    # def plotPointCloudAndGriddedForAGUPoster():
        # fig = plt.figure(figsize=(3.4,3))
        # plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
        # plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
        # plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
        # plt.rc('legend', fontsize=16) 
        # ax1 = plt.axes([0.25,0.62,0.6,0.35])
        # ax1.scatter(pcData[100000:150000]['Y'],pcData[100000:150000]['X'],.1,pcData[100000:150000]['Z'],vmin=2,vmax=3)
        # ax1.invert_yaxis()
        # ax1.set_xlim(21.97378467517421, 34.34252082231811)
        # ax1.set_ylim(91.74174349320019, 80.11063950014464)
        # ax1.set_xticklabels([])
        # ax1.text(22,80.5,'Point cloud',ha='left',va='top',fontweight='bold')
        
        # ax2 = plt.axes([0.25,0.25,0.6,0.35])
        # h = ax2.pcolor(yy,xx,np.transpose(dsm),vmin=2,vmax=3)
        # ax2.invert_yaxis()
        # ax2.set_xlim(21.97378467517421, 34.34252082231811)
        # ax2.set_ylim(91.74174349320019, 80.11063950014464)
        # ax2.plot((25,25),(91.74174349320019, 80.11063950014464),'k--')
        # ax2.plot((30,30),(91.74174349320019, 80.11063950014464),'k--')
        # ax2.text(22,80.5,'Gridded DSM + profiles',ha='left',va='top',fontweight='bold')
        # fig.text(0.05,0.6,'x-shore (m)',rotation=90,va='center',ha='center',fontsize=24)
        # fig.text(0.5,0.05,'alongshore (m)',va='center',ha='center',fontsize=24)
        
        # cbax = plt.axes([0.87,0.4,0.05,0.3])
        # fig.colorbar(h,cbax)
        # fig.text(0.899,0.55,'Elev (m)',rotation=90,va='center',ha='center',fontsize=14)
    


    # Plot the changes #
    def plotBothSurfacesAndChangeFig():
        fig = plt.figure(figsize=(6.5,4))
        ax = fig.subplots(3,sharex=True)
        plt.rcParams.update({'font.size': 8})
            
        count = -1
        titles = [date_pre+' (pre-storm)',date_post+' (post-storm)']
        for axis in ax[0:-1]:
            count+=1
            h = axis.pcolor(yy,xx,np.transpose(dsms[count]),vmin=0,vmax=8,cmap='gist_earth')
            axis.invert_yaxis()
            axis.set_title(titles[count],fontsize=8,pad=0.1)
            axis.set_xticks([])

        h2 = ax[2].pcolor(yy,xx,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')
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

        plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_bothSurfacesAndChangeFig.png',dpi=350)
##        plt.show()

    def changesOverviewFig():

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

        h = ax1.pcolor(yy,xx,np.transpose(dsms[0]),vmin=0,vmax=8,cmap='gist_earth')
        ax1.invert_yaxis()

        ax2.pcolor(yy,xx,np.transpose(dsms[1]),vmin=0,vmax=8,cmap='gist_earth')
        ax2.invert_yaxis()

        ax3.pcolor(yy,xx,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')    
        ax3.invert_yaxis()


        T_nums = [ [150,160,170,180],[350,360,370,380],[500,510,520,530],[650,660,670,680]]#,[800,810,820,830],[940,950,960,970] ]
        for i in range(0,4):
            for t in range(0,len(T_nums[i])):
                axx = eval('axt'+str(i+1)+'_'+str(t+1))
                axx.plot(T[0][T_nums[i][t]]['X'],T[0][T_nums[i][t]]['Z'],'k-')
                axx.plot(T[1][T_nums[i][t]]['X'],T[1][T_nums[i][t]]['Z'],'--',color='grey')
                axx.set_ylim(0,9)
                ax3.plot(T[0][T_nums[i][t]]['Y'],T[0][T_nums[i][t]]['X'],'k',linewidth=0.5) 

        plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_changesOverview.png',dpi=350)
##        plt.show()

    plotBothSurfacesAndChangeFig()
    changesOverviewFig()








   

### Examine results #
##regions = scarpLab.scarpRegions
##
##
##
##
##
##
### Save to catelog #
##f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/ScarpRetreatCatelog.pkl','rb'); catelog = pickle.load(f)
##data_dict = {'DSMs':dsms,'Transects':T,'Scarp Transects':scarpLab.T_scarps,'Bf':scarpLab.Bf,'BT':scarpLab.BT,'BT/Bf':scarpLab.BTBf}
##catelog[date_pre+'-'+date_post] = data_dict
##with open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/ScarpRetreatCatelog.pkl','wb') as f:
##    pickle.dump(catelog,f)




##
##
##def plotBTOverBf():
##
##    plt.rcParams.update({'font.size': 8})
##    fig = plt.figure(figsize=(3.5,4.5))
##    ax1 = plt.axes([.15,.7,.8,.25])
##    ax2 = plt.axes([.15,.4,.8,.25])
##    ax3 = plt.axes([.35,.1,.3,.2])
##    
##    ax1.plot(scarpLab.scarpToes[0][0][:,1],scarpLab.BT[0],'-s')
##    ax1.plot(scarpLab.scarpToes[0][0][:,1],scarpLab.Bf[0],'-s')
##    ax1.set_ylim(0,0.3)
##    ax1.legend(['$B_T$','$B_f$'])
##    ax1.set_xticklabels([])
####    ax1.set_title('Method: '+method,fontweight='normal',pad=0)
##    
##    ax2.plot(scarpLab.scarpToes[0][0][:,1],scarpLab.BTBf[0],'-s')
##    ax2.set_ylim(-1,3)
##    ax2.text(90,2.5,'$B_T/B_f$')
##    ax2.set_xlabel('FRF Y (m)')
##
##    ax3.hist(scarpLab.BTBf[0],bins=np.arange(-1,3,0.2),density=True)
##    ax3.set_xlabel('$B_T/B_f$')
##    ax3.set_xticks([0,1,2])
##    ax3.set_ylabel('Density')
##    plt.show()


