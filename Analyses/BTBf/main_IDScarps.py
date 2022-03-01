# Standard library imports #
import os

# 3rd party imports #
import matplotlib.pyplot as plt
import numpy as np
import pickle

# Project imports
import pyCLARIS.pyCLARISAnalysis as claris


# # Uncomment to create your own pyBeach machine learning dune toe classification algorithm.
# # To do so, you will need to have IDd some scarps using the manual_transects method
# # so that your identified toes can be passed to the classifier creator. Once created,
# # you can change the clf arg in calcBTOverBf to 'custom_clf' #
# m = claris.scarpManager()
# m.create_pyBeachClassifier('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/')
###########################

for date in [ ['20210318','20210323'],['20200910','20200925'],['20191114','20191119'],
                ['20191011','20191015'],['20190910','20190924'],['20190904','20190910'],
                ['20180303','20180309'],['20170922','20170929'],['20130306','20130320'] ]:
    
    date_pre = date[0]
    date_post = date[1]
    analysisLen = 5
    direcs = ['/Users/frfuser/Documents/pyCLARIS_project/data/'+i for i in sorted(os.listdir('/Users/frfuser/Documents/pyCLARIS_project/data')) if date_pre in i or date_post in i]
    thresh_vertChangeU = 0.5
    scarpToeMethod = 'mc_supervised'
    slopeMethod = 'ep'
    DEMGridSize = '0.5' # 0.5 or 1 #
    
    if scarpToeMethod == 'manual_pc':
        if os.path.exists('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_scarpResults_manual.pkl'):
            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_scarpResults_manual.pkl','rb')
            savedFile = pickle.load(f)
            savedToes = savedFile.scarpToes
        else: 
            savedToes=None
    elif scarpToeMethod == 'manual_transects':
        if os.path.exists('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_scarpResults_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m_withoutSlopeCriteria.pkl'):
            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_scarpResults_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m_withoutSlopeCriteria.pkl','rb')
            savedFile = pickle.load(f)
            savedToes = savedFile.scarpToes
        else: 
            savedToes=None
    else:
        savedToes=None
        savedFile=None
        
        
        
    if analysisLen == 1:
        xx = np.arange(50,150,0.5)
        yy = np.arange(0,1000,0.5)
    elif analysisLen == 5:
        xx = np.arange(-200,150,0.5)
        yy = np.arange(0,4000,0.5)
    
    
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_dsms.pkl','rb'); dsms = pickle.load(f)
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_T.pkl','rb'); T = pickle.load(f)
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_regions_agg_'+str(thresh_vertChangeU)+'m.pkl','rb');regions_agg = pickle.load(f)
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_regions_deg_'+str(thresh_vertChangeU)+'m.pkl','rb');regions_deg = pickle.load(f)
    
    scarpLab = claris.scarpManager(thresh_vertChange=thresh_vertChangeU,thresh_longshoreContinuity=25,thresh_minimumElev=2,thresh_slope_after=0)
    BTBf = scarpLab.calcBTOverBf(xx,yy,dsms[0],dsms[1],T,scarpToeMethod,slopeMethod,
                                     regions_agg=regions_agg,regions_deg=regions_deg,
                                     file_pre=direcs[0]+'/FRF_'+str(analysisLen)+'km.las',
                                     file_post=direcs[1]+'/FRF_'+str(analysisLen)+'km.las',
                                     savedFile=savedFile,
                                     clf='mixed_clf')
    
    
    if len(scarpLab.BTBf)>0:
        with open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_scarpResults_'+scarpToeMethod+'.pkl','wb') as f:
            pickle.dump(scarpLab,f)
    
    
    
    
    
    # regions = scarpLab.scarpRegions
    # T_scarp = scarpLab.T_scarps
    # toes = scarpLab.scarpToes
    
    def plotExampleRegionForAGUPoster():
        regions = scarpLab.scarpRegions
        
        fig = plt.figure(figsize=(3.4,2.6))
        plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
        plt.rc('legend', fontsize=16) 
        ax = plt.axes([0.25,0.48,0.58,0.38])
        
        h = ax.pcolor(yy,xx,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')
        ax.invert_yaxis()  
        ax.plot(regions[8][:,1],regions[8][:,0],'k')
        ax.set_xlim(min(regions[8][:,1])-5,max(regions[8][:,1])+5)
        for tt in T[0]:
            if min(regions[8][:,1])<=tt['Y'][0]<=max(regions[8][:,1]):
                ax.plot(tt['Y'],tt['X'],'k--',linewidth=1)
        ax.set_xlabel('alongshore (m)')
        ax.set_ylabel('x-shore (m)')
        
        cbax = plt.axes([0.2,0.14,0.4,0.025])
        fig.colorbar(h,cbax,orientation='horizontal')
        fig.text(.62,.13,'$\Delta Z$ (m)',ha='left',va='center',fontsize=18)
        
        axsub = plt.axes([.83,.48,.15,.38],sharey=ax)
        axsub.plot(T[0][786]['Z'],T[0][786]['X'],'k',linewidth=2.5)
        axsub.yaxis.tick_right()
        axsub.xaxis.tick_top()
        
    
    # def plotIDdScarps():
    
    #     for i in range(0,len(regions)):
    
    #         fig = plt.figure(figsize=(6.5,6.5))
    #         plt.rcParams.update({'font.size': 8})
    
    #         ax1 = plt.axes([0.15,0.75,0.8,0.2])
    #         axt1 = plt.axes([0.15,0.46,0.3,0.18])
    #         axt2 = plt.axes([0.15,0.28,0.3,0.18],sharex=axt1)
    #         axt3 = plt.axes([0.15,0.1,0.3,0.18],sharex=axt1)
    
    #         axBox = plt.axes([0.5,0.05,0.45,0.6]);axBox.set_xticks([]);axBox.set_yticks([])
    #         axf1 = plt.axes([0.6,0.45,0.32,0.15])
    #         axf2 = plt.axes([0.6,0.3,0.32,0.15])
    #         axf3 = plt.axes([0.6,0.11,0.32,0.15])
    
    
    #         ax1.pcolor(yy,xx,np.transpose(dsms[1]-dsms[0]),vmin=-1,vmax=1,cmap='seismic_r')
    #         ax1.invert_yaxis()
    #         ax1.set_xlabel('FRF Y (m)')
    #         ax1.set_ylabel('FRF X (m)')
    #         ax1.set_title(date_pre+'-'+date_post+', method='+scarpToeMethod+', scarp '+str(i+1)+' of '+str(len(regions)),fontsize=8,pad=0.1,loc='left')
    #         ax1.set_xlim(min(regions[i][:,1])-100,max(regions[i][:,1])+100)
    #         ax1.set_ylim(max(regions[i][:,0])+20,min(regions[i][:,0])-20)
    #         ax1.plot(regions[i][:,1],regions[i][:,0],'m')
    
    #         if len(T_scarp[i][0])>0:
    #             i_T = np.round(np.linspace(0,len(T_scarp[i][0])-2,3)).astype(int)
    #             axt1.plot(T_scarp[i][0][i_T[0]]['X'],T_scarp[i][0][i_T[0]]['Z'],'k')
    #             axt1.plot(T_scarp[i][1][i_T[0]]['X'],T_scarp[i][1][i_T[0]]['Z'],'--',color='grey')
    #             axt1.plot(toes[i][0][i_T[0]][0],toes[i][0][i_T[0]][2],'k*')
    #             axt1.plot(toes[i][1][i_T[0]][0],toes[i][1][i_T[0]][2],'*',color='grey')
    #             ax1.plot(T_scarp[i][0][i_T[0]]['Y'],T_scarp[i][0][i_T[0]]['X'],'k',linewidth=0.5)
    #             axt2.plot(T_scarp[i][0][i_T[1]]['X'],T_scarp[i][0][i_T[1]]['Z'],'k')
    #             axt2.plot(T_scarp[i][1][i_T[1]]['X'],T_scarp[i][1][i_T[1]]['Z'],'--',color='grey')
    #             axt2.plot(toes[i][0][i_T[1]][0],toes[i][0][i_T[1]][2],'k*')
    #             axt2.plot(toes[i][1][i_T[1]][0],toes[i][1][i_T[1]][2],'*',color='grey')
    #             ax1.plot(T_scarp[i][0][i_T[1]]['Y'],T_scarp[i][0][i_T[1]]['X'],'k',linewidth=0.5)
    #             axt3.plot(T_scarp[i][0][i_T[2]]['X'],T_scarp[i][0][i_T[2]]['Z'],'k')
    #             axt3.plot(T_scarp[i][1][i_T[2]]['X'],T_scarp[i][1][i_T[2]]['Z'],'--',color='grey')
    #             axt3.plot(toes[i][0][i_T[2]][0],toes[i][0][i_T[2]][2],'k*')
    #             axt3.plot(toes[i][1][i_T[2]][0],toes[i][1][i_T[2]][2],'*',color='grey')            
    #             ax1.plot(T_scarp[i][0][i_T[2]]['Y'],T_scarp[i][0][i_T[2]]['X'],'k',linewidth=0.5)
    #             axt3.set_xlabel('FRF X (m)')
    
    
    #         axf1.plot(toes[i][0][:,1],scarpLab.BT[i],'-s')
    #         axf1.plot(toes[i][0][:,1],scarpLab.Bf[i],'-s')
    #         axf1.set_ylim(0,0.3)
    #         axf1.legend(['$B_T$','$B_f$'])
    #         axf1.set_xticklabels([])
            
    #         axf2.plot(toes[i][0][:,1],scarpLab.BTBf[i],'-s')
    #         axf2.set_ylim(-1,3)
    # ##        axf2.text(min(toes[i][0][:,1])+5,2.5,'$B_T/B_f$')
    # ##        axf2.text(min(toes[i][0][:,1])+5,-0.9,'FRF Y (m)')
    
    #         axf3.hist(scarpLab.BTBf[i],bins=np.arange(-1,3,0.2),density=True)
    #         axf3.set_xlabel('$B_T/B_f$')
    #         axf3.set_xticks([0,1,2])
    #         axf3.set_ylabel('Density')
    
    
    # #         plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/figs/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_scarp'+str(i+1)+'of'+str(len(regions))+'_'+scarpToeMethod+'.png',
    # #                     dpi=350)
    # #         plt.show()
    
    # plotIDdScarps()



