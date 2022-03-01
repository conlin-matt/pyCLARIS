# Standard library imports #
import os

# 3rd party inputs #
import matplotlib.pyplot as plt
import numpy as np
import pickle

# Project imports
import pyCLARIS.pyCLARISAnalysis as claris

iRow=-1
vals = np.empty([3,4],dtype=object)
for thresh_lc in [10,25,50]:
    iRow+=1
    iCol=-1
    for thresh_min in [0,1.5]:
        iCol+=1
        print('thresh_min = '+str(thresh_min))

        dataDict = {}
        for date in [ ['20210318','20210323'],['20200910','20200925'],['20200305','20200310'],['20191114','20191119'],
                        ['20191011','20191015'],['20190910','20190924'],['20190904','20190910'],['20180303','20180309'],
                        ['20170922','20170929'],['20170918','20170922'],['20130306','20130320']]:#,
                       # ['20110825','20110829'] ]:
            
            date_pre = date[0]#'20191114'
            date_post = date[1]#'20191119'
            analysisLen = 5
            direcs = ['/Users/frfuser/Documents/pyCLARIS_project/data/'+i for i in sorted(os.listdir('/Users/frfuser/Documents/pyCLARIS_project/data')) if date_pre in i or date_post in i]
            scarpToeMethod = 'ml'
            slopeMethod = 'ep'
            DEMGridSize = '0.5' # 0.5 or 1 #
            savedToes = None
            savedFile = None
            thresh_vc = 0.5
                
                                
            if analysisLen == 1:
                xx = np.arange(50,150,0.5)
                yy = np.arange(0,1000,0.5)
            elif analysisLen == 5:
                xx = np.arange(-200,150,0.5)
                yy = np.arange(-1000,4000,0.5)
            
            
            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_dsms.pkl','rb'); dsms = pickle.load(f)
            f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_T.pkl','rb'); T = pickle.load(f)
            try:
                f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_regions_agg_'+str(thresh_vc)+'m.pkl','rb');regions_agg = pickle.load(f)
                f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_regions_deg_'+str(thresh_vc)+'m.pkl','rb');regions_deg = pickle.load(f)
            except:
                regions_agg = None
                regions_deg = None
            scarpLab = claris.scarpManager(thresh_vertChange=thresh_vc,thresh_longshoreContinuity=thresh_lc,thresh_minimumElev=thresh_min,thresh_slope_after=0)
            try:
                BTBf = scarpLab.calcBTOverBf(xx,yy,dsms[0],dsms[1],T,scarpToeMethod,slopeMethod,
                                                 regions_agg=regions_agg,regions_deg=regions_deg,
                                                 file_pre=direcs[0]+'/FRF_'+str(analysisLen)+'km.las',
                                                 file_post=direcs[1]+'/FRF_'+str(analysisLen)+'km.las',
                                                 savedFile=savedFile,
                                                 clf='mixed_clf')
            except:
                pass
            dataDict[date[0]] = scarpLab
        vals[iRow,iCol] = dataDict
            
            
            
            
            
import datetime
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from matplotlib import patches


yy = np.array([datetime.datetime(2011,1,1)+datetime.timedelta(days=i) for i in range(0,(365*11)+20,100)])
xx1 = np.arange(-1000,4100,100)         
            
fig,axx = plt.subplots(3,2,figsize=(12,7))
for iRow in range(0,len(axx)):
    for iCol in range(0,len(axx)-1):
        
        ax = axx[iRow,iCol]
        
        dates_storms = [datetime.datetime(2011,8,27),datetime.datetime(2013,3,13),
        datetime.datetime(2017,9,20),datetime.datetime(2017,9,25),
        datetime.datetime(2018,3,6),datetime.datetime(2019,9,7),
        datetime.datetime(2019,9,17),datetime.datetime(2019,10,13),
        datetime.datetime(2019,11,16),datetime.datetime(2020,9,17,12),
        datetime.datetime(2021,3,20,12)]
        names_storms = ['Irene','NorEaster','Jose','Maria','Riley','Dorian','Humberto','NorEaster','NorEaster','Teddy','NorEaster']
        for i in range(0,len(names_storms)):
            ax.plot((min(xx1),max(xx1)),(dates_storms[i],dates_storms[i]),'k--',linewidth=0.5)
            
        val = vals[iRow,iCol]
        
        for key in val.keys():
            
            data = val[key]
            
            scarpToes = data.scarpToes
            BT = data.BT
            T_scarps = data.T_scarps
            numScarps = len(T_scarps)
            
            date_pre = datetime.datetime(int(key[0:4]),int(key[4:6]),int(key[6:8]))
            date_post = date_pre+datetime.timedelta(days=6)
            date_plot = date_pre+((date_post-date_pre)/2)
            
            for scarp in range(0,numScarps):
                try:
                    loc = np.arange(min(scarpToes[scarp][0][:,1]),max(scarpToes[scarp][0][:,1]))
                    BT_region = np.mean(BT[scarp])
                    xx = loc
                    yy = np.array([date_plot+datetime.timedelta(i) for i in range(-20,21)])
                    R = patches.Rectangle((xx[0],yy[0]),max(xx)-min(xx),max(yy)-min(yy),facecolor='k',edgecolor='k')
                    ax.add_artist(R)
                    # a = ax.pcolor(xx,yy,np.tile(BT_region,[len(yy),len(xx)]),facecolor='k',edgecolor='k')
                    
                except:
                    pass
        ax.set_xlim(-1000,4000)
axx[2,0].set_xlabel('FRF X (m)')
axx[2,1].set_xlabel('FRF X (m)')

        







           