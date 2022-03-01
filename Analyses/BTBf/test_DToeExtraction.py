# Standard library imports #
import os

# 3rd party inputs #
import matplotlib.pyplot as plt
import numpy as np
import pickle
from pybeach.beach import Profile
from scipy.signal import detrend

# Project imports
import pyCLARIS.pyCLARISAnalysis as claris
import pyCLARIS.coastalGeoUtils as utils

dates = [ ['20210318','20210323'],['20200910','20200925'],['20200305','20200310'],['20191114','20191119'],
                ['20191011','20191015'],['20190910','20190924'],['20190904','20190910'],['20180303','20180309'],
                ['20170922','20170929'],['20170918','20170922'],['20130306','20130320'],
                ['20110825','20110829'] ]
    
date_pre = dates[1][0]
date_post = dates[1][1]
analysisLen = 5
direcs = ['/Users/frfuser/Documents/pyCLARIS_project/data/'+i for i in sorted(os.listdir('/Users/frfuser/Documents/pyCLARIS_project/data')) if date_pre in i or date_post in i]
thresh_vertChangeU = 0.5
scarpToeMethod = 'mc'
slopeMethod = 'ep'
DEMGridSize = '0.5' # 0.5 or 1 #
savedFile = None
savedToes = None
        
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

scarpLab = claris.scarpManager(thresh_vertChange=thresh_vertChangeU,thresh_longshoreContinuity=25,thresh_minimumElev=1.5,thresh_slope_after=0)
BTBf = scarpLab.calcBTOverBf(xx,yy,dsms[0],dsms[1],T,scarpToeMethod,slopeMethod,
                                 regions_agg=regions_agg,regions_deg=regions_deg,
                                 file_pre=direcs[0]+'/FRF_'+str(analysisLen)+'km.las',
                                 file_post=direcs[1]+'/FRF_'+str(analysisLen)+'km.las',
                                 savedFile=savedFile,
                                 clf='mixed_clf')
T_scarps = scarpLab.T_scarps



# Calculate BT values using both pd and mc on very limited profiles and compare results #
for pp in range(0,2):
    toes_pd = []
    toes_mc = []
    toes_ml = []
    toes_rr = []
    y = []
    x_pd = []
    x_mc = []
    x_ml = []
    x_rr = []
    z_pd = []
    z_mc = []
    z_ml = []
    z_rr = []
    T_all = []
    for scarp in range(0,len(T_scarps)):
        for t_num in range(0,len(T_scarps[scarp][pp])):
            
            T = T_scarps[scarp][pp][t_num]
            pb = Profile(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)],T['Z'][np.logical_and(T['Z']>=1,T['Z']<=5)])
          
            toe_pd = pb.predict_dunetoe_pd(dune_crest=None,shoreline=None)
            x_pd1 = T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_pd]
            z_pd1 = T['Z'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_pd]
            
            toe_ml = pb.predict_dunetoe_ml(clf_name='mixed_clf',dune_crest=None)[0]
            x_ml1 = T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_ml]
            z_ml1 = T['Z'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_ml]
            
            toe_rr = pb.predict_dunetoe_rr()
            x_rr1 = T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_rr]
            z_rr1 = T['Z'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_rr]
           
            # fig,ax = plt.subplots(1)
            # ax.plot(T['X'],T['Z'],'k',linewidth=3)
            # h1 = ax.plot(T['X'][np.logical_and(T['Z']>=2,T['Z']<=4)][toe_mc],T['Z'][np.logical_and(T['Z']>=2,T['Z']<=4)][toe_mc],'o')
            
            # plt.pause(1)
            # plt.show()
            # plt.close('all')
            
            # Manual mc method #
            # Max curvature #
            def smooth(x,window_len=7,window='hanning'):
                s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
                #print(len(s))
                if window == 'flat': #moving average
                    w=np.ones(window_len,'d')
                else:
                    w=eval('np.'+window+'(window_len)')
            
                y=np.convolve(w/w.sum(),s,mode='valid')
                return y[round(window_len/2-1):-round(window_len/2)]   
            
            xx = T['X']
            zz = T['Z']
            
            xx = xx[~np.isnan(zz)]
            zz = zz[~np.isnan(zz)]  
            
            xi = np.arange(min(xx),max(xx),1)
            zi = np.interp(xi,xx,zz)
            xx = xi
            zz = zi
            
            xx1 = xx#[np.logical_and(zz>=1,zz<=5)]
            zz1 = zz#[np.logical_and(zz>=1,zz<=5)]
            zz_smooth = zz1#smooth(zz1)
            iSub = np.logical_and(zz_smooth>=1,zz_smooth<=5)
            dx = np.gradient(xx1, xx1)  # first derivatives
            dy = np.gradient(zz_smooth, xx1)                
            d2x = np.gradient(dx, xx1)  # second derivatives
            d2y = np.gradient(dy, xx1)                
            curvature = d2y / ((1 + dy ** 2)) ** 1.5  # curvature     
            curvature_sub = curvature[iSub]
            toe_mc = np.where(curvature==max(curvature_sub))[0]  
            x_mc1 = xx1[toe_mc]
            z_mc1 = zz1[toe_mc]
            
            
            # JMull method #
            xx = T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)]
            zz = T['Z'][np.logical_and(T['Z']>=1,T['Z']<=5)]     
            p = np.polyfit(xx,zz,3)
            fit = p[0]*xx**3+p[1]*xx**2+p[2]*xx+p[3]
            fit_detrend = detrend(fit,axis=0)
            toe_mull = np.where(fit_detrend==min(fit_detrend))[0]
            
            
            # Brodie and Spore method #
            xx = xx[~np.isnan(zz)]
            zz = zz[~np.isnan(zz)]
            zz = smooth(zz)
            x_high = xx[np.where(zz==max(zz))[0][0]]
            z_high = zz[np.where(zz==max(zz))[0][0]]
            try:
                x_low = utils.transectElevIntercept(0.5,xx,zz)
                z_low = 0.5
            except IndexError:
                x_low = xx[-1]
                z_low = zz[-1]
            m = (z_high-z_low)/(x_high-x_low)
            b = z_high-(m*x_high)
            yhat = (m*xx[np.where(zz==z_high)[0][0]:np.where(zz==z_low)[0][0]])+b           
            dist = yhat-zz[np.where(zz==z_high)[0][0]:np.where(zz==z_low)[0][0]]
            firstGuess = np.where(dist==max(dist))[0]
            iSub = np.where(abs(xx-xx[firstGuess])<=2.5)[0]
          
            dx = np.gradient(xx, xx)  # first derivatives
            dy = np.gradient(zz, xx)
            d2x = np.gradient(dx, xx)  # second derivatives
            d2y = np.gradient(dy, xx)
            curvature = d2y / ((1 + dy ** 2)) ** 1.5  # curvature
            curvature_sub = curvature[iSub]
            iToe = np.where(curvature==max(curvature_sub))[0]

            # fig,ax = plt.subplots(2,1,sharex=True)
            # ax[0].plot(xx,zz)
            # ax[0].plot(xx[np.where(zz==z_high)[0][0]:np.where(zz==z_low)[0][0]],yhat,'r')
            # ax[0].plot(xx[firstGuess],zz[firstGuess],'r.')
            # ax[0].plot((xx[firstGuess]-5,xx[firstGuess]+5),(zz[firstGuess],zz[firstGuess]),'r--')
            # ax[1].plot(xx,curvature)
            # ax[1].plot(xx[iSub],curvature[iSub])
            # ax[1].plot(xx[iToe],curvature[iToe],'g.')
            # ax[0].plot(xx[iToe],zz[iToe],'g.')

            
            
            
            
            toes_pd.append(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_pd])
            toes_mc.append(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_mc])
            toes_ml.append(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_ml])
            toes_rr.append(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_rr])           
            y.append(T['Y'][0])
            x_pd.append(x_pd1)
            x_mc.append(x_mc1)
            x_ml.append(x_ml1)  
            x_rr.append(x_rr1)                       
            z_pd.append(z_pd1)
            z_mc.append(z_mc1)
            z_ml.append(z_ml1)
            z_rr.append(z_rr1)            
            T_all.append(T)
    
            
    if pp == 0:
        toes_pd_pre = toes_pd
        toes_mc_pre = toes_mc
        x_pd_pre = x_pd
        x_mc_pre = x_mc
        x_ml_pre = x_ml
        x_rr_pre = x_rr      
        z_pd_pre = z_pd
        z_mc_pre = z_mc
        z_ml_pre = z_ml
        z_rr_pre = z_rr        
        T_all_pre = T_all
    elif pp==1:
        toes_pd_post = toes_pd
        toes_mc_post = toes_mc
        x_pd_post = x_pd
        x_mc_post = x_mc
        x_ml_post = x_ml
        x_rr_post = x_rr        
        z_pd_post = z_pd
        z_mc_post = z_mc 
        z_ml_post = z_ml
        z_rr_post = z_rr        
        T_all_post = T_all


ex = 87
fig,ax = plt.subplots(1)
hh1 = ax.plot(T_all_pre[ex]['X'],T_all_pre[ex]['Z'],'k',linewidth=2)
hh2 = ax.plot(T_all_post[ex]['X'],T_all_post[ex]['Z'],'grey',linewidth=2)
manuall = plt.ginput(2)
h1 = ax.plot(x_pd_pre[ex],z_pd_pre[ex],'ro')
h2 = ax.plot(x_mc_pre[ex],z_mc_pre[ex],'bo')
h3 = ax.plot(x_ml_pre[ex],z_ml_pre[ex],'go')
h4 = ax.plot(x_rr_pre[ex],z_rr_pre[ex],'mo')
h5 = ax.plot(manuall[0][0],manuall[0][1],'ko')
ax.plot(x_pd_post[ex],z_pd_post[ex],'ro')
ax.plot(x_mc_post[ex],z_mc_post[ex],'bo')
ax.plot(x_ml_post[ex],z_ml_post[ex],'go')
ax.plot(x_rr_post[ex],z_rr_post[ex],'mo')
ax.plot(manuall[1][0],manuall[1][1],'ko')
ax.legend([h1[0],h2[0],h3[0],h4[0],h5[0]],['pd','mc','ml','rr','manual'],loc=1)


files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in i])       
manual_file = [iii for iii in files if date_pre in iii and date_post in iii and '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in iii][0]
f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+manual_file,'rb')
manual = pickle.load(f) 
BT_manual = manual.BT
BT_manual_all = []
for i in range(0,len(BT_manual)):
    for ii in range(0,len(BT_manual[i])):
        BT_manual_all.append(BT_manual[i][ii])



fig,ax = plt.subplots(1)
h1 = ax.hist((np.array(z_pd_post)-np.array(z_pd_pre))/-(np.array(x_pd_post)-np.array(x_pd_pre)),bins=list(np.linspace(-2,1,100))+[np.inf],density=True,cumulative=True,histtype='step',color='r',label='pd')
h2 = ax.hist((np.array(z_mc_post)-np.array(z_mc_pre))/-(np.array(x_mc_post)-np.array(x_mc_pre)),bins=list(np.linspace(-2,1,100))+[np.inf],density=True,cumulative=True,histtype='step',color='b',label='mc')
h3 = ax.hist(BT_manual_all,bins=list(np.linspace(-2,1,100))+[np.inf],density=True,cumulative=True,histtype='step',color='k',label='manual')
ax.legend(['pd','mc','manual'])
ax.plot((0,0),(0,1),'--',color='grey',linewidth=1)
ax.set_xlim(-1,1.5)
ax.set_xlabel('BT')




# Do the limited mc and pd comparison for all storms #
fig,ax = plt.subplots(2,4)
axx = []
for iRow in range(0,2):
    for iCol in range(0,4):
        axx.append(ax[iRow,iCol])

names = ['2021NE','Teddy','2019NE2','2019NE1','Dorian','Riley','Jose','2013NE']
ic = -1
for d in dates:
    
    date_pre = d[0]
    date_post = d[1]
    
    
    files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in i])       
    manual_file = [iii for iii in files if date_pre in iii and date_post in iii and '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in iii]

    if len(manual_file)>0:
        ic+=1
        manual_file = manual_file[0]
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+manual_file,'rb')
        manual = pickle.load(f) 
        BT_manual = manual.BT
        T_scarps = manual.T_scarps
        BT_manual_all = []
        for i in range(0,len(BT_manual)):
            for ii in range(0,len(BT_manual[i])):
                BT_manual_all.append(BT_manual[i][ii])
        
        
        
        for pp in range(0,2):
            toes_pd = []
            toes_mc = []
            toes_ml = []
            y = []
            x_pd = []
            x_mc = []
            x_ml = []
            z_pd = []
            z_mc = []
            z_ml = []
            T_all = []
            for scarp in range(0,len(T_scarps)):
                for t_num in range(0,len(T_scarps[scarp][pp])):
                    
                    T = T_scarps[scarp][pp][t_num]
                    pb = Profile(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)],T['Z'][np.logical_and(T['Z']>=1,T['Z']<=5)])
                  
                    toe_pd = pb.predict_dunetoe_pd(dune_crest=None,shoreline=None)
                    # toe_mc = pb.predict_dunetoe_mc()
                    toe_ml = pb.predict_dunetoe_ml(clf_name='mixed_clf',dune_crest=None)[0]
        
                   
                    # fig,ax = plt.subplots(1)
                    # ax.plot(T['X'],T['Z'],'k',linewidth=3)
                    # h1 = ax.plot(T['X'][np.logical_and(T['Z']>=2,T['Z']<=4)][toe_mc],T['Z'][np.logical_and(T['Z']>=2,T['Z']<=4)][toe_mc],'o')
                    
                    # plt.pause(1)
                    # plt.show()
                    # plt.close('all')
                    
                    xx = T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)]
                    zz = T['Z'][np.logical_and(T['Z']>=1,T['Z']<=5)]
                    dx = np.diff(zz)/np.diff(xx)
                    dx2 = np.diff(dx)/np.diff(xx[1:len(xx)])
                    iMax = np.where(dx2==max(dx2))[0]
                    iMax = iMax+2
                    toe_mc = iMax
                    
                    
                    toes_pd.append(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_pd])
                    toes_mc.append(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_mc])
                    toes_ml.append(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_ml])
                    y.append(T['Y'][0])
                    x_pd.append(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_pd])
                    x_mc.append(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_mc])
                    x_ml.append(T['X'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_ml])            
                    z_pd.append(T['Z'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_pd])
                    z_mc.append(T['Z'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_mc])
                    z_ml.append(T['Z'][np.logical_and(T['Z']>=1,T['Z']<=5)][toe_ml])
                    T_all.append(T)
            
                    
            if pp == 0:
                toes_pd_pre = toes_pd
                toes_mc_pre = toes_mc
                x_pd_pre = x_pd
                x_mc_pre = x_mc
                x_ml_pre = x_ml
                z_pd_pre = z_pd
                z_mc_pre = z_mc
                z_ml_pre = z_ml
                T_all_pre = T_all
            elif pp==1:
                toes_pd_post = toes_pd
                toes_mc_post = toes_mc
                x_pd_post = x_pd
                x_mc_post = x_mc
                x_ml_post = x_ml
                z_pd_post = z_pd
                z_mc_post = z_mc 
                z_ml_post = z_ml
                T_all_post = T_all
                
        
        axx[ic].hist((np.array(z_pd_post)-np.array(z_pd_pre))/-(np.array(x_pd_post)-np.array(x_pd_pre)),bins=list(np.linspace(-2,1,100))+[np.inf],density=True,cumulative=True,histtype='step',color='r',label='pd')
        axx[ic].hist((np.array(z_mc_post)-np.array(z_mc_pre))/-(np.array(x_mc_post)-np.array(x_mc_pre)),bins=list(np.linspace(-2,1,100))+[np.inf],density=True,cumulative=True,histtype='step',color='b',label='mc')
        axx[ic].hist((np.array(z_ml_post)-np.array(z_ml_pre))/-(np.array(x_ml_post)-np.array(x_ml_pre)),bins=list(np.linspace(-2,1,100))+[np.inf],density=True,cumulative=True,histtype='step',color='g',label='ml')
        axx[ic].hist(BT_manual_all,bins=list(np.linspace(-2,1,100))+[np.inf],density=True,cumulative=True,histtype='step',color='k',label='manual')
        if ic==0:
            fig.legend(['pd','mc','ml','manual'],loc='upper center',ncol=4)
        axx[ic].plot((0,0),(0,1),'--',color='grey',linewidth=1)
        axx[ic].set_xlim(-1,1.5)
        axx[ic].set_xlabel('BT')
        axx[ic].set_title(names[ic])
        if ic!=0 and ic!=4:
            axx[ic].set_yticklabels([])

