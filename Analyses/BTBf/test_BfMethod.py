'''
Code to test the impact of method to calculate Bf: end point vs. linear regression. For each
saved trasect, the slope is recalculated using the lr method, and saved as another variable.
To toggle between making the figure using ep or lr method, on line 78 change the Bf variable 
(used for ep method) to Bf_lr (used for linear regression).
''' 


import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
# Project imports
import pyCLARIS.pyCLARISAnalysis as claris
import pyCLARIS.coastalGeoUtils as utils
 
DEMGridSize = '0.5' # '1' or '0.5' #      
files = sorted([i for i in os.listdir('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/') if '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in i])
dates_b = [20130306, 20170918, 20170922, 20180303, 20190904, 20190910, 20191011, 20191114, 20200910, 20210318]
dates_e = [20130320, 20170922, 20170929, 20180309, 20190910, 20190924, 20191015, 20191119, 20200925, 20210323]
names = ['2013NE','Jose','Maria','Riley','Dorian','Humberto','2019NE1','2019NE2','Teddy','2021NE']
results_byStorm = pd.DataFrame(columns=['Name','results'])


Bf_all = []
Bf_lr_all = []
for i in range(0,len(dates_b)):
        
    results = pd.DataFrame(columns=['BT','Bf','$Z_{Toe}$ (m)','Post-storm toe elev (m)','Toe dx (m)','Toe dz (m)',
                            '$\B_{Dune}$',
                            'F (m)','Beach Width (m)','Impact Duration (hr)',
                            '$\Delta V_B$ ($m^3/m$)','$\Delta V_D$ ($m^3/m$)','Dune Height (m)'])
    bdate = dates_b[i]
    edate = dates_e[i]
    ws = '/Users/frfuser/Documents/pyCLARIS_project/data/FRFWaves/'+str(bdate)[0:4]+'_17m.txt'
    hydro = utils.hydroLab(bdate,edate,station_wl=8651370,station_waves=ws,buoyDepth=17.8)
    
    file = [iii for iii in files if str(dates_b[i])+'-' in iii and '_manual_transects_'+DEMGridSize+'mx'+DEMGridSize+'m' in iii]
    if len(file)>0:
        file = file[0]
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+file,'rb')
        scarpResults=pickle.load(f)
        T_scarps = scarpResults.T_scarps
        scarpToes = scarpResults.scarpToes
        BT = scarpResults.BT
        Bf = scarpResults.Bf
        
        # Recalculate Bf for each transect using linear regression
        Bf_lr = []
        for scarp in range(0,len(T_scarps)):
                Bf_lr_1 = []
                for t in range(0,len(T_scarps[scarp][0])):
                    x = T_scarps[scarp][0][t]['X']
                    z = T_scarps[scarp][0][t]['Z']
                    toe = scarpToes[scarp][0][t]
                    b,yhat,r2,p_values = utils.linearRegression(x[x>toe[0]],z[x>toe[0]])
                    m = -b[1]
                    Bf_lr_1.append(m)
                Bf_lr_1 = np.array(Bf_lr_1)
                Bf_lr.append(Bf_lr_1)
        
        Bf_storm = np.empty([0,])
        Bf_lr_storm = np.empty([0,])
        for ii in range(0,len(Bf)):
            Bf1 = Bf[ii]
            Bf_lr1 = Bf_lr[ii]
            Bf_storm = np.hstack([Bf_storm,Bf1])
            Bf_lr_storm = np.hstack([Bf_lr_storm,Bf_lr1])
        
        Bf_all.append(Bf_storm)
        Bf_lr_all.append(Bf_lr_storm)
        
        
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
                                          '$Z_{Toe}$ (m)':toe_pre,
                                          'Post-storm toe elev (m)':toe_post,
                                          'Toe dx (m)':dx,
                                          'Toe dz (m)':dz,
                                          '$B_{Dune}$':slope,
                                          'F (m)':F,
                                          'Beach Width (m)':width,
                                          'Impact Duration (hr)':impact_dur,
                                          '$\Delta V_B$ ($m^3/m$)':vol_beach,
                                          '$\Delta V_D$ ($m^3/m$)':vol_dune,
                                          'Dune Height (m)':dHeight},
                                         ignore_index=True)
                
        name = names[i]
        results_byStorm = results_byStorm.append({'Name':name,'results':results},ignore_index=True)
    else:
        pass
        
BT_all1 = []
ii = -1
for i in names:
    ii+=1
    try:
        vals = results_byStorm.loc[results_byStorm['Name'] == i]['results'][ii]['BT']
    except:
        ii-=1
    else:
        BT_all1.append(vals)
BT_all = [item for sublist in BT_all1 for item in sublist]


varsToPlot = ['Bf']  
fig,ax = plt.subplots(3,len(varsToPlot),sharey=True,figsize=(2.75,5))
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
        except:
            ii-=1
        else:
            vals_all1.append(vals)
    vals_all = [item for sublist in vals_all1 for item in sublist]
    
    for s in range(0,len(vals_all1)):
        ax[1,iii].plot(vals_all1[s],BT_all1[s],'.',markersize=2)
    ax[1,iii].set_xlabel(colorVal)
    ax[1,iii].set_ylabel('BT')
    
    # coefs,yhat,r2,p_values = utils.linearRegression(np.array(vals_all)[~np.isnan(vals_all)],np.array(BT_all)[~np.isnan(vals_all)])
    
    from sklearn.linear_model import LinearRegression
    from scipy import stats
    reg = LinearRegression().fit(np.array(vals_all)[~np.isnan(vals_all)].reshape(-1,1),
                                 np.array(BT_all)[~np.isnan(vals_all)].reshape(-1,1))
    R2 = reg.score(np.array(vals_all)[~np.isnan(vals_all)].reshape(-1,1),
                                 np.array(BT_all)[~np.isnan(vals_all)].reshape(-1,1))
    coef = reg.coef_
    inter = reg.intercept_
    # Calculate p-values. Method taken from top answer to this SO question https://stackoverflow.com/questions/27928275/find-p-value-significance-in-scikit-learn-linearregression #
    x = np.array(vals_all)[~np.isnan(vals_all)].reshape(-1,1)
    y = np.array(BT_all)[~np.isnan(vals_all)].reshape(-1,1)
    params = np.append(inter[0],coef[0])
    predictions = reg.predict(x)
    newX = np.append(np.ones((len(x),1)), x.reshape(-1,1), axis=1)
    MSE = (sum((y-predictions)**2))/(len(newX)-len(newX[0]))
    var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
    sd_b = np.sqrt(var_b)
    ts_b = params/ sd_b
    p_values =[2*(1-stats.t.cdf(np.abs(i),(len(newX)-len(newX[0])))) for i in ts_b] 
    
    
    if p_values[1]>+0.01:
        t = ax[1,iii].text(min(ax[1,iii].get_xlim()),-0.5,'$R^2=$'+str(round(R2,2))+'\n$p=$'+str(round(p_values[1],3)),ha='left')
    else:
        a = p_values[1]
        t = ax[1,iii].text(min(ax[1,iii].get_xlim()),-0.5,'$R^2=$'+str(round(R2,2))+'\n'+r"$p = {0:0.2e}$".format(a),ha='left')
    t.set_bbox(dict(facecolor='red', alpha=0.3))

    
        
    
    
    from matplotlib.patches import Rectangle
    binsize = 0.1
    binedges = np.arange(-.6,.5,binsize)
    iBins = np.digitize(np.array(BT_all),binedges)
    iPos = np.where(np.array(BT_all)>0)[0]
    iNeg = np.where(np.array(BT_all)<0)[0]
           
    val_pos_mean = np.nanmean(np.array(vals_all)[iPos].astype(float))
    val_neg_mean = np.nanmean(np.array(vals_all)[iNeg].astype(float))
    
    stat,p_val = stats.ks_2samp(np.array(vals_all)[iPos][~np.isnan(np.array(vals_all)[iPos])],np.array(vals_all)[iNeg][~np.isnan(np.array(vals_all)[iNeg])])
     
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

            
    ax[0,iii].set_xlim(ax[1,iii].get_xlim())
    m1 = ax[0,iii].plot(val_pos_mean,min(binedges),'b*',alpha=1)
    m2 = ax[0,iii].plot(val_neg_mean,min(binedges),'r*',alpha=1)
    m1[0].set_clip_on(False)
    m2[0].set_clip_on(False)
    fig.legend([m1[0],m2[0]],['Mean for (+)BT','Mean for (-)BT'],loc='upper center',frameon=False,ncol=2,handletextpad=0.1)
    
    a = p_val
    t = ax[0,iii].text(min(ax[0,iii].get_xlim()),-0.5,r"$p = {0:0.2e}$".format(a),ha='left')
    t.set_bbox(dict(facecolor='red', alpha=0.3))
    
    
    # Do the storm-average plot #   
    h = []
    BT_mean_all = []
    val_mean_all = []
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
        
        ax[2,iii].plot( (val_mean,val_25),(BT_mean,BT_mean),'k',linewidth=1 )
        ax[2,iii].plot( (val_mean,val_75),(BT_mean,BT_mean),'k',linewidth=1 )
        ax[2,iii].plot( (val_mean,val_mean),(BT_mean,BT_25),'k',linewidth=1 )
        ax[2,iii].plot( (val_mean,val_mean),(BT_mean,BT_75),'k',linewidth=1 )
        hh = ax[2,iii].plot(val_mean,BT_mean,'s',alpha=alpha_scaled,markersize=14)
        h.append(hh[0])
        ax[2,iii].set_xlabel(colorVal)
        ax[2,iii].set_ylabel('BT')
        # if iii==0:
        # ax[0,iii].set_xlim(ax[1,iii].get_xlim())
        # else:
        ax[2,iii].set_xlim(ax[1,iii].get_xlim())
        
        val_mean_all.append(val_mean)
        BT_mean_all.append(BT_mean)
        
    # coefs,yhat,r2,p_values = utils.linearRegression(val_mean_all,BT_mean_all)
    
    reg = LinearRegression().fit(np.array(val_mean_all).reshape(-1,1)[~np.isnan(val_mean_all)],np.array(BT_mean_all).reshape(-1,1)[~np.isnan(val_mean_all)])
    R2 = reg.score(np.array(val_mean_all).reshape(-1,1)[~np.isnan(val_mean_all)],np.array(BT_mean_all).reshape(-1,1)[~np.isnan(val_mean_all)])
    coef = reg.coef_
    inter = reg.intercept_
    # Calculate p-values. Method taken from top answer to this SO question https://stackoverflow.com/questions/27928275/find-p-value-significance-in-scikit-learn-linearregression #
    x = np.array(val_mean_all).reshape(-1,1)[~np.isnan(val_mean_all)]
    y = np.array(BT_mean_all).reshape(-1,1)[~np.isnan(val_mean_all)]
    params = np.append(inter[0],coef[0])
    predictions = reg.predict(x)
    newX = np.append(np.ones((len(x),1)), x.reshape(-1,1), axis=1)
    MSE = (sum((y-predictions)**2))/(len(newX)-len(newX[0]))
    var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
    sd_b = np.sqrt(var_b)
    ts_b = params/ sd_b
    p_values =[2*(1-stats.t.cdf(np.abs(i),(len(newX)-len(newX[0])))) for i in ts_b] 
    
    
    if p_values[1]>+0.01:
        t = ax[2,iii].text(min(ax[1,iii].get_xlim()),-0.5,'$R^2=$'+str(round(R2,2))+'\n$p=$'+str(round(p_values[1],3)),ha='left')
    else:
        a = p_values[1]
        t = ax[2,iii].text(min(ax[1,iii].get_xlim()),-0.5,'$R^2=$'+str(round(R2,2))+'\n'+r"$p = {0:0.2e}$".format(a),ha='left')
    t.set_bbox(dict(facecolor='red', alpha=0.3))
        
    fig.legend(h,['2013NE','Jose','Riley','Dorian','2019NE1','2019NE2','Teddy','2021NE'],loc='lower center',
               frameon=False,handletextpad=0.1,fontsize=8,ncol=3)
        
                
        


        
Bf = [item for sublist in Bf_all for item in sublist]       
Bf_lr = [item for sublist in Bf_lr_all for item in sublist]       

fig,ax = plt.subplots(1)
ax.plot(Bf,Bf_lr,'k.',markersize=2)   
ax.plot(ax.get_xlim(),ax.get_xlim(),'r')
plt.axis('equal')
ax.set_xlabel('Bf via end-point slope')
ax.set_ylabel('Bf via linear regression')
ax.set_yticks(ax.get_xticks())
   
        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        