# Standard library imports #

# 3rd party inputs #
import matplotlib.pyplot as plt
import pickle
from pybeach.beach import Profile
import numpy as np
from scipy.interpolate import interp1d

# Project imports
import pyCLARIS.coastalGeoUtils as utils
import pyCLARIS.pyCLARISAnalysis as claris


###############################################
### Manually extract scarp positions #
###############################################
##pc = claris.pcManager('/Users/frfuser/Documents/pyCLARIS_project/data/20210318_FRF5km/FRF.las')
##scarp_bef = pc.grabFeatures(1,'rgb')


#############################################
# Examine and quantify changes to scarp at southern end of FRF property #
#############################################

y_min = 90 # Extents for the southern FRF scarp #
y_max = 230

f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/BeachSurfaceChange/data/T_5km.pkl','rb')
T = pickle.load(f)

T1_small = [i for i in T[0] if y_min<=i['Y'][0]<=y_max]
T2_small = [i for i in T[1] if y_min<=i['Y'][0]<=y_max]
T_small = [T1_small,T2_small]

us = []
sigmas = []
for method in ['manual','mc','pd']:

    
    if method is "manual":
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/ScarpChange/data/scarp_bef.pkl','rb')
        scarp_bef = pickle.load(f)
        f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/ScarpChange/data/scarp_aft.pkl','rb')
        scarp_aft = pickle.load(f)
        
        scarps = [scarp_bef,scarp_aft]

        xis = []
        zis = []
        yys = [T_small[0][i]['Y'][0] for i in range(0,len(T_small[0]))]
        toes = []
        for scarp in scarps:
            
            f = interp1d(scarp[0][:,1],scarp[0][:,0],bounds_error=False,fill_value=np.nan)
            xi = f(yys)

            f = interp1d(scarp[0][:,1],scarp[0][:,2],bounds_error=False,fill_value=np.nan)
            zi = f(yys)

            toe = np.transpose(np.vstack([xi,yys,zi]))
            toes.append(toe)


    else:
        
        toes = [np.empty([0,3]),np.empty([0,3])]
        for day in range(0,2):
            for t in range(0,len(T_small[0])):
                x = T_small[day][t]['X'][30:-1]
                z = T_small[day][t]['Z'][30:-1]
                pb = Profile(x,z)
                toe = eval('pb.predict_dunetoe_'+method+'()')
                toe_x = x[toe[0]]
                toe_y = T_small[day][t]['Y'][0]
                toe_z = z[toe[0]]
                toes[day] = np.vstack([toes[day],np.hstack([toe_x,toe_y,toe_z])])

    # Plot the IDd scarp toe #
##    for i in range(0,len(T_small[0])):
##        fig,ax = plt.subplots(1)
##        ax.plot(T_small[0][i]['X'],T_small[0][i]['Z'],'k')
##        ax.plot(toes[0][i,0],toes[0][i,2],'k*')
##        
##        ax.plot(T_small[1][i]['X'],T_small[1][i]['Z'],'grey')
##        ax.plot(toes[1][i,0],toes[1][i,2],'*',color='grey')
##        ax.set_title(method)
##        plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/ScarpChange/figs/ScarpBasePickValidation/'+method+'/'+str(i)+'.png')
##        fig.show()
##        plt.pause(1)
##        plt.close('all')


    # Calc the scarp retreat angle #
    beta_T = -((toes[1][:,2]-toes[0][:,2])/(toes[1][:,0]-toes[0][:,0]))
    
    # Calc the pre-storm slope of the profile from IDd dune toe to end offshore profile extent #
    beta_f = []
    for t in range(0,len(T_small[0])):
        x_beta = T_small[0][t]['X'][T_small[0][t]['X']>=toes[0][t,0]]
        z_beta = T_small[0][t]['Z'][T_small[0][t]['X']>=toes[0][t,0]]
        icept,m = utils.linearRegression(x_beta[~np.isnan(z_beta)],z_beta[~np.isnan(z_beta)])
        beta_f.append(-m)

    # Calc the ratio #
    rBeta = beta_T/beta_f

    # Calc stats of the ratio distribution #
    u = np.mean(rBeta)
    us.append(u)
    sigma = np.std(rBeta)
    sigmas.append(sigma)

    def plotResults():
        plt.rcParams.update({'font.size': 8})
        fig = plt.figure(figsize=(3.5,4.5))
        ax1 = plt.axes([.15,.7,.8,.25])
        ax2 = plt.axes([.15,.4,.8,.25])
        ax3 = plt.axes([.35,.1,.3,.2])
        
        ax1.plot(toes[0][:,1],beta_T,'-s')
        ax1.plot(toes[0][:,1],beta_f,'-s')
        ax1.set_ylim(0,0.3)
        ax1.legend(['$B_T$','$B_f$'])
        ax1.set_xticklabels([])
        ax1.set_title('Method: '+method,fontweight='normal',pad=0)
        
        ax2.plot(toes[0][:,1],rBeta,'-s')
        ax2.set_ylim(-1,3)
        ax2.text(90,2.5,'$B_T/B_f$')
        ax2.set_xlabel('FRF Y (m)')

        ax3.hist(rBeta,bins=np.arange(-1,3,0.2),density=True)
        ax3.set_xlabel('$B_T/B_f$')
        ax3.set_xticks([0,1,2])
        ax3.set_ylabel('Density')
        plt.savefig('/Users/frfuser/Documents/pyCLARIS_project/Analyses/March2021Storm/ScarpChange/figs/rBeta_'+method+'.png')
        plt.show()

    plotResults()
        







    





