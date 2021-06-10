# Standard library imports #
import os

# 3rd party inputs #
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.interpolate import griddata

# Project imports
import pyCLARIS.pyCLARISAnalysis as claris



date_pre = '20210318'
date_post = '20210323'
savedVars = False
analysisLen = 1
direcs = ['/Users/frfuser/Documents/pyCLARIS_project/data/'+i for i in sorted(os.listdir('/Users/frfuser/Documents/pyCLARIS_project/data')) if date_pre in i or date_post in i]

if analysisLen == 1:
    xx = np.arange(50,120,1)
    yy = np.arange(0,1000,1)
elif analysisLen == 5:
    xx = np.arange(-200,120,1)
    yy = np.arange(-1000,4000,1)
    

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
                croperU = None
            claris.createFRFLas(direc,croper=croperU)
            os.rename(direc+'/FRF.las',direc+'/FRF_'+str(analysisLen)+'km.las')

##        # Create the PC object #
##        pc = claris.pcManager(direc+'/FRF_'+str(analysisLen)+'km.las') 
##
##        # Grid the PC into a DSM and a colored image #
##        print('Making DSM')
##        dsm = pc.gridPC_Slocum(xx,yy,z_val='z',function='min')
##
##        # Pull xshore transects from the PC #
##        print('Making transects')
##        transects = pc.createTransects(dy=5)
##
##        # Save the DSM and transects #
##        dsms.append(dsm)
##        T.append(transects)
##
##    for val in ['dsms','T']:
##        with open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_'+val+'.pkl','wb') as f:
##            pickle.dump(eval(val),f)


else:
    
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_dsms.pkl','rb'); dsms = pickle.load(f)
    f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/'+date_pre+'-'+date_post+'_'+str(analysisLen)+'km_T.pkl','rb'); T = pickle.load(f)
    



# ID scarps #
print('IDing scarps')
scarpLab = claris.scarpManager(thresh_vertChange=0.5,thresh_longshoreContinuity=100,thresh_minimumElev=1.5,thresh_slope_after=35)
BTBf = scarpLab.calcBTOverBf(xx,yy,dsms[0],dsms[1],T,'pd',
                             regions_agg=None,regions_deg=None,file_pre=direcs[0]+'/FRF_'+str(analysisLen)+'km.las',file_post=direcs[1]+'/FRF_'+str(analysisLen)+'km.las')


# Save to catelog #
f = open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/ScarpRetreatCatelog.pkl','rb'); catelog = pickle.load(f)
data_dict = {'DSMs':dsms,'Transects':T,'Scarp Transects':scarpLab.T_scarps,'Bf':scarpLab.Bf,'BT':scarpLab.BT,'BT/Bf':scarpLab.BTBf}
catelog[date_pre+'-'+date_post] = data_dict
with open('/Users/frfuser/Documents/pyCLARIS_project/Analyses/BTBf/data/ScarpRetreatCatelog.pkl','wb') as f:
    pickle.dump(catelog,f)






def plotBTOverBf():

    plt.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(3.5,4.5))
    ax1 = plt.axes([.15,.7,.8,.25])
    ax2 = plt.axes([.15,.4,.8,.25])
    ax3 = plt.axes([.35,.1,.3,.2])
    
    ax1.plot(scarpLab.scarpToes[0][0][:,1],scarpLab.BT[0],'-s')
    ax1.plot(scarpLab.scarpToes[0][0][:,1],scarpLab.Bf[0],'-s')
    ax1.set_ylim(0,0.3)
    ax1.legend(['$B_T$','$B_f$'])
    ax1.set_xticklabels([])
##    ax1.set_title('Method: '+method,fontweight='normal',pad=0)
    
    ax2.plot(scarpLab.scarpToes[0][0][:,1],scarpLab.BTBf[0],'-s')
    ax2.set_ylim(-1,3)
    ax2.text(90,2.5,'$B_T/B_f$')
    ax2.set_xlabel('FRF Y (m)')

    ax3.hist(scarpLab.BTBf[0],bins=np.arange(-1,3,0.2),density=True)
    ax3.set_xlabel('$B_T/B_f$')
    ax3.set_xticks([0,1,2])
    ax3.set_ylabel('Density')
    plt.show()










