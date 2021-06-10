import datetime
import matplotlib.pyplot as plt
import numpy as np
from pyCLARIS import coastalGeoUtils as utils


waveObj = utils.NDBCWaveRecord(44056,2020)
waves = waveObj.get()

bdate = 20130306
edate = 20130320

waves_small1 = waves.loc[waves['mo'] == int(str(bdate)[4:6])]
waves_small = waves_small1.loc[np.logical_and(waves_small1['day']>=int(str(bdate)[6:8]),waves_small1['day']<= int(str(edate)[6:8]))]
tWave = [datetime.datetime(int(i[0]),int(i[1]),int(i[2]),int(i[3]),int(i[4])) for i in np.array(waves_small)]


plt.plot(tWave,waves_small['wvht (m)']),plt.show()
