# Standard library imports #
import math

# Third-party imports #
import calendar
import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize
from scipy.interpolate import interp1d



class NOAAWaterLevelRecord():
    
    '''
    Class to interact with NOAA water level records. Class is initialized with desired station and record time.
    
    Inputs:
        station: (int) NDBC station from which to extract data
        bdate: (int or str): the date for the start of the record in the form yyyymmddHHMMSS
        ebdate: (int or str): the date for the end of the record in the form yyyymmddHHMMSS
    
    Methods:
        get(): Download water level record of initialized duration from initialized station. 
        atTime(): Get the water level at a specific time. Must have used get() first.
        
    '''
    
    def __init__(self,station,bdate,edate):
        self.station=station
        self.bdate = bdate
        self.edate = edate
        
    
    def get(self):
       
        product = 'hourly_height&application=UF'
        product2 = 'predictions&application=UF'
        
        api = 'https://tidesandcurrents.noaa.gov/api/datagetter?product='+product+'&begin_date='
        url = api+str(self.bdate)+'&end_date='+str(self.edate)+'&datum=NAVD&station='+str(self.station)+'&time_zone=lst_ldt&units=metric&format=csv'
        api2 = 'https://tidesandcurrents.noaa.gov/api/datagetter?product='+product2+'&begin_date='
        url2 = api2+str(self.bdate)+'&end_date='+str(self.edate)+'&datum=NAVD&station='+str(self.station)+'&time_zone=lst_ldt&units=metric&interval=h&format=csv'
        
        
        dat_obs = pd.read_csv(url)
        dat_pred = pd.read_csv(url2)
        
        wl = pd.DataFrame({'Time':dat_obs['Date Time'],'wl_obs':dat_obs[' Water Level'],'wl_pred':dat_pred[' Prediction']})
        
        self.wl = wl
        return wl
    
    def atTime(self,wl,date):
        d = pd.to_datetime(self.wl['Time'])
        # Turn the date times into numeric values (seconds since the first observation)
        time_int = []
        for i in d:
            td = i-d[0]
            time_int.append(td.total_seconds())
            
        # Get the numeric value for the desired date #
        td = pd.to_datetime(date)-d[0]
        timeWant = td.total_seconds()
        
        # Interpolate wl at that time #
        f = interp1d(time_int,self.wl['wl_obs'])
        wli = float(f(timeWant))
        
        return wli




class NDBCWaveRecord():
    
    '''
    Class to auto-download data from NDBC station. Class can return downloaded data and
    can extract a parameter value at a user-input time.
    
    Inputs:
        station: (int) NDBC station from which to extract data
        years: (int or list of ints) Year(s) for which to extract data
        
    Methods:
        get(): Downlod the record.
        atTime(): Get the value of a specified parameter at a certain time.
    '''
    
    def __init__(self,station,years):
        self.station = station
        self.years = years
        
    def get(self):
        '''
        Returns a dataframe containing the data.
        '''
        
        if type(self.years) is int:
            self.years = [self.years]
        
        allDat = None  
        for yr in self.years:
            
            if yr != datetime.datetime.now().year:
                url = 'https://www.ndbc.noaa.gov/view_text_file.php?filename='+str(self.station)+'h'+str(yr)+'.txt.gz&dir=data/historical/stdmet/'
                dat = pd.read_csv(url,sep=' ',delimiter=' ',header=2,index_col=False,usecols=[0,1,2,3,4,9,11,13,15],
                                  names=['yr','mo','day','hr','mm','wvht (m)','DPD (sec)','APD (sec)','MWD (degT)'])
            else:
                for i in range(1,datetime.datetime.now().month):
                    try:
                        datetime_object = datetime.datetime.strptime(str(i), "%m")
                        month_name = datetime_object.strftime("%b")
                        
                        url = 'https://www.ndbc.noaa.gov/view_text_file.php?filename='+str(self.station)+str(i)+str(yr)+'.txt.gz&dir=data/stdmet/'+month_name+'/'
               
                        dat1 = pd.read_csv(url,sep=' ',delimiter=' ',header=1,index_col=False,usecols=[0,1,2,3,4,8,9,10,11],skipinitialspace=True,
                                           names=['yr','mo','day','hr','mm','wvht (m)','DPD (sec)','APD (sec)','MWD (degT)'])
                        if i == 1:
                            dat = dat1
                        else:
                            dat = dat.append(dat1)
                        
                    except:
                        break
                
            if allDat is not None:
                allDat = allDat.append(dat)
            else:
                allDat = dat
         
        allDat.set_index(np.arange(0,len(allDat)),inplace=True)    
        return allDat
        
    def atTime(self,dat,date,param):
        
        if param==1:
            p = 'wvht (m)'
        elif param==2:
            p = 'DPD (sec)'
        elif param==3:
            p = 'APD (sec)'
        elif param==4:
            p = 'MWD (degT)'
        
        dtimes = []
        for i in range(0,len(dat)):
            if dat['mo'][i]<10:
                dummo = '0'
            else:
                dummo = ''
            if dat['day'][i]<10:
                dumday = '0'
            else:
                dumday = ''
            if dat['hr'][i]<10:
                dumhr = '0'
            else:
                dumhr = ''
            if dat['mm'][i]<10:
                dummm = '0'
            else:dummm = ''
            
            dtimes.append(str(dat['yr'][i])+dummo+str(dat['mo'][i])+dumday+str(dat['day'][i])+dumhr+str(dat['hr'][i])+dummm+str(dat['mm'][i]))
            
        d = pd.to_datetime(dtimes)
        # Turn the date times into numeric values (seconds since the first observation)
        time_int = []
        for i in d:
            td = i-d[0]
            time_int.append(td.total_seconds())
            
        # Get the numeric value for the desired date #
        td = pd.to_datetime(date)-d[0]
        timeWant = td.total_seconds()
        
        # Interpolate wl at that time #
        f = interp1d(time_int,dat[p])
        pi = float(f(timeWant))
        
        return pi
    
    
def calcTWL(bdate,edate,B,station_wl,station_waves):
    
    # Get the water level record #
    wlObj = NOAAWaterLevelRecord(station_wl,bdate,edate)
    wl = wlObj.get()
    
    # Get the wave record and slice it to only include the desired time #
    bdate_yr = int(str(bdate)[0:4])
    edate_yr = int(str(edate)[0:4])
    if int(edate_yr)-int(bdate_yr)==0:
        years = bdate_yr
    else:
        years = [i for i in np.arange(bdate_yr,edate_yr+1)]
    waveObj = NDBCWaveRecord(station_waves,years)
    waves = waveObj.get()
    tWave = [datetime.datetime(int(i[0]),int(i[1]),int(i[2])) for i in np.array(waves)]
    bdate_datestr = datetime.datetime(int(str(bdate)[0:4]),int(str(bdate)[4:6]),int(str(bdate)[6:8]))
    edate_datestr = datetime.datetime(int(str(edate)[0:4]),int(str(edate)[4:6]),int(str(edate)[6:8]))
    waves = waves.loc[np.logical_and(np.array(tWave)>=bdate_datestr,np.array(tWave)<=edate_datestr)]

    # Interpolate the hourly water level observations to the same sampling times as the waves #
    def toTimestamp(d):
        return calendar.timegm(d.timetuple())
    
    tWL = [datetime.datetime(int(i[0:4]),int(i[5:7]),int(i[8:10]),int(i[11:13]),int(i[14:16])) for i in np.array(wl['Time'])]
    tWave = [datetime.datetime(int(i[0]),int(i[1]),int(i[2]),int(i[3]),int(i[4])) for i in np.array(waves)]

    tsWL = [toTimestamp(i) for i in tWL]
    tsWave = [toTimestamp(i) for i in tWave]
    wli = np.interp(tsWave,tsWL,list(wl['wl_obs']))

    # Calculate R2% #
    H = waves['wvht (m)']
    T = waves['DPD (sec)']
    L = (9.81*np.square(T))/(2*math.pi) # APPROXIMATION #

    setup = 0.35*beta*np.sqrt(H*L)
    swash = np.sqrt(H*L*((0.563*(beta**2))+0.004))/2
    TWL_mag = wli+(1.1*(setup+swash))

    TWL = np.transpose(np.vstack([tWave,list(TWL_mag)]))

    return TWL
                      

def linearRegression(x,y):
    '''
    b1 = slope
    b0 = intercept
    '''
    
    x = np.array(x)
    y = np.array(y)
    
    n = np.size(x)
    m_x,m_y = np.mean(x),np.mean(y)
    SS_xy = np.sum(y*x) - n*m_y*m_x
    SS_xx = np.sum(x*x) - n*m_x*m_x
    b_1 = SS_xy/SS_xx
    b_0 = m_y-b_1*m_x

    return b_0,b_1






