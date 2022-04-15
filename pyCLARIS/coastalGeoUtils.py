# Standard library imports #
import math

# Third-party imports #
import calendar
import csv
import datetime
import numpy as np
import pandas as pd
from scipy import optimize,stats
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
    
    def __init__(self,station,bdate,edate):
        self.station = station
        self.bdate=bdate
        self.edate=edate
        
        nyr = int(str(self.edate)[0:4])-int(str(self.bdate)[0:4])+1
        if nyr==1:
            self.years = int(str(self.edate)[0:4])
        else:
            self.years = list(np.arange(int(str(self.bdate)[0:4]),int(str(self.edate)[0:4])+1))
            
        
    def get(self):
        '''
        Returns a dataframe containing the data.
        '''
        
        if type(self.years) is int:
            self.years = [self.years]
        
        allDat = None  
        for yr in self.years:
            
            if yr != 2021:#datetime.datetime.now().year:
                url = 'https://www.ndbc.noaa.gov/view_text_file.php?filename='+str(self.station)+'h'+str(yr)+'.txt.gz&dir=data/historical/stdmet/'
                dat = pd.read_csv(url,sep=' ',delimiter=' ',header=2,index_col=False,usecols=[0,1,2,3,4,5,7,9,11,13,15],
                                  names=['yr','mo','day','hr','mm','wdir (degT)','wspeed (m/s)','wvht (m)','DPD (sec)','APD (sec)','MWD (degT)'])
            else:
                for i in range(1,12):#datetime.datetime.now().month):
                    try:
                        datetime_object = datetime.datetime.strptime(str(i), "%m")
                        month_name = datetime_object.strftime("%b")
                        
                        url = 'https://www.ndbc.noaa.gov/view_text_file.php?filename='+str(self.station)+str(i)+str(yr)+'.txt.gz&dir=data/stdmet/'+month_name+'/'
               
                        dat1 = pd.read_csv(url,sep=' ',delimiter=' ',header=1,index_col=False,usecols=[0,1,2,3,4,5,6,8,9,10,11],skipinitialspace=True,
                                           names=['yr','mo','day','hr','mm','wdir (degT)','wspeed (m/s)','wvht (m)','DPD (sec)','APD (sec)','MWD (degT)'])
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
        # allDat = allDat[allDat['wvht (m)']<90]
        
        # Truncate the dataframe to the bdate and edate #
        waves_d = []
        for i in range(0,len(allDat)):
            yr = str(allDat['yr'][i])
            mo = str(allDat['mo'][i])
            if len(mo)==1:
                mo = '0'+mo
            day = str(allDat['day'][i])
            if len(day)==1:
                day = '0'+day
            ds = yr+mo+day
            waves_d.append(int(ds))
        allDat = allDat.loc[np.logical_and(waves_d>=np.int64(self.bdate),waves_d<=np.int64(self.edate))]
        allDat = allDat.reset_index()
        allDat = allDat.drop('index',axis=1)
        
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
    
 
class FRFWaveRecord():
    
    def __init__(self,file,bdate,edate):
        self.bdate = bdate
        self.edate = edate
        self.file = file
        
    def get(self):
        '''
        Returns a dataframe containing the data.
        '''
        
        data = pd.read_excel(self.file,engine='openpyxl')
        data = data.rename(columns={'Year':'yr',
                              'Month':'mo',
                              'Day':'day',
                              'Hour':'hr',
                              'Min':'mm',
                              'Sec':'sec',
                              'SigWaveHeight_m':'wvht (m)',
                              'PeakWavePeriod_s':'Tp (s)',
                              'MeanWaveDirection_deg':'MWD (degT)'})

        
        # data = pd.read_csv(self.file)
        
        # rows = []
        # csv_header = ['Year','Mo','Da','Time','Gauge','dir','HmO','Period']
        # frame_header = ['yr','mo','day','hrmm','Gauge','MWD (degT)','wvht (m)','DPD (sec)']
        
        # with open(self.file, 'rt') as f_input:
        #     for row in csv.DictReader(f_input, delimiter=' ', fieldnames=csv_header[:-1], restkey=csv_header[-1], skipinitialspace=True):
        #         try:
        #             rows.append([row['Year'],row['Mo'],row['Da'],row['Time'],row['Gauge'],row['dir'],row['HmO'],row['Period']])
        #         except KeyError:
        #             rows.append([row['Year'],row['Mo'],row['Da'],row['Time'],row['Gauge'],row['dir'],row['HmO'],row['Period']])
        
        # data = pd.DataFrame(rows, columns=frame_header)
        # data = data.iloc[5:-1]
        # data = data.reset_index()
        # data = data.drop(columns='index')
        # # Fix period column #
        # period_fill = [i[0] for i in data['DPD (sec)']]
        # data['DPD (sec)'] = period_fill
        # # Create separate hour and minute columns #
        # hr= [int(str(i)[0:2]) for i in data['hrmm']]
        # mm= [int(str(i)[2:4]) for i in data['hrmm']]
        # data = data.drop(columns='hrmm')
        # data.insert(3,'hr',hr)
        # data.insert(4,'mm',mm)
        
        # Truncate the dataframe to the bdate and edate #
        waves_d = []
        for i in range(0,len(data)):
            yr = str(data['yr'][i])
            mo = str(data['mo'][i])
            if len(mo)==1:
                mo = '0'+mo
            day = str(data['day'][i])
            if len(day)==1:
                day = '0'+day
            ds = yr+mo+day
            waves_d.append(int(ds))
        data = data.loc[np.logical_and(waves_d>=np.int64(self.bdate),waves_d<=np.int64(self.edate))]
        data = data.reset_index()
        data = data.drop('index',axis=1)
        
        return data
    
    def atTime(self,dat,date,param):
        
        if param==1:
            p = 'wvht (m)'
        elif param==2:
            p = 'Tp (s)'
        elif param==3:
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
                     


class hydroLab():
    def __init__(self,bdate,edate,station_wl=None,station_waves=None,buoyDepth=None):

        if not station_wl and not station_waves:
            raise ValueError('Please provide a water level station ID, wave buoy number, or both')

        self.buoyDepth=buoyDepth

        if station_wl:
            wlObj = NOAAWaterLevelRecord(station_wl,bdate,edate)
            self.wl = wlObj.get()
            
        if station_waves:
            if type(station_waves) is int: #NDBC#
                waveObj = NDBCWaveRecord(station_waves,bdate,edate)
                self.waves = waveObj.get()
            elif type(station_waves) is str: #FRF file#
                waveObj = FRFWaveRecord(station_waves,bdate,edate)
                self.waves = waveObj.get()
            

    def calcTWL(self,beta,exceedance=2):

        wl = self.wl
        waves = self.waves

        # Interpolate the hourly water level observations to the same sampling times as the waves #
        def toTimestamp(d):
            return calendar.timegm(d.timetuple())
        
        tWL = [datetime.datetime(int(i[0:4]),int(i[5:7]),int(i[8:10]),int(i[11:13]),int(i[14:16])) for i in np.array(wl['Time'])]
        tWave = [datetime.datetime(int(i[0]),int(i[1]),int(i[2]),int(i[3]),int(i[4])) for i in np.array(waves)]

        tsWL = [toTimestamp(i) for i in tWL]
        tsWave = [toTimestamp(i) for i in tWave]
        wli = np.interp(tsWave,tsWL,list(wl['wl_obs']))

        # Calculate R2% #
        H = waves['wvht (m)'].astype(float)
        T = waves['Tp (s)'].astype(float)
        L = []
        for t in T:
            k = newtRaph(t,self.buoyDepth)
            L.append((2*math.pi)/k)
        L = pd.Series(L)
        
        if exceedance==2:
            n_sigma=2
        elif exceedance==7:
            n_sigma=1.5
        elif exceedance==16:
            n_sigma=1
        elif exceedance==23:
            n_sigma=0.75
        elif exceedance==30:
            n_sigma=0.5
        elif exceedance==50:
            n_sigma=0
            
        setup = 0.35*math.tan(beta)*np.sqrt(H*L)
        swash = (np.sqrt(H*L*((0.563*(beta**2))+0.004))/2)*(n_sigma/2)
        
        # R2_manual = 1.1*(setup+swash)
        
        # from py_wave_runup import models
        # beta = np.tile(beta,len(L))
        # for i in range(0,len(L)):
        #     model_sto06 = models.Stockdon2006(Hs=H[i], Lp=L[i], beta=beta[i])   
        #     R2_stockdon = model_sto06.R2
        #     pow18 = models.Power2018(Hs=H[i], Lp=L[i], beta=beta[i], r=0.00075)
        #     R2_power = pow18.R2
        #     hol86 = models.Holman1986(Hs=H[i], Lp=L[i], beta=beta[i])
        #     R2_holman = hol86.R2
        #     niel09 = models.Nielsen2009(Hs=H[i], Lp=L[i], beta=beta[i])
        #     R2_nielson = niel09.R2
        #     rug01 = models.Ruggiero2001(Hs=H[i], Lp=L[i], beta=beta[i])
        #     R2_prugg = rug01.R2
        #     vou12 = models.Vousdoukas2012(Hs=H[i], Lp=L[i], beta=beta[i])
        #     R2_vous = vou12.R2
        #     atk17 = models.Atkinson2017(Hs=H[i], Lp=L[i], beta=beta[i])
        #     R2_atk = atk17.R2
        #     sen11 = models.Senechal2011(Hs=H[i], Lp=L[i], beta=beta[i])
        #     R2_sen = sen11.R2
            
        # fig,ax = plt.subplots(1)
        # ax.bar(['Stockdon','Power','Holman','Nielson','Prugg','Vousdoukas','Atkinson','Senechal'],
        #         [R2_stockdon,R2_power,R2_holman,R2_nielson,R2_prugg,R2_vous,R2_atk,R2_sen])
        # ax.set_title('Hs='+str(round(H[i],2))+',Lp='+str(round(L[i],2))+',Tp='+str(round(T[i],2))+',beta='+str(beta[i]))
        # ax.set_ylabel('R2%')
            
       
        TWL_mag = wli+(1.1*(setup+swash))

        TWL = np.transpose(np.vstack([tWave,list(TWL_mag)]))

        return TWL
    

    def transectTWLIntercept(self,TWL,x,z):

        xx = x[~np.isnan(z)]
        zz = z[~np.isnan(z)]
        TWL = np.zeros(len(xx))+TWL
        idx = np.argwhere(np.diff(np.sign(TWL - zz))).flatten()
        x_intercept = xx[idx]
        return x_intercept[0]
        


def newtRaph(T,h):
    
    '''
    Function to determine k from dispersion relation given period (T) and depth (h) using
    the Newton-Raphsun method.
    '''
    
    if not np.isnan(T):
        L_not = (9.81*(T**2))/(2*math.pi)
        k1 = (2*math.pi)/L_not
    
        def fk(k):
            return (((2*math.pi)/T)**2)-(9.81*k*math.tanh(k*h))
        def f_prime_k(k):
            return (-9.81*math.tanh(k*h))-(9.81*k*(1-(math.tanh(k*h)**2)))
    
        k2 = 100
        i = 0
        while abs((k2-k1))/k1 > 0.01:
              i+=1
              if i!=1:
                  k1=k2
              k2 = k1-(fk(k1)/f_prime_k(k1))
    else:
        k2 = np.nan

    return k2


def calcTransectSlope(x,z):
    '''
    Function to calculate the slope to use in an R2% calculation given a transect.
    Slope is calculate by performing a linear regression of all points below 4 m
    elevation (approx the dune toe?)
    '''
    b0,m = linearRegression(x[z<=3],z[z<=3])

    return(-m)
    


def linearRegression(x,y):
    '''
    b1 = slope
    b0 = intercept
    '''
    
    x = np.array(x)
    y = np.array(y)
    
    # Perform the regression #
    n = np.size(x)
    m_x,m_y = np.mean(x),np.mean(y)
    SS_xy = np.sum(y*x) - n*m_y*m_x
    SS_xx = np.sum(x*x) - n*m_x*m_x
    b_1 = SS_xy/SS_xx
    b_0 = m_y-b_1*m_x
    
    # Calculate R2 #
    m = b_1
    b = b_0
    
    yhat = (x*m)+b
    sst = np.sum((y-np.mean(y))**2)
    ssr = np.sum((yhat-np.mean(y))**2)
    r2 = ssr/sst
    
    # Calculate p-values. Method taken from top answer to this SO question https://stackoverflow.com/questions/27928275/find-p-value-significance-in-scikit-learn-linearregression #
    params = np.append(b,m)
    predictions = yhat
    
    newX = np.append(np.ones((len(x),1)), x.reshape(-1,1), axis=1)
    MSE = (sum((y-predictions)**2))/(len(newX)-len(newX[0]))
    
    
    var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
    sd_b = np.sqrt(var_b)
    ts_b = params/ sd_b
    
    p_values =[2*(1-stats.t.cdf(np.abs(i),(len(newX)-len(newX[0])))) for i in ts_b]    

    return (b_0,b_1),yhat,r2,p_values


def transectElevIntercept(elev,x,z):

    xx = x[~np.isnan(z)]
    zz = z[~np.isnan(z)]
    elev = np.zeros(len(xx))+elev
    idx = np.argwhere(np.diff(np.sign(elev - zz))).flatten()
    x_intercept = xx[idx]
    return x_intercept[0]


def FRF2LL(xFRF,yFRF):
    '''
    Function to convert FRF coordinates to lat/lon. This is a partial translation of
    Kent Hathaway frfCoord.m script.

    args:
        xFRF: (array-like) FRF X values
        yFRF: (array-like) FRF y values
        
    returns:
        lat: latitude values (decimal degrees)
        lon: longitude values (decimal degrees)
    '''
    
    r2d = 180.0 / math.pi
    ALat0=36.1775975              # Origin Lat minutes
    ALon0=75.7496860              # Origin Lon minutes
    DegLat = 110963.35726         # m/deg Lat
    DegLon = 89953.36413          # m/deg long
    GridAngle=18.1465 /r2d

    R = np.sqrt(xFRF**2 + yFRF**2)
    Ang1 = np.arctan2(xFRF, yFRF)          # CW from Y
    Ang2 = Ang1 - GridAngle            
    ALatLeng = np.multiply(R,np.cos(Ang2))
    ALonLeng = np.multiply(R,np.sin(-Ang2))     # neg to go west
    ALat = ALatLeng/DegLat + ALat0    # was 60 * ALatLeng./DegLat + ALat0;
    ALon = ALonLeng/DegLon + ALon0
    
    return ALat,ALon
    
    


