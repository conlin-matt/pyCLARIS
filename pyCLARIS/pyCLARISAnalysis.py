

# Standard library imports #
import math
import os

# 3rd party imports #
import alphashape
import matplotlib.pyplot as plt
import numpy as np
import numpy_groupies as npg
import pdal
import pickle
import pptk
from pybeach.beach import Profile
from scipy import ndimage
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import shapely

# Project imports #
from . import coastalGeoUtils as utils

def createFRFLas(lasDirec_or_lasFile,croper='frf'):
    '''
    Function to take las CLARIS data tile(s) and create a las point cloud file
    where the point cloud is in FRF coordinates.
    If multiple tiles, pass a directory containing (only) the tiles. If a single
    tile or file, pass the full file name (directory+file name). Creates an FRF.las
    point cloud.
    args:
        lasDirec_or_lasFile: Directory containing multiple tiles or full file name of single
                             file. Pass as a string.
        croper: A variable describing how you would like the point cloud cropped. Options are:
            'frf': crop the point cloud to the FRF property (i.e. 0 to 1000 m FRF Y)
            '5km': crop the point cloud to the typical 5 km CLARIS stretch (~-1000-4000 FRF Y)
            None: do not crop the point cloud
            [minX,maxX,minY,maxY]: a list of bounding coordinates in FRF coordinates
    returns:
        None, but will create an FRF.las point cloud in the directory specified or in the
              directory of the file that was passed.
    '''


    dirname, filename = os.path.split(os.path.abspath(__file__)) # Get the directory to this file #
    dirname = dirname.replace('\\','/')
    
    if not os.path.isfile(lasDirec_or_lasFile):
        files = [i for i in sorted(os.listdir(lasDirec_or_lasFile)) if '.las' in i]
        files = [lasDirec_or_lasFile+'/'+i for i in files if "Store" not in i and "2" in i]
        multiple = True
    else:
        files = [lasDirec_or_lasFile]
        multiple = False

    if croper:
        if croper=='frf':
            bounds = [-100,10000,0,1000]
        elif croper=='5km':
            # bounds = [-100,10000,-800,3995]
            bounds = [-100,10000,0,4000]
            
        else:
            bounds = croper
    else:
        bounds = [-100,10000,-800,3995]
        
    # Apply the FRF transformation pipeline to each tile, only keep tiles that cover the FRF #
    for lasFile in files:
    
        ## Crop the point cloud to only include points within FRF property and save as new point cloud ##
        json = """
        [
            {
                "type":"readers.las",
                "filename":"""+'"'+lasFile+'"'+""" 
            },
            {
                "type":"filters.python",
                "script":"""+'"'+dirname+'/rotatePC.py'+'"'+""",
                "function":"rotatePC",
                "module":"foo_bar"
            },
            {
                "type":"filters.crop",
                "bounds":"(["""+str(bounds[0])+""","""+str(bounds[1])+"""],["""+str(bounds[2])+""","""+str(bounds[3])+"""])"
            },
            {
                "type":"writers.las",
                "filename":"""+'"'+lasFile.split('.')[0]+'_frf.las'+'"'+"""
            }
        ]
        """

            
        pipeline = pdal.Pipeline(json)
        numPts = pipeline.execute()
        
        
        if numPts == 0: # If no points in the file, remove the newly created file #
            os.remove(lasFile.split('.')[0]+'_frf.las')

        del json
        del pipeline
        del numPts

    # Merge all the _FRF tiles (if there are multiple) into a single final FRF.las point cloud #
    if multiple:
        files_frf = sorted(os.listdir(lasDirec_or_lasFile))
        files_frf = [lasDirec_or_lasFile+'/'+i for i in files_frf if '_frf' in i]        

        if len(files_frf)>0:
            if len(files_frf)>1:
                string_block = ''
                for i in files_frf:
                    if i == files_frf[0]:
                        string_block += i+'",\n            '
                    elif i != files_frf[-1]:
                        string_block += '"'+i+'",\n        '
                    else:
                        string_block += '"'+i        
        
                json = """
                [
                    """+'"'+string_block+'",'+"""
                    {
                        "type":"filters.merge"
                    },
                    """+'"'+lasDirec_or_lasFile+'/FRF.las'+'"'+"""
                ]
                """
        
                pipeline = pdal.Pipeline(json)
                pipeline.execute()
        
                for i in files_frf:
                    os.remove(i)
            else:
                os.rename(files_frf[0],files_frf[0].replace(files_frf[0].split('/')[-1],'FRF.las'))
                
        else:
            print('WARNING: no data were found in desired region- no FRF file was created')
        
    else:
        os.rename(lasFile.split('.')[0]+'_frf.las',lasDirec_or_lasFile.replace(lasDirec_or_lasFile.split('/')[-1],'FRF.las'))
        


class pcManager():

    '''
    Class to efficiently and easily(?) manage large CLARIS point clouds. Importantly, the PC data is only loaded
    into memory if the loadData() method is used- no other methods require loading the data into memory.
    Methods:
        loadData()
        viewPC()
        grabFeatures()
            lasViewer
            idFeature
        gridPC_Slocum
        createTransects
    Attributes:
        lasFile: the point cloud file used to create this object
        pipeline: TODO NEED TO REMOVE
    Created by Matt Conlin, 05/2021, UF and FRF
    '''
        
    
    def __init__(self,lasFile):
        self.lasFile = lasFile

        self.run()
        
        
    def run(self):

        '''
        Creates the point cloud object upon calling the class
        '''

        # Read the file in using a simple PDAL pipeline #
        json = """
        [
            {
                "type":"readers.las",
                "filename":"""+'"'+self.lasFile+'"'+""" 
            }
        ]
        """

        self.pipeline = pdal.Pipeline(json)
        self.pipeline.execute()
         

    def loadData(self):
        '''
        Method to load the point cloud data into a variable.
        args:
            None
        returns:
            data: A numpy array containing the data. Fields can be accessed via string indicies,
            e.g. x = data['X'], red = data['Red']
        '''
        
        data = self.pipeline.arrays
              
        return data[0]
    
    
    def viewPC(self,rgb_or_I='rgb'):
        '''
        Method to view the point cloud in a viewer window. Method leverages PDAL's viewer capability.
        To navigate the viewer:
        click+drag to rotate, shift+click+drag to translate, scroll wheel to zoom in/out
        
        args:
            rgb_or_I: Switch to color the point cloud via RGB values ('rgb') or Intensity ('I') values.
                        Note: Not all point clouds have RGB info associated with them. If the 'rgb' flag is passed
                              for a pc that does not contain RGB info, the point cloud will be shown colored by
                              elevation by default.
        
        returns:
            None
        
        '''
        
        # View the point cloud #
        data = self.pipeline.arrays 
        self.lasViewer(data,rgb_or_I)
        

    def grabFeatures(self,numFeatures,rgb_or_I):
        '''
        Method to manually extract points from point cloud via clicking. Method leverages PDAL's viewer and
        point selector capabilities. You can extract as many features as you want at once.
        To use:
        cntrl+click on a series of points to represent a feature (e.g. dune toe) and press Enter. Then repear for numFeatures times.
        
        args:
            numFeatures: number of features to select.
            rgb_or_I: Switch to color the point cloud via RGB values ('rgb') or Intensity ('I') values
                        Note: Not all point clouds have RGB info associated with them. If the 'rgb' flag is passed
                              for a pc that does not contain RGB info, the point cloud will be shown colored by
                              elevation by default.
        
        returns:
            feat: a list of numFeatures (e.g. 3) nx3 arrays of point coordinates
        '''
        
        self.numFeatures = numFeatures

        data = self.pipeline.arrays 
        self.lasViewer(data,rgb_or_I)

        # Save the XYZ of identified features #
        self.features = []
        for i in range(0,self.numFeatures):
            self.idFeature()
            self.features.append(self.feat_xyz)

        feat = self.features
        return feat
            

    def gridPC_Slocum(self,xx,yy,z_val='z',function='min'):
        '''
        Function to grid the classified point cloud using simple nearest neighbor interpolation following the strategy
        devised by Richie Slocum. This function is essentially a direct Python translation of his Matlab function
        for this specific 2d case (his function allows for arbitrary number of dimensions). Folks at the FRF know
        where to find the Matlab function. The gridding strategy contained herein requires seconds for multi-million
        point point clods, as opposed to hours or more (or a computer crash) for something like griddata.
        args:
            xx: vector of grid x values, created by e.g. np.arange(minx,maxx,dx)
            yy: vector of grid y values, created by e.g. np.arange(miny,maxy,dy)
            z_val: the field to grid. Options are: 'z','rgb','intensity', or 'time'
            function: the function to calculate the gridded value from the group of data points for each grid note.
                  Can be a variety of things, see (https://github.com/ml31415/numpy-groupies#available-functions) for all available
        returns:
            z_grid: Interpolated field values over the grid. If z_val is 'z' or 'intensity', z_grid is a 2d array. If
                    z_val is "rgb", z_grid is a 3-page 2d array (i.e. dimension m,n,3) for use in pyplot.imshow.
        example (elevation):
            xx = np.arange(0,300,1)
            yy = np.arange(-100,1000,1)
            dsm = pc.gridPC_Slocum(xx,yy,z_val='z',function='min') --> pc is a point cloud object created by calling this class on a las/laz file i.e. pc = pyCLARISAnalysis.pcManager(lasFile)
            matplotlib.pyplot.pcolor(yy,xx,np.transpose(dsm),vmin=-1,vmax=8,cmap='gist_earth')
        example (rgb):
            xx = np.arange(0,300,1)
            yy = np.arange(-100,1000,1)            
            im = pc.gridPC_Slocum(xx,yy,z_val='rgb',function='mean') --> pc is a point cloud object created by calling this class on a las/laz file i.e. pc = pyCLARISAnalysis.pcManager(lasFile)
            matplotlib.pyplot.imshow(im,extent=(min(xx),max(xx),min(yy),max(yy)))
        '''

##        xx = np.arange(-50,300,1)
##        yy = np.arange(-100,1000,1)
            
##        # Example test data #
##        x = np.random.randint(0,100,(1000000,))
##        y = np.random.randint(0,100,(1000000,))
##        z = np.multiply(x,y)
##        X = np.transpose(np.vstack([x,y,z]))
##
##        xgi = np.arange(0,101)
##        ygi = np.arange(0,101)
##        Xg = [xgi,ygi]

        # Load the data #
        data = self.pipeline.arrays[0]
        data = data[data['Classification']==0] # Filter the data to only include ground points #
        if len(data)==0: # If there is no classification, just use all points #
            data = self.pipeline.arrays[0]

        # Determine the field(s) to grid based on user input #
        if z_val=='z':
            fieldL = ['Z']
        elif z_val=='rgb':
            fieldL = ['Red','Green','Blue']
        elif z_val=='intensity':
            fieldL=['Intensity']
        elif z_val=='time':
            fieldL=['GpsTime']
        else:
            raise ValueError("z_val must be one of 'z','rgb','intensity',or 'time'")

        all_vals = []
        for field in fieldL:
            
            X = np.transpose(np.vstack([data['X'],data['Y'],data[field]]))
            Xg = [xx,yy]

            # Get rid of data points outside the grid region #
            Xfilt = X
            for iDim in range(0,2):
                dx = np.mean(np.diff(Xg[iDim])) # Grid spacing in dimension #
                lowXg = Xg[iDim][0] # lowest val #
                highXg = Xg[iDim][-1] # Highest val #
                Xfilt = Xfilt[np.logical_and(Xfilt[:,iDim]>lowXg-(dx/2),Xfilt[:,iDim]<highXg+(dx/2)),:] # Keep only point between lowest and highest vals in the dimension #

            # Interpolate data point indicies to grid verticies #
            Xind = Xfilt[:,0:2]
            for iDim in range(0,2):
                f = interp1d(Xg[iDim],np.arange(0,len(Xg[iDim])),kind='nearest',fill_value='extrapolate')
                Xind[:,iDim] = f(Xfilt[:,iDim]).astype(int)
            Xind = Xind.astype(int)
            sizeXg = (len(Xg[0]),len(Xg[1]))

            # Get the interpolated values #
            vals1 = npg.aggregate(np.transpose(Xind),Xfilt[:,2],func=function,size=(sizeXg[0],sizeXg[1]),fill_value=np.nan)
            all_vals.append(vals1)

##        # Visually check accuracy #
##        fig,ax = plt.subplots(1)
##        ax.pcolor(xx,yy,np.transpose(vals1),vmin=-2,vmax=8)
##        ax.scatter(X[0:-1:100,0],X[0:-1:100,1],5,X[0:-1:100,2],vmin=-2,vmax=8,edgecolor='k',linewidth=0.2)
##        fig.show()

        if len(all_vals)==1:
            z_grid = np.transpose(all_vals[0])
        else:
            z_grid = np.empty([len(yy),len(xx),3])
            for ar in range(0,3):
                z_grid[:,:,ar] = np.transpose(all_vals[ar])
            z_grid = z_grid/256/255
                

        return z_grid


    def createTransects_fromDSM(self,xx,yy,dsm,dy=5):
        
        minY = min(yy)
        maxY = max(yy)
        dy_dsm = abs(round((minY-maxY)/len(yy),1))
        spacing = dy/dy_dsm
        iUse = np.arange(0,len(yy),spacing)
        
        transects = []
        for i in iUse:
            data_t_f= np.zeros((len(xx)), dtype=[("X","<f8"),("Y","<f8"),("Z","<f8")])
            data_t_f['X'] = xx
            data_t_f['Y'] = np.tile(yy[int(i)],np.shape(data_t_f['X']))
            data_t_f['Z'] = dsm[int(i),:]
            
            transects.append(data_t_f)
            
        return transects
            
    
    def createTransects(self,dy=5,y_min=None,y_max=None):

        '''
        Create cross-shore transects directly from the classified point cloud at specified longshore spacing.
        args:
            dy: longshore transect spacing
        returns:
            transects: a list where each entry is a transect. The array is structured the same way as the data returned from the
                       pdal pipeline.
        example:
            transects = pc.createTransects(dy=5) -> pc is a point cloud object created by calling this class on a las/laz file i.e. pc = pyCLARISAnalysis.pcManager(lasFile)
            fig = plt.figure()
            plt.plot(transects[51]['X'],transects[51]['Z'],'k'),plt.show()
            
        TODO: correct RGB value calc
        TODO: make faster
        
        NOTE: Each smoothed entry in the returned transect list is a data array
        with only fields: X,Y,Z,Intensity,Red,Green,Blue. I also had to change the dtype of the Intensity,R,G,B columns to float
        rather than int to accomodate NaNs. Could pose an issue later on?
        '''

        data = self.pipeline.arrays

        if y_min and y_max:
            yy = np.arange(y_min,y_max,dy)
        else:
            yy = np.arange(min(data[0]['Y']),max(data[0]['Y']),dy)

        transects = []
        x_window = .25
        for y in yy:
            
            data_t = data[0][abs(data[0]['Y']-y)<dy/2]
            data_t = data_t[data_t['Classification']==0]
            if len(data_t)==0: # If there is no classification, just use all points #
                data_t = data[0][abs(data[0]['Y']-y)<dy/2]
            try:
                xs = np.arange(min(data_t['X']),max(data_t['X']),x_window)
                if len(xs)==0:
                    raise ValueError
            except:
                data_t_f= np.zeros((1), dtype=[("X","<f8"),("Y","<f8"),("Z","<f8"),("Intensity","<f8"),("Red","<f8"),("Green","<f8"),("Blue","<f8")])
                data_t_f['X'] = np.nan
                data_t_f['Y'] = np.nan               
                for field in ['Z','Intensity','Red','Green','Blue']:
                    data_t_f[field] = np.nan
            else:
                data_t_f= np.zeros((len(xs)), dtype=[("X","<f8"),("Y","<f8"),("Z","<f8"),("Intensity","<f8"),("Red","<f8"),("Green","<f8"),("Blue","<f8")])
                data_t_f['X'] = [xs[ii] for ii in range(0,len(xs))]
                data_t_f['Y'] = np.tile(y,np.shape(data_t_f['X']))
                for field in ['Z','Intensity','Red','Green','Blue']:
                    data_t_f[field] = [np.mean(data_t[field][abs(data_t['X']-xs[ii])<x_window/2]) for ii in range(0,len(xs))]
                    
    ##            data_t_f['Z'] = [np.mean(data_t['Z'][abs(data_t['X']-xs[ii])<x_window/2]) for ii in range(0,len(xs))]
    ##            data_t_f['Intensity'] = [np.mean(data_t['Intensity'][abs(data_t['X']-xs[ii])<x_window/2]) for ii in range(0,len(xs))]
    ##            data_t_f['Red'] = [np.mean(data_t['Red'][abs(data_t['X']-xs[ii])<x_window/2]) for ii in range(0,len(xs))]
    ##            data_t_f['Green'] = [np.mean(data_t['Green'][abs(data_t['X']-xs[ii])<x_window/2]) for ii in range(0,len(xs))]
    ##            data_t_f['Blue'] = [np.mean(data_t['Blue'][abs(data_t['X']-xs[ii])<x_window/2]) for ii in range(0,len(xs))]

            transects.append(data_t_f)

        return transects



    def lasViewer(self,data,rgb_or_I):
        self.xyz = np.transpose(np.array([data[0]['X'],data[0]['Y'],data[0]['Z']]))

        if rgb_or_I == 'rgb':
            if np.all(data[0][1:1000]['Blue']==65280) or np.all(data[0][1:1000]['Blue']==0): # No RGB info #
                self.rgb = data[0]['Z']
                self.v = pptk.viewer(self.xyz,self.rgb)
                self.v.set(color_map_scale=(0,8))
            else:
                self.rgb = np.transpose(np.array([data[0]['Red']/256/255,data[0]['Green']/256/255,data[0]['Blue']/256/255]))
                self.v = pptk.viewer(self.xyz,self.rgb)
        elif rgb_or_I == 'I':
            self.rgb = data[0]['Intensity']
            self.v = pptk.viewer(self.xyz,self.rgb)
        else:
            raise ValueError("Please input an arg for rgb_or_I: 'rgb' if you want the point cloud colored by RGB values, 'I' if by intensity values")
            

        self.v.set(point_size=.02)

    def idFeature(self):
        self.v.wait()
        feat_i = self.v.get('selected')
        self.feat_xyz = self.xyz[feat_i,:]



def extractChangeAreas(xx,yy,dsm_pre,dsm_post,thresh=0.1):
    '''
    Function to automatically extract contiguous regions of aggradation/degradation between two gridded surfaces
    (e.g. pre and post storm surfaces). Function is a bit slow because the concave hull calculation for each
    region can take some time.
    args:
        xx: the x-values of the grid as a vector (shape=(m,)) from e.g. xx = np.linspace(min_x,max_x,num_x)
        yy: the y-values of the grid as a vector (shape=(m,)) from e.g. yy = np.linspace(min_y,max_y,num_y)
        dsm_pre: first (in time) gridded surface (e.g. pre-storm)
        dsm_post: second (in time) gridded surface (e.g. post-storm)
        thresh: The elevation change threshold which must be exceeded to be considered a change region.
    returns:
        regions_agg: a list of regions of aggradation. Each region is an nx2 array of x,y coordinates representing the boundary of the region.
        regions_deg: a list of regions of degradation. Each region is an nx2 array of x,y coordinates representing the boundary of the region.
    example:
        direcs = ["direc1","direc2"]
        dsms = []
        for direc in direcs:
            pc =  pyCLARIS.pyCLARISAnalysis.pcManager(direc+'/FRF.las') --> See pcManager class above 
            dsm = pc.gridPC_Slocum(numpy.arange(50,120,1),numpy.arange(0,1000,1),z_val='z',function='min')
            dsms.append(dsm)
            regions_agg,regions_deg = pyCLARIS.pyCLARISAnalysis.extractChangeAreas(dsms[0],dsms[1],thresh=0.1)
    
    '''
    
    # Extract areas of change #
    dod = dsm_post-dsm_pre

    blobs = [dod > thresh , dod < -thresh]
    regions_both = []
    for blob in blobs:
        labels, nlabels = ndimage.label(blob)
        regions = []
        for num in range(1,nlabels+1):
            if len(dod[labels==num])>3:
                locs_i = np.where(labels==num)
                locs_x = xx[locs_i[1]]
                locs_y = yy[locs_i[0]]
                locs = np.transpose(np.vstack([locs_x,locs_y]))
                if np.sum(np.diff(locs[:,0])) != 0 and np.sum(np.diff(locs[:,1])): # Errors thrown when locs are a line #
                    
                    alpha_shape = alphashape.alphashape(
                                    locs,
                                    lambda ind, r: 1.0 + any(np.array(locs)[ind][:,0] == 0.0))
                    try:
                        verts1 = alpha_shape.exterior.coords.xy
                    except:
                        alpha_shape = alphashape.alphashape(locs)
                        verts1 = alpha_shape.exterior.coords.xy

                    verts = np.transpose(np.vstack([verts1[0],verts1[1]]))
                
                    regions.append(verts)

##                fig,ax = plt.subplots(1)
##                ax.plot(locs[:,0],locs[:,1],'.')
##                ax.plot(verts[:,0],verts[:,1],'k')
##                fig.show()
                    
        regions_both.append(regions)

    regions_agg = regions_both[0]
    regions_deg = regions_both[1]

    return regions_agg,regions_deg




def exportNetCDF(lasDirec,date,croper='5km'):    
    '''
    Function to export CLARIS data as a NetCDF file using Spicer's pyMakeNetCDF library.
    
    args:
        lasDirec: directory containing (only) .las CLARIS data tiles
        date: The date of the survey. Format should be "yyyy-mm-dd"
        croper: How to crop the data, as specified in creteFRFLas
        
    returns:
        Nothing, but creates a NetCDF called FRFnc.nc in the directory containing the las tiles
        
    NOTE:
        To use this, you will need to git clone https://github.com/SBFRF/pyMakeNetCDF into your site-packages folder (or anywhere on your python path)
    '''
    
    
    from datetime import datetime
    import numpy.ma as ma
    import netCDF4 as nc
    from pyMakeNetCDF import py2netCDF as p2nc
    import pytz
    
    if croper=='5km':
        xFRF = np.arange(-100,200,0.25)
        yFRF = np.arange(-1000,4000,0.25)
    elif croper=='frf':
        xFRF = np.arange(-100,200,0.25)
        yFRF = np.arange(0,1000,0.25)
    else:
        xFRF = np.arange(croper[0],croper[1],0.25)
        yFRF = np.arange(croper[1],croper[3],0.25)
         
    xx,yy = np.meshgrid(xFRF,yFRF)
    
    # Create a point cloud in FRF coords cropped to desired extent, and then call it #
    createFRFLas(lasDirec,croper)
    pc = pcManager(lasDirec+'/FRF.las')
    
    # Get the time #
    local = pytz.timezone("America/New_York")
    naive = datetime.strptime(date+" 12:00:00", "%Y-%m-%d %H:%M:%S")
    local_dt = local.localize(naive, is_dst=None) # This may be an hour off during DST, but doesn't mater much? #
    utc_dt = local_dt.astimezone(pytz.utc)
    time = float(nc.date2num(utc_dt,'seconds since 1970-01-01'))
        
    # Grid the elevations #
    elev = pc.gridPC_Slocum(xFRF,yFRF,z_val='z',function='min')
    elev[np.isnan(elev)] = int(-999)
    # elev = ma.array(elev,mask=np.isnan(elev),fill_value=-999.99) # Make masked array #
    elev = np.reshape(elev,[1,len(elev[:,0]),len(elev[0,:])]) # Reshape to add time dimension #
    
    # Convert FRF xy to lat lon #
    lat,lon = utils.FRF2LL(xx,yy)
    
    # Create the dictionary to be exported #
    data = {'time':time,
            'xFRF':xFRF,
            'yFRF':yFRF,
            'latitude':lat,
            'longitude':lon,
            'elevation':elev}
    
    # Get the yaml metadata files #
    dirname, filename = os.path.split(os.path.abspath(__file__)) # Get the directory to this file #
    dirname = dirname.replace('\\','/')
    dirname = dirname.rpartition('/')[0]+'/yaml'
    varYaml = dirname+'/claris_grid_var.yml'
    globalYaml = dirname+'/claris_Global.yml'
    
    # Create the NetCDF file #
    p2nc.makenc_generic(lasDirec+'/FRF_'+date+'.nc', globalYaml, varYaml, data)
    
 
    
def calcAverageProfile(direcs,analysisLen=5,xi=np.arange(0,1.01,0.01)):
    
    def cropNormalizeGrid(transects1,xi):
        
        import copy       
        
        # Crop, normalize, and grid transects #
        transects = []
        transects_c = []
        transects_c_n = []
        transects_c_n_i = []
        for t in transects1:
            if max(t['Z'])>3.5 and min(t['Z'])<1: # Make sure looking at full profiles #
                t_c = t[np.logical_and(t['Z']>=1,t['Z']<=3)] # Crop the transect #
                
                dx = np.diff(t_c['X'][~np.isnan(t_c['Z'])])
                if len(dx)>0:
                    if max(dx)<0.3: # Make sure there is not a bunch of missing data in the profile #
                        t_c_n = copy.copy(t_c) # Normalize the transect #
                        t_c_n['X'] = (t_c_n['X']-min(t_c_n['X']))/(max(t_c_n['X'])-min(t_c_n['X']))
                        
                        zi = np.interp(xi,t_c_n['X'],t_c_n['Z'])
                        t_c_n_i = np.zeros([len(xi)],dtype=[('X','<f8'),('Z','<f8')])
                        t_c_n_i['X'] = xi
                        t_c_n_i['Z'] = zi
                        
                        
                        transects.append(t)
                        transects_c.append(t_c)
                        transects_c_n.append(t_c_n)
                        transects_c_n_i.append(t_c_n_i)
 
            else:
                pass
            
            
        return transects,transects_c,transects_c_n,transects_c_n_i
    
    
    transects_all = []
    for direc in direcs:
        # Create the cropped point cloud if it doesn't exist #
        if not os.path.exists(direc+'/FRF_'+str(analysisLen)+'km.las'):
            if analysisLen == 1:
                croperU = 'frf'
            else:
                croperU = '5km'
            createFRFLas(direc,croper=croperU)
            os.rename(direc+'/FRF.las',direc+'/FRF_'+str(analysisLen)+'km.las')
            
        # Pull transects #
        file = direc+'/FRF_'+str(analysisLen)+'km.las'
        pc = pcManager(file)
        transects1 = pc.createTransects(dy=5)
        
        transects,transects_c,transects_c_n,transects_c_n_i = cropNormalizeGrid(transects1,xi)
        
        # Average all profiles for this survey #
        transect_mean = np.zeros([len(xi)],dtype=[('X','<f8'),('Z','<f8')])
        transect_mean['X'] = xi
        for x in range(0,len(xi)):
            vals = [i['Z'][x] for i in transects_c_n_i]
            transect_mean['Z'][x] = np.mean(vals)
   
                
        transects_all.append(transect_mean)  
    
    # Average all the mean profiles to get an overall mean profile #
    transect_mean = np.zeros([len(xi)],dtype=[('X','<f8'),('Z','<f8')])
    transect_mean['X'] = xi
    for x in range(0,len(xi)):
        vals = [i['Z'][x] for i in transects_all]
        transect_mean['Z'][x] = np.mean(vals)
    meanProf = transect_mean
    
    return meanProf


def calcAverageProfiles(direcs,analysisLen=5,xi=np.arange(0,1.01,0.01)):
    
    def cropNormalizeGrid(transects1,xi):
        
        import copy       
        
        # Crop, normalize, and grid transects #
        transects = []
        transects_c = []
        transects_c_n = []
        transects_c_n_i = []
        for t in transects1:
            if max(t['Z'])>3.5 and min(t['Z'])<1: # Make sure looking at full profiles #
                t_c = t[np.logical_and(t['Z']>=1,t['Z']<=3)] # Crop the transect #
                
                dx = np.diff(t_c['X'][~np.isnan(t_c['Z'])])
                if len(dx)>0:
                    if max(dx)<0.3: # Make sure there is not a bunch of missing data in the profile #
                        t_c_n = copy.copy(t_c) # Normalize the transect #
                        t_c_n['X'] = (t_c_n['X']-min(t_c_n['X']))/(max(t_c_n['X'])-min(t_c_n['X']))
                        
                        zi = np.interp(xi,t_c_n['X'],t_c_n['Z'])
                        t_c_n_i = np.zeros([len(xi)],dtype=[('X','<f8'),('Z','<f8')])
                        t_c_n_i['X'] = xi
                        t_c_n_i['Z'] = zi
                        
                        
                        transects.append(t)
                        transects_c.append(t_c)
                        transects_c_n.append(t_c_n)
                        transects_c_n_i.append(t_c_n_i)
 
            else:
                pass
            
            
        return transects,transects_c,transects_c_n,transects_c_n_i
    
    # Make a Nx1000 array for the profiles at each longshore location (cols) for each survey (rows)
    transects_all = np.empty([0,1000])
    for direc in direcs:
        # Create the cropped point cloud if it doesn't exist #
        if not os.path.exists(direc+'/FRF_'+str(analysisLen)+'km.las'):
            if analysisLen == 1:
                croperU = 'frf'
            else:
                croperU = '5km'
            createFRFLas(direc,croper=croperU)
            os.rename(direc+'/FRF.las',direc+'/FRF_'+str(analysisLen)+'km.las')
            
        # Pull transects #
        file = direc+'/FRF_'+str(analysisLen)+'km.las'
        pc = pcManager(file)
        yy = np.arange(-1000,4000,0.5)
        transects1 = pc.createTransects(dy=5,y_min=min(yy),y_max=max(yy))
        
        transects_all = np.vstack([transects_all,np.array(transects1)])
    
    # Find the mean profile at each longshore location (col) #
    meanProfs = []
    for t in range(0,len(transects1)):
        transects_here = [transects_all[ii,t] for ii in range(0,len(transects_all[:,0]))]
        transects,transects_c,transects_c_n,transects_c_n_i = cropNormalizeGrid(transects_here,xi)
            
        # Average all profiles at this location #
        transect_mean = np.zeros([len(xi)],dtype=[('X','<f8'),('Z','<f8')])
        transect_mean['X'] = xi
        for x in range(0,len(xi)):
            vals = [i['Z'][x] for i in transects_c_n_i]
            transect_mean['Z'][x] = np.mean(vals)
            
        meanProfs.append(transect_mean)
        
    return meanProfs
   
              

class scarpManager():
    
    def __init__(self,thresh_vertChange=0.5,thresh_longshoreContinuity=100,thresh_minimumElev=1.5,thresh_slope_after=35):
        self.thresh_vertChange = thresh_vertChange
        self.thresh_longshoreContinuity = thresh_longshoreContinuity
        self.thresh_minimumElev=thresh_minimumElev
        self.thresh_slope_after = thresh_slope_after

    
    def create_pyBeachClassifier(self,data_direc):
        
        dat = sorted([i for i in os.listdir(data_direc) if 'scarpResults_manual_transects' in i])
        
        xi = np.arange(-100,200,0.1)
        data_z = np.empty([len(xi),0])
        data_toe=[]
        
        for file in dat:
            f = open(data_direc+'/'+file,'rb')
            scarpResults = pickle.load(f)
            for region in range(0,len(scarpResults.scarpRegions)):
                T = scarpResults.T_scarps[region]
                toes = scarpResults.scarpToes[region]
                for date in range(0,len(T)):
                    for tran in range(0,len(T[date])):
                        f = interp1d(T[date][tran]['X'],T[date][tran]['Z'],fill_value=np.nan,bounds_error=False)
                        zi = f(xi)
                        data_z = np.hstack([data_z,np.reshape(zi,(-1,1))])
                        
                        toe_loc = toes[date][tran][0]
                        i = np.where(abs(xi-toe_loc)==min(abs(xi-toe_loc)))[0][0]
                        data_toe.append(i)
        
        # Create classifier
        from pybeach.support import classifier_support as cs
        clf = cs.create_classifier(xi, np.transpose(data_z), np.array(data_toe), window=40, min_buffer=40, max_buffer=200)
        
        # save classifier to package directory
        import joblib
        import pybeach
        path = os.path.dirname(pybeach.__file__) + '/classifiers/custom_clf.joblib'
        with open(path,'wb') as f:
            joblib.dump(clf, f, protocol=pickle.HIGHEST_PROTOCOL)


    def idDuneToe(self,T,method,clf='mixed_clf'):
        '''
        Method to extract dune toes from transects using a specified automated method from pyBeach.
        Method is pretty much the same as the idScarpToesOnTransects() function in the calcBTOverBf method, 
        which is lame, but oh well.
        
        args:
            T: transects e.g. from pcManager.createTransects()
            method: Method to extract dToe. Options are: 'ml','mc','rr','pd'
            clf: optional, only used if method='ml'. pyBeach clf to use with ml method.
            
        returns:
            toes: a list containing the dune toe locations.
        
        '''

        toes = [np.empty([0,3]),np.empty([0,3])]
        for day in range(0,1):
            for t in range(0,len(T[0])):
                
                xx = T[day][t]['X']
                zz = T[day][t]['Z']
                
                if len(xx)>1:
                
                    try:
                        pb = Profile(xx[np.logical_and(zz>=1,zz<=5)],zz[np.logical_and(zz>=1,zz<=5)])
                    except:
                        toe_x = np.nan
                        toe_y = np.nan
                        toe_z = np.nan
                    else:                                          
                        if method=='ml':
                            toe = pb.predict_dunetoe_ml(clf)
                        elif method=='contour': # Find where transect intersects 3 m contour #
                            xi = np.linspace(min(xx),max(xx),1000)
                            z = np.interp(xi,xx[~np.isnan(zz)],zz[~np.isnan(zz)])
                            x = xi
                            toe = np.where(abs(z-3) == min(abs(z-3)))[0]
                        elif method=='pd':
                            toe = pb.predict_dunetoe_pd(dune_crest=None,shoreline=None)
                        elif method=='mc':                           
                            def smooth(x,window_len=7,window='hanning'):
                                s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
                                #print(len(s))
                                if window == 'flat': #moving average
                                    w=np.ones(window_len,'d')
                                else:
                                    w=eval('np.'+window+'(window_len)')
                            
                                y=np.convolve(w/w.sum(),s,mode='valid')
                                return y[round(window_len/2-1):-round(window_len/2)]
                            
                            xx1 = xx[np.logical_and(zz>=1,zz<=5)]
                            zz1 = zz[np.logical_and(zz>=1,zz<=5)]
                            zz_smooth = smooth(zz1)
                            dx = np.diff(zz1)/np.diff(xx1)
                            dx2 = np.diff(dx)/np.diff(xx1[1:len(xx1)])
                            dx_smooth = np.diff(zz_smooth)/np.diff(xx1)
                            dx2_smooth = np.diff(dx_smooth)/np.diff(xx1[1:len(xx1)])
                            
                            i_subset = np.where(np.logical_and(zz1>=1,zz1<=3))[0]
                            i_subset=i_subset+2
                                                   
                            iMax = np.where(dx2==max(dx2[i_subset]))[0]
                            iMax_smooth = np.where(dx2_smooth==max(dx2_smooth[i_subset]))[0]
                            iMax_use = iMax+2
                            iMax_use_smooth = iMax_smooth+2
                            toe = iMax_use_smooth
        
                            # toe = pb.predict_dunetoe_mc()
                        
                        elif method=='mc_supervised':
                                
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
                              
                            xx = xx[~np.isnan(zz)]
                            zz = zz[~np.isnan(zz)]  
                            
                            xi = np.arange(min(xx),max(xx),1)
                            zi = np.interp(xi,xx,zz)
                            xx = xi
                            zz = zi
                            
                            xx1 = xx#[np.logical_and(zz>=1,zz<=5)]
                            zz1 = zz#[np.logical_and(zz>=1,zz<=5)]
                            zz_smooth = zz1#smooth(zz1)
                            iSub = np.logical_and(zz_smooth>=1.5,zz_smooth<=5)
                            dx = np.gradient(xx1, xx1)  # first derivatives
                            dy = np.gradient(zz_smooth, xx1)                
                            d2x = np.gradient(dx, xx1)  # second derivatives
                            d2y = np.gradient(dy, xx1)                
                            curvature = d2y / ((1 + dy ** 2)) ** 1.5  # curvature     
                            curvature_sub = curvature[iSub]
                            toe = np.where(curvature==max(curvature_sub))[0]
                            
                            # fig,ax = plt.subplots(2,1,sharex=True,figsize=(5,6))
                            # ax[0].plot(xx,zz)
                            # ax[0].plot(xx1,zz_smooth)
                            # ax[0].plot(xx1[toe],zz_smooth[toe],'k.')            
                            # ax[1].plot(xx1[iSub],curvature[iSub])
                            # ax[1].plot(xx1[toe],curvature[toe],'k.')  
                            # plt.pause(0.1)
                            # fig.show()
                            # yn = input('Satisfied with this dune toe pick? y or n: ')
                            # if yn=='y':
                            #     pass
                            # else:
                            #     toe_manual1 = plt.ginput(1)
                            #     toe = np.where(abs(xx1-toe_manual1[0][0])==min(abs(xx1-toe_manual1[0][0])))[0]
                            # plt.close('all')
                        else:
                            toe = eval('pb.predict_dunetoe_'+method+'()')
                        
                    toe_x = xx[toe]#[np.logical_and(zz>=1,zz<=5)][toe]
                    toe_y = T[day][t]['Y'][0]
                    toe_z = zz[toe]#[np.logical_and(zz>=1,zz<=5)][toe]
                    toes[day] = np.vstack([toes[day],np.hstack([toe_x,toe_y,toe_z])])
                else:
                    toe_x = np.nan
                    toe_y = np.nan
                    toe_z = np.nan
                    toes[day] = np.vstack([toes[day],np.hstack([toe_x,toe_y,toe_z])])


        return toes
    
    

    def calcBTOverBf(self,xx,yy,dsm_pre,dsm_post,T,IDmethod,slopeMethod,regions_agg=None,regions_deg=None,file_pre=None,file_post=None,savedFile=None,clf='mixed_clf'):
        '''
        Method to calculate the BT/Bf value for a scarp. Method has a number of sub-methods to do this in parts; the results
        of each sub-method can be accessed as part of the class after this method is called.
        Sub-methods:
            idScarpRegions(): ID regions of likely scarp retreat from a DEM as lonshore-continuous regions of pronounced degradation above a contour
            extractScarpTransects(): Find transects (previously created) that intersect the scarp retreat regions, and keep only those with a steep slope (scarp)
            idScarpToesOnTransects(): Idenitify the toe of the scarp on each transect
            calcBf(): Calculate the beachface slope in front of the scarp toe in the pre DSM
            calcBT(): Calculate the retreat trajectory
            calcRatio(): Calculate the BT/Bf ratio
        args:
            xx: vector of x grid values used to create the DSMs
            yy: vector of y grid values used to create the DSMs
            dsm_pre: The pre-storm dsm
            dsm_post: The post-storm dsm
            T: The transect list creted from the two point clouds
            IDmethod: The method to use for scarp toe identification. Options are 'manual','manual_transects','ml','mc','rr','pd'
            slopeMethod: Method to comput beachface slope. Either linear regression ('lr') or end point ('ep')
            regions_agg and _deg: (optional) The extraction of change regions from the DoD can be quite time consuming, so if you have previously
                                  created and saved these regions, you can input them and skip the step of computing them here.
            file_pre and _post: (optional, only used if IDmethod="manual") Full path+filename to las point cloud files pre and post-storm
            savedFile: (optional, only used if IDmethod='manual') saved file returned by running this class previously. 
                                                                   Useful if changing something and want to retain the dtoe picks already completed.
            clf: (optional, only used if IDmethod='ml') The pyBeach classifier type to use in the ml dune toe prediction. Options are
                                                        'barrier_island_clf','wave_embayed_clf','mixed_clf',or xxxxx_clf (if you have
                                                        made your own using the create_pyBeachClassifier function)
        
        returns:
            BTBf: The BT/Bf ratio computed at each transect within the change region(s). This is a list of len=number of change regions.
        '''

        def idScarpRegions(regions_agg,regions_deg):
            '''
            Identify regions of scarp retreat from a DEM of Difference using a rules-based approach.
            '''
            
            # Find change regions above the minimum vertChange threshold #
            if regions_agg is not None and regions_deg is not None:
                pass
            else:
                regions_agg,regions_deg = extractChangeAreas(xx,yy,dsm_pre,dsm_post,thresh=self.thresh_vertChange)

            # Regions that meet the longshore continuity threshold #
            regions_deg1 = []
            for reg in regions_deg:
                longshoreDist = np.max(reg[:,1])-np.min(reg[:,1])
                if longshoreDist>=self.thresh_longshoreContinuity:
                    regions_deg1.append(reg)

            # Regions that meet elevation threshold #
            regions_deg2 = []
            for reg in regions_deg1:
                z = []
                for i in reg:
                    z.append(dsm_pre[np.where(yy==i[1])[0][0],np.where(xx==i[0])[0][0]])
                if any(np.array(z)[~np.isnan(z)]>self.thresh_minimumElev):
                    regions_deg2.append(reg)

            return regions_deg2


        def extractScarpTransects(scarpRegions):

            # Extract transects in the IDd scarp region that meet the angle threshold #
            T_scarps = []
            for scarp in scarpRegions:
                i_int = []
                for i in range(0,len(T[1])):
                    if min(scarp[:,1])<=T[1][i]['Y'][0]<=max(scarp[:,1]) and min(scarp[:,1])<=T[0][i]['Y'][0]<=max(scarp[:,1]):
                        i_int.append(i)
                    T_int = [[T[0][i] for i in i_int],[T[1][i] for i in i_int]]
                    
                i_scarp = []
                for i in range(0,len(T_int[0])):
                    i_use = np.logical_and(T_int[1][i]['X']>=min(scarp[:,0]),T_int[1][i]['X']<=max(scarp[:,0]))
                    x = T_int[1][i]['X'][i_use]
                    z = T_int[1][i]['Z'][i_use]
                    dzdx = np.diff(z)/np.diff(x)
                    ang = np.degrees(np.tan(dzdx))
                    if np.nanmax(abs(ang))>self.thresh_slope_after:
                        i_scarp.append(i_int[i])

                T_scarp = [[T[0][i] for i in i_scarp],[T[1][i] for i in i_scarp]]
                T_scarps.append(T_scarp)
                
            self.scarpRegions = [self.scarpRegions[i] for i in range(0,len(T_scarps)) if len(T_scarps[i][0])>0]
            T_scarps = [T_scarps[i] for i in range(0,len(T_scarps)) if len(T_scarps[i][0])>0]
                

            return T_scarps


        def idScarpToesOnTransects(T_scarps,method):
            
            toes_all = []
            for T_scarp in T_scarps:
            
                toes = [np.empty([0,3]),np.empty([0,3])]
                for day in range(0,2):
                    for t in range(0,len(T_scarp[0])):
                        
                        xx = T_scarp[day][t]['X']
                        zz = T_scarp[day][t]['Z']
                        
                        pb = Profile(xx[np.logical_and(zz>=1,zz<=5)],zz[np.logical_and(zz>=1,zz<=5)])
                                       
                        if method=='ml':
                            toe = pb.predict_dunetoe_ml(clf)
                            toe = toe[0]
                        elif method=='contour': # Find where transect intersects 3 m contour #
                            xi = np.linspace(min(xx),max(xx),1000)
                            z = np.interp(xi,xx[~np.isnan(zz)],zz[~np.isnan(zz)])
                            x = xi
                            toe = np.where(abs(z-3) == min(abs(z-3)))[0]
                        elif method=='pd':
                            toe = pb.predict_dunetoe_pd(dune_crest=None,shoreline=None)
                        elif method=='mc':                           
                            
                            # Brodie and Spore #
                              def smooth(x,window_len=7,window='hanning'):
                                  s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
                                  #print(len(s))
                                  if window == 'flat': #moving average
                                      w=np.ones(window_len,'d')
                                  else:
                                      w=eval('np.'+window+'(window_len)')
                            
                                  y=np.convolve(w/w.sum(),s,mode='valid')
                                  return y[round(window_len/2-1):-round(window_len/2)]
                            
                              xx = xx[~np.isnan(zz)]
                              zz = zz[~np.isnan(zz)]
                              zz = smooth(zz)
                              x_high = xx[np.where(zz==max(zz))[0][0]]
                              z_high = zz[np.where(zz==max(zz))[0][0]]
                              try:
                                x_low = utils.transectElevIntercept(0.5,xx,zz)
                                z_low = np.interp(x_low,xx,zz)
                              except IndexError:
                                x_low = xx[-1]
                                z_low = zz[-1]
                              m = (z_high-z_low)/(x_high-x_low)
                              b = z_high-(m*x_high)
                              yhat = (m*xx[np.where(zz==z_high)[0][0]:np.where(zz==z_low)[0][0]])+b           
                              dist = yhat-zz[np.where(zz==z_high)[0][0]:np.where(zz==z_low)[0][0]]
                              firstGuess = np.where(dist==max(dist))[0]
                              iSub = np.where(abs(xx-xx[firstGuess])<=5)[0]
                          
                              dx = np.gradient(xx, xx)  # first derivatives
                              dy = np.gradient(zz, xx)
                              d2x = np.gradient(dx, xx)  # second derivatives
                              d2y = np.gradient(dy, xx)
                              curvature = d2y / ((1 + dy ** 2)) ** 1.5  # curvature
                              curvature_sub = curvature[iSub]
                              toe = np.where(curvature==max(curvature_sub))[0]
                              
                              
                              # fig,ax = plt.subplots(2,1,sharex=True)
                              # ax[0].plot(xx,zz)
                              # ax[0].plot(xx[toe],zz[toe],'k.')  
                              # ax[1].plot(xx,curvature)
                              # ax[1].plot(xx[iSub],curvature_sub)
                              # ax[1].plot(xx[toe],curvature[toe],'k.')  
                              # plt.pause(0.1)
                              # fig.show()
                              # yn = input('Satisfied with this dune toe pick? y or n: ')
                              # if yn=='y':
                              #     pass
                              # else:
                              #     toe_manual1 = plt.ginput(1)
                              #     toe = np.where(abs(xx-toe_manual1[0][0])==min(abs(xx-toe_manual1[0][0])))[0]
                              # plt.close('all')
                            
                            
                            # # Max curvature #
                            # def smooth(x,window_len=7,window='hanning'):
                            #     s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
                            #     #print(len(s))
                            #     if window == 'flat': #moving average
                            #         w=np.ones(window_len,'d')
                            #     else:
                            #         w=eval('np.'+window+'(window_len)')
                            
                            #     y=np.convolve(w/w.sum(),s,mode='valid')
                            #     return y[round(window_len/2-1):-round(window_len/2)]   
                            
                            # xx1 = xx[np.logical_and(zz>=1,zz<=4)]
                            # zz1 = zz[np.logical_and(zz>=1,zz<=4)]
                            # zz_smooth = smooth(zz1)
                            # # iSub = np.logical_and(zz_smooth>=2,zz_smooth<=4)
                            # dx = np.gradient(xx1, xx1)  # first derivatives
                            # dy = np.gradient(zz_smooth, xx1)                
                            # d2x = np.gradient(dx, xx1)  # second derivatives
                            # d2y = np.gradient(dy, xx1)                
                            # curvature = d2y / ((1 + dy ** 2)) ** 1.5  # curvature     
                            # # curvature_sub = curvature[iSub]
                            # toe = np.where(curvature==max(curvature))[0]
                            
                         
                        
                        elif method=='mc_supervised':
                            
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
                              
                            xx = xx[~np.isnan(zz)]
                            zz = zz[~np.isnan(zz)]  
                            
                            xi = np.arange(min(xx),max(xx),1)
                            zi = np.interp(xi,xx,zz)
                            xx = xi
                            zz = zi
                            
                            xx1 = xx#[np.logical_and(zz>=1,zz<=5)]
                            zz1 = zz#[np.logical_and(zz>=1,zz<=5)]
                            zz_smooth = zz1#smooth(zz1)
                            iSub = np.logical_and(zz_smooth>=2,zz_smooth<=5)
                            dx = np.gradient(xx1, xx1)  # first derivatives
                            dy = np.gradient(zz_smooth, xx1)                
                            d2x = np.gradient(dx, xx1)  # second derivatives
                            d2y = np.gradient(dy, xx1)                
                            curvature = d2y / ((1 + dy ** 2)) ** 1.5  # curvature     
                            curvature_sub = curvature[iSub]
                            toe = np.where(curvature==max(curvature_sub))[0]
                            
                            # fig,ax = plt.subplots(2,1,sharex=True,figsize=(5,6))
                            # ax[0].plot(xx,zz)
                            # ax[0].plot(xx1,zz_smooth)
                            # ax[0].plot(xx1[toe],zz_smooth[toe],'k.')            
                            # ax[1].plot(xx1[iSub],curvature[iSub])
                            # ax[1].plot(xx1[toe],curvature[toe],'k.')  
                            # plt.pause(0.1)
                            # fig.show()
                            # yn = input('Satisfied with this dune toe pick? y or n: ')
                            # if yn=='y':
                            #     pass
                            # else:
                            #     toe_manual1 = plt.ginput(1)
                            #     toe = np.where(abs(xx1-toe_manual1[0][0])==min(abs(xx1-toe_manual1[0][0])))[0]
                            # plt.close('all')
                            
                        else:
                            toe = eval('pb.predict_dunetoe_'+method+'()')
                        
                        toe_x = xx[toe]#[np.logical_and(zz>=1,zz<=5)][toe]
                        toe_y = T_scarp[day][t]['Y'][0]
                        toe_z = zz[toe]#[np.logical_and(zz>=1,zz<=5)][toe]
                        toes[day] = np.vstack([toes[day],np.hstack([toe_x,toe_y,toe_z])])

                toes_all.append(toes)

            return toes_all


        def idScarpToesManually(file_pre,file_post,scarpRegions,T_scarps):
         
            pc_bef = pcManager(file_pre)
            scarp_bef = pc_bef.grabFeatures(len(scarpRegions),'rgb')

            pc_aft = pcManager(file_post)
            scarp_aft = pc_aft.grabFeatures(len(scarpRegions),'rgb')

            toes_all = []
            for ii in range(0,len(scarpRegions)):
                scarps = [scarp_bef[ii],scarp_aft[ii]]
                xis = []
                zis = []
                yys = [T_scarps[ii][0][i]['Y'][0] for i in range(0,len(T_scarps[ii][0]))]
                toes = []
                for scarp in scarps:
                    
                    f = interp1d(scarp[:,1],scarp[:,0],bounds_error=False,fill_value=np.nan)
                    xi = f(yys)

                    f = interp1d(scarp[:,1],scarp[:,2],bounds_error=False,fill_value=np.nan)
                    zi = f(yys)

                    toe = np.transpose(np.vstack([xi,yys,zi]))
                    toes.append(toe)
                toes_all.append(toes)

            return toes_all
                
        
        def idScarpToesManually_OnTransects(scarpRegions,T_scarps,savedFile=None):
            toes_all = []          
            for ii in range(0,len(scarpRegions)): 
                if savedFile:
                    regions_saved = savedFile.scarpRegions
                    for rr in range(0,len(regions_saved)):
                        if len(scarpRegions[ii])==len(regions_saved[rr]):
                            if sum(scarpRegions[ii].flatten() - regions_saved[rr].flatten())==0:
                                test = 1
                                break
                            else:
                                test = 0
                        else:
                            test = 0
                else:
                    test==0
                        
                if test == 1:
                    toes_all.append(savedFile.scarpToes[rr])
                    self.T_scarps[ii] = savedFile.T_scarps[rr]
                else:
                    toe_pre = np.empty([0,3])
                    toe_post = np.empty([0,3])
                    for i_t in range(0,len(T_scarps[ii][0])):
                        fig,ax = plt.subplots(1)
                        ax.plot(T_scarps[ii][0][i_t]['X'],T_scarps[ii][0][i_t]['Z'])
                        ax.set_title('Pre, line '+str(i_t+1)+' of '+str(len(T_scarps[ii][0])))
                        toe_pre1_1 = plt.ginput(n=1,timeout=-1)
                        plt.close('all')
                        toe_pre1 = [toe_pre1_1[0][0],T_scarps[ii][0][i_t]['Y'][0],toe_pre1_1[0][1]]
                        toe_pre = np.vstack([toe_pre,toe_pre1])
    
                        fig,ax = plt.subplots(1)
                        ax.plot(T_scarps[ii][1][i_t]['X'],T_scarps[ii][1][i_t]['Z'])
                        ax.set_title('Post, line '+str(i_t+1)+' of '+str(len(T_scarps[ii][0])))                    
                        toe_post1_1 = plt.ginput(n=1,timeout=-1)
                        plt.close('all')
                        toe_post1 = [toe_post1_1[0][0],T_scarps[ii][0][i_t]['Y'][0],toe_post1_1[0][1]]
                        toe_post = np.vstack([toe_post,toe_post1])   
                    toes = [toe_pre,toe_post]
                    toes_all.append(toes)
                    del toes
            
            return toes_all
                
        
        def checkRetreatRegions(scarpRegions,T_scarps,toes):
            
            # Check to make sure erosion region is landward of
            # the pre-storm dune toe at each transect. If erosion region
            # is seaward of the pre-storm toe, this is beach erosion
            # and we should throw out that region. #
            iBad = []
            for i in range(0,len(scarpRegions)):
                iBad_t = []
                for ii in range(0,len(T_scarps[i][0])):
                    xx = T_scarps[i][0][ii]['X']
                    yy = T_scarps[i][0][ii]['Y']
                    toe_x = toes[i][0][ii][0]
                    toe_x_post = toes[i][1][ii][0]

                    toe_dx = toe_x_post-toe_x
                    
                    coords = [(xx[i],yy[i]) for i in range(0,len(xx))]
                    transect = shapely.geometry.LineString(coords)                    
                    region = shapely.geometry.Polygon(scarpRegions[i])
                    cross = transect.intersection(region)
                    try: # Throws an error when the region intersects the transect multiple times #
                        cross_minx = min(cross.coords.xy[0])
                    except:
                        mins = [min(cross[i].coords.xy[0]) for i in range(0,len(cross))]
                        cross_minx = min(mins)
        
                    if cross_minx>toe_x or toe_dx>-0.5:
                        iBad_t.append(ii)
     

                for iii in sorted(iBad_t, reverse=True):
                    T_scarps[i][0] = np.delete(T_scarps[i][0],iii)
                    T_scarps[i][1] = np.delete(T_scarps[i][1],iii)
                    toes[i][0] = np.delete(toes[i][0],iii,axis=0)
                    toes[i][1] = np.delete(toes[i][1],iii,axis=0)
                    
            
                
                    
                prop = len(toes[i][0])#len(iBad_t)/len(T_scarps[i][0])
                if prop<5:#>0.5:    
                    iBad.append(i)
                else:
                    pass
            
            if len(iBad)>0:
                scarpRegions = [scarpRegions[i] for i in range(0,len(scarpRegions)) if i not in iBad]
                T_scarps = [T_scarps[i] for i in range(0,len(T_scarps)) if i not in iBad]
                toes = [toes[i] for i in range(0,len(toes)) if i not in iBad]
            else:
                pass
                        
            
            return scarpRegions,T_scarps,toes
            

        def calcBf(T_scarps,toes,method):
            
            Bf_all = []
            for ii in range(0,len(T_scarps)):

                Bf = []
                for t in range(0,len(T_scarps[ii][0])):
                    if method=='lr':
##                        x_beta = T_scarps[ii][0][t]['X'][T_scarps[ii][0][t]['X']>=toes[ii][0][t,0]]
##                        z_beta = T_scarps[ii][0][t]['Z'][T_scarps[ii][0][t]['X']>=toes[ii][0][t,0]]
                        x_beta = T_scarps[ii][0][t]['X'][np.logical_and(T_scarps[ii][0][t]['X']>=toes[ii][0][t,0],T_scarps[ii][0][t]['Z']>=0.36)]
                        z_beta = T_scarps[ii][0][t]['Z'][np.logical_and(T_scarps[ii][0][t]['X']>=toes[ii][0][t,0],T_scarps[ii][0][t]['Z']>=0.36)]
                        try:
                            b,yhat,r2,p_values = utils.linearRegression(x_beta[~np.isnan(z_beta)],z_beta[~np.isnan(z_beta)])
                            m=b[1]
                        except:
                            m = np.nan
                        Bf.append(-m)
                    elif method=='ep':
##                        x_beta = T_scarps[ii][0][t]['X'][T_scarps[ii][0][t]['X']>=toes[ii][0][t,0]]
##                        z_beta = T_scarps[ii][0][t]['Z'][T_scarps[ii][0][t]['X']>=toes[ii][0][t,0]]                        
                        x_beta = T_scarps[ii][0][t]['X'][np.logical_and(T_scarps[ii][0][t]['X']>=toes[ii][0][t,0],T_scarps[ii][0][t]['Z']>=0.36)]
                        z_beta = T_scarps[ii][0][t]['Z'][np.logical_and(T_scarps[ii][0][t]['X']>=toes[ii][0][t,0],T_scarps[ii][0][t]['Z']>=0.36)]
                        try:
                            m = abs((z_beta[-1]-z_beta[0])/(x_beta[-1]-x_beta[0]))
                        except:
                            m = np.nan
                        Bf.append(m)

                Bf_all.append(Bf)

            return Bf_all


        def calcBT(toes):

            BT_all = []
            for ii in range(0,len(toes)):
                BT = -((toes[ii][1][:,2]-toes[ii][0][:,2])/(toes[ii][1][:,0]-toes[ii][0][:,0]))
                BT_all.append(BT)
  
            return BT_all
        
        
        def checkForAnomalies():
            
            # Flag any transects where the dune toe is extracted to have advanced (probably a bad pick) #
            forwardFlag = []
            for ii in range(0,len(self.scarpToes)):
                
                id_bad = np.where(self.scarpToes[ii][1][:,0]-self.scarpToes[ii][0][:,0]>0)[0]
                forwardFlag.append(id_bad)
                
            return forwardFlag
                
                

        def calcRatio(Bf,BT):

            BTBf_all = []
            for ii in range(0,len(Bf)):
                BTBf = BT[ii]/Bf[ii]
                BTBf_all.append(BTBf)

            return BTBf_all

            
        self.scarpRegions = idScarpRegions(regions_agg,regions_deg)
        self.T_scarps = extractScarpTransects(self.scarpRegions)
        
        if IDmethod == 'manual':
            self.scarpToes = idScarpToesManually(file_pre,file_post,self.scarpRegions,self.T_scarps)
        elif IDmethod == 'manual_transects':
            self.scarpToes = idScarpToesManually_OnTransects(self.scarpRegions,self.T_scarps,savedFile)
        else:
            self.scarpToes = idScarpToesOnTransects(self.T_scarps,IDmethod)  
        
        scarpRegions_new,T_scarps_new,toes_new = checkRetreatRegions(self.scarpRegions,self.T_scarps,self.scarpToes) 
        self.scarpRegions = scarpRegions_new; self.T_scarps = T_scarps_new; self.scarpToes = toes_new
        
        self.Bf = calcBf(self.T_scarps,self.scarpToes,slopeMethod)
        self.BT = calcBT(self.scarpToes)
        self.forwardFlag = checkForAnomalies()
        self.BTBf = calcRatio(self.Bf,self.BT)
        return self.BTBf

