

# Standard library imports #
import math
import os

# 3rd party imports #
import numpy as np
import numpy_groupies as npg
import pandas as pd
import pdal
import pptk
from scipy.interpolate import interp1d



class pcManager():

    '''
    Class to efficiently and easily(?) manage large CLARIS point clouds. Importantly, the PC data is only loaded
    into memory if the loadData() method is used- no other methods require loading the data into memory.

    Methods:
        loadData()
        viewPC()
        grabFeatures()
            _lasViewer
            _idFeature
        gridPC_Slocum

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

##            json = """
##            [
##                {
##                    "type":"readers.las",
##                    "filename":"""+'"'+self.lasFile+'"'+""" 
##                },
##                {
##                    "type":"writers.text",
##                    "format":"csv",
##                    "order":"X,Y,Z,Red,Green,Blue,Intensity,Classification",
##                    "keep_unspecified":"false",
##                    "filename":"""+'"'+self.lasFile.split('.')[0]+'.csv'+'"'+"""
##                }
##            ]
##            """

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
    
    
    def viewPC(self):
        '''
        Method to view the point cloud in a viewer window. Method leverages PDAL's viewer capability.

        To navigate the viewer:
        click+drag to rotate, shift+click+drag to translate, scroll wheel to zoom in/out

        args:
            None

        returns:
            None
        
        '''
        
        # View the point cloud #
        data = self.pipeline.arrays 
        self.lasViewer(data)
        

    def grabFeatures(self,numFeatures):
        '''
        Method to manually extract points from point cloud via clicking. Method leverages PDAL's viewer and
        point selector capabilities. You can extract as many features as you want at once.

        To use:
        cntrl+click on a series of points to represent a feature (e.g. dune toe) and press Enter. Then repear for numFeatures times.

        args:
            numFeatures: number of features to select.

        returns:
            feat_xyz: a list of numFeatures (e.g. 3) nx3 arrays of point coordinates
        '''
        
        self.numFeatures = numFeatures

        self.lasViewer(data)

        # Save the XYZ of identified features #
        self.features = []
        for i in range(0,self.numFeatures):
            self.idFeature()
            self.features.append(self.feat_xyz)
            

    def gridPC_Slocum(self,xx,yy,function='min'):
        '''
        Function to grid the point cloud using simple nearest neighbor interpolation following the strategy
        devised by Richie Slocum. This function is essentially a direct Python translation of his Matlab function
        for this specific 2d case (his function allows for arbitrary number of dimensions). Folks at the FRF know
        where to find the Matlab function. The gridding strategy contained herein requires seconds for multi-million
        point point clods, as opposed to hours or more (or a computer crash) for something like griddata.

        args:
            xx: vector of grid x values, created by e.g. np.arange(minx,maxx,dx)
            yy: vector of grid y values, created by e.g. np.arange(miny,maxy,dy)
            function: the function to calculate the gridded value from the group of data points for each grid note.
                  Can be a variety of things, see (https://github.com/ml31415/numpy-groupies#available-functions) for all available

        returns:
            vals1: Interpolated Z values over the grid

        example:
            xx = np.arange(0,300,1)
            yy = np.arange(-100,1000,1)
            dsm = pc.gridPC_Slocum(xx,yy) --> pc is a point cloud object created by calling this class on a las/laz file i.e. pc = pyCLARISAnalysis.pcManager(lasFile)
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
        
        data = self.pipeline.arrays 
        X = np.transpose(np.vstack([data[0]['X'],data[0]['Y'],data[0]['Z']]))
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

##        fig,ax = plt.subplots(1)
##        ax.pcolor(xx,yy,np.transpose(vals1),vmin=-2,vmax=8)
##        ax.scatter(X[0:-1:100,0],X[0:-1:100,1],5,X[0:-1:100,2],vmin=-2,vmax=8,edgecolor='k',linewidth=0.2)
##        fig.show()

        return np.transpose(vals1)


        
    def createTransects(self,dy=5):

        '''
        Create cross-shore transects directly from the point cloud at specified longshore spacing.

        args:
            dy: longshore transect spacing

        returns:
            transects: a list where each entry is a transect. Each entry is a 2-element list, where the first element
                       is the raw point data (though only points with classification=1 i.e. ground points) used to create
                       the transect, which is all points within dy/2 of the transect's longshore location. The second element
                       is the smoothed transect. Both elements are arrays structured the same way was the data returned from the
                       pdal pipeline.

        example:
            transects = pc.createTransects(dy=5) -> pc is a point cloud object created by calling this class on a las/laz file i.e. pc = pyCLARISAnalysis.pcManager(lasFile)
            fig = plt.figure()
            plt.plot(transects[51][0]['X'],transects[51][0]['Z'],'.'),plt.plot(transects[51][1]['X'],transects[51][1]['Z'],'k'),plt.show()
            

        TODO: potentially calculate transect RGB values by dividing the Intensity values directly, i.e.:
        data_t_f['Red'] = [np.mean(data_t['Red'][abs(data_t['X']-xs[ii])<x_window/2]/data_t['Intensity'][abs(data_t['X']-xs[ii])<x_window/2]) for ii in range(0,len(xs))]

        TODO: make faster
        
        NOTE: Each smoothed entry in the returned transect list (i.e. element 2 of each transect entry) is a data array
        with only fields: X,Y,Z,Intensity,Red,Green,Blue. I also had to change the dtype of the Intensity,R,G,B columns to float
        rather than int to accomodate NaNs. Could pose an issue later on?

        '''

        data = self.pipeline.arrays
        
        yy = np.arange(min(data[0]['Y']),max(data[0]['Y']),dy)

        transects = []
        x_window = .5
        for y in yy:
            
            data_t = data[0][abs(data[0]['Y']-y)<dy/2]
            data_t = data_t[data_t['Classification']==0]
            
            xs = np.arange(min(data_t['X']),max(data_t['X']),x_window)

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

            transects.append([data_t,data_t_f])

        return transects





            

            

        
            
            


 

    def lasViewer(self,data):
        self.xyz = np.transpose(np.array([data[0]['X'],data[0]['Y'],data[0]['Z']]))
        self.rgb = np.transpose(np.array([data[0]['Red']/data[0]['Intensity'],data[0]['Green']/data[0]['Intensity'],data[0]['Blue']/data[0]['Intensity']]))

        self.v = pptk.viewer(self.xyz,self.rgb)
        self.v.set(point_size=.02)

    def idFeature(self):
        self.v.wait()
        feat_i = self.v.get('selected')
        self.feat_xyz = self.xyz[feat_i,:]


    
def createFRFLas(lasDirec_or_lasFile):
    '''
    Function to take las CLARIS data tile(s) and create a las point cloud file
    where the point cloud covers only the FRF property and is in FRF coordinates.
    If multiple tiles, pass a directory containing (only) the tiles. If a single
    tile or file, pass the full file name (directory+file name). Creates an FRF.las
    point cloud.

    args:
        lasDirec_or_lasFile: Directory containing multiple tiles or full file name of single
                             file. Pass as a string.

    returns:
        None, but will create an FRF.las point cloud in the directory specified or in the
              directory of the file that was passed.
    '''

    if not os.path.isfile(lasDirec_or_lasFile):
        files = sorted(os.listdir(lasDirec_or_lasFile))
        files = [lasDirec_or_lasFile+'/'+i for i in files]
        multiple = True
    else:
        files = [lasDirec_or_lasFile]
        multiple = False
        
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
                "script":"C:/Users/conli/Documents/FRF/rotatePC.py",
                "function":"rotatePC",
                "module":"foo_bar"
            },  
            {
                "type":"filters.crop",
                "bounds":"([-100,10000],[0,1000])"
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
    else:
        os.rename(lasFile.split('.')[0]+'_frf.las',lasDirec_or_lasFile.replace(lasDirec_or_lasFile.split('/')[-1],'FRF.las'))
        




