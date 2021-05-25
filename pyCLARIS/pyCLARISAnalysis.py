

# Standard library imports #
import math
import os

# 3rd party imports #
import alphashape
import numpy as np
import numpy_groupies as npg
import pdal
import pptk
from scipy import ndimage
from scipy.interpolate import interp1d



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
        files = [lasDirec_or_lasFile+'/'+i for i in files if "Store" not in i]
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
            rgb_or_I: Switch to color the point cloud via RGB values ('rgb') or Intensity ('I') values

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
            z_val: the field to grid. Options are: 'z','rgb', or 'intensity'
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

        # Determine the field(s) to grid based on user input #
        if z_val=='z':
            fieldL = ['Z']
        elif z_val=='rgb':
            fieldL = ['Red','Green','Blue']
        elif z_val=='intensity':
            fieldL=['Intensity']
        else:
            raise ValueError("z_val must be one of 'z','rgb',or 'intensity'")

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


        
    def createTransects(self,dy=5):

        '''
        Create cross-shore transects directly from the classified point cloud at specified longshore spacing.

        args:
            dy: longshore transect spacing

        returns:
            transects: a list where each entry is a transect. Each entry is a 2-element list, where the first element
                       is the raw point data (though only points with classification=0 i.e. ground points) used to create
                       the transect, which is all points within dy/2 of the transect's longshore location. The second element
                       is the smoothed transect. Both elements are arrays structured the same way as the data returned from the
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



    def lasViewer(self,data,rgb_or_I):
        self.xyz = np.transpose(np.array([data[0]['X'],data[0]['Y'],data[0]['Z']]))

        if rgb_or_I == 'rgb':
            self.rgb = np.transpose(np.array([data[0]['Red']/256/255,data[0]['Green']/256/255,data[0]['Blue']/256/255]))
        elif rgb_or_I == 'I':
            self.rgb = data[0]['Intensity']
        else:
            raise ValueError("Please input an arg for rgb_or_I: 'rgb' if you want the point cloud colored by RGB values, 'I' if by intensity values")
            
        self.v = pptk.viewer(self.xyz,self.rgb)
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




