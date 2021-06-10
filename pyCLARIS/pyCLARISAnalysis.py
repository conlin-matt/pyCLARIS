

# Standard library imports #
import math
import os

# 3rd party imports #
import alphashape
import numpy as np
import numpy_groupies as npg
import pdal
import pptk
from pybeach.beach import Profile
from scipy import ndimage
from scipy.interpolate import interp1d

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
            None: do not crop the point cloud
            [minX,maxX,minY,maxY]: a list of bounding coordinates in FRF coordinates

    returns:
        None, but will create an FRF.las point cloud in the directory specified or in the
              directory of the file that was passed.
    '''


    dirname, filename = os.path.split(os.path.abspath(__file__)) # Get the directory to this file #
    
    if not os.path.isfile(lasDirec_or_lasFile):
        files = sorted(os.listdir(lasDirec_or_lasFile))
        files = [lasDirec_or_lasFile+'/'+i for i in files if "Store" not in i and "2" in i]
        multiple = True
    else:
        files = [lasDirec_or_lasFile]
        multiple = False

    if croper:
        if croper=='frf':
            bounds = [-100,10000,0,1000]
        elif croper=='5km':
            bounds = [-100,10000,-800,3995]
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
        
        yy = np.arange(min(data[0]['Y']),max(data[0]['Y']),dy)

        transects = []
        x_window = .25
        for y in yy:
            
            data_t = data[0][abs(data[0]['Y']-y)<dy/2]
            data_t = data_t[data_t['Classification']==0]
            try:
                xs = np.arange(min(data_t['X']),max(data_t['X']),x_window)
            except:
                data_t_f = []
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


class scarpManager():
    
    def __init__(self,thresh_vertChange=0.5,thresh_longshoreContinuity=100,thresh_minimumElev=1.5,thresh_slope_after=35):
        self.thresh_vertChange = thresh_vertChange
        self.thresh_longshoreContinuity = thresh_longshoreContinuity
        self.thresh_minimumElev=thresh_minimumElev
        self.thresh_slope_after = thresh_slope_after

    def calcBTOverBf(self,xx,yy,dsm_pre,dsm_post,T,IDmethod,regions_agg=None,regions_deg=None,file_pre=None,file_post=None):
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
            IDmethod: The method to use for scarp toe identification. Options are 'manual','ml','mc','rr','pd'
            regions_agg and _deg: (optional) The extraction of change regions from the DoD can be quite time consuming, so if you have previously
                                  created and saved these regions, you can input them and skip the step of computing them here.
            file_pre and _post: (optional, only needed if IDmethod="manual") Full path+filename to las point cloud files pre and post-storm

        returns:
            BTBf: The BT/Bf ration computed at each transect within the change region(s). This is a list of len=number of change regions.
        '''

        def idScarpRegions(regions_agg,regions_deg):
            '''
            Identify regions of scarp retreat from a DEM of Difference using a rules-based approach.
            '''
            
            # Find change regions above the minimum vertChange threshold #
            if regions_agg and regions_deg:
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
                if all(np.array(z)>self.thresh_minimumElev):
                    regions_deg2.append(reg)

            return regions_deg2


        def extractScarpTransects(scarpRegions):

            # Extract transects in the IDd scarp region that meet the angle threshold #
            T_scarps = []
            for scarp in scarpRegions:
                i_int = []
                for i in range(0,len(T[1])):
                    if min(scarp[:,1])<=T[1][i]['Y'][0]<=max(scarp[:,1]):
                        i_int.append(i)
                    T_int = [[T[0][i] for i in i_int],[T[1][i] for i in i_int]]
                    
                i_scarp = []
                for i in range(0,len(T_int[0])):
                    i_use = np.logical_and(T_int[1][i]['X']>=min(scarp[:,0]),T_int[1][i]['X']<=max(scarp[:,0]))
                    x = T_int[1][i]['X'][i_use]
                    z = T_int[1][i]['Z'][i_use]
                    dzdx = np.diff(z)/np.diff(x)
                    ang = np.degrees(np.tan(dzdx))
                    if max(abs(ang))>50:
                        i_scarp.append(i_int[i])

                T_scarp = [[T[0][i] for i in i_scarp],[T[1][i] for i in i_scarp]]
                T_scarps.append(T_scarp)

            return T_scarps


        def idScarpToesOnTransects(T_scarps,method):

            toes_all = []
            for T_scarp in T_scarps:
            
                toes = [np.empty([0,3]),np.empty([0,3])]
                for day in range(0,2):
                    for t in range(0,len(T_scarp[0])):
                        x = T_scarp[day][t]['X'][30:-1]
                        z = T_scarp[day][t]['Z'][30:-1]
                        pb = Profile(x,z)
                        toe = eval('pb.predict_dunetoe_'+method+'()')
                        toe_x = x[toe[0]]
                        toe_y = T_scarp[day][t]['Y'][0]
                        toe_z = z[toe[0]]
                        toes[day] = np.vstack([toes[day],np.hstack([toe_x,toe_y,toe_z])])

                toes_all.append(toes)

            return toes_all


        def idScarpToesManually(file_pre,file_post,scarpRegions,T_scarps):
            for i in scarpRegions:
                print('IDd scarp region is y='+str(min(reg[i][:,1]))+' to y='+str(max(reg[i][:,1]))+'. Click on the scarp toe in the window and press Enter, then do it again in the new window.')

            pc_bef = pcManager(file_pre)
            scarp_bef = pc_bef.grabFeatures(len(scarpRegions),'rgb')

            pc_aft = pcManager(file_post)
            scarp_aft = pc_aft.grabFeatures(len(scarpRegions),'rgb')

            toes_all = []
            for ii in range(0,len(scarpRegions)):
                scarps = [scarp_bef[ii],scarp_aft[ii]]
                xis = []
                zis = []
                yys = [T_scarps[0][0][i]['Y'][0] for i in range(0,len(T_scarps[0][0]))]
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
                
        

        def calcBf(T_scarps,toes):
            
            Bf_all = []
            for ii in range(0,len(T_scarps)):

                Bf = []
                for t in range(0,len(T_scarps[ii][0])):
                    x_beta = T_scarps[ii][0][t]['X'][T_scarps[ii][0][t]['X']>=toes[ii][0][t,0]]
                    z_beta = T_scarps[ii][0][t]['Z'][T_scarps[ii][0][t]['X']>=toes[ii][0][t,0]]
                    icept,m = utils.linearRegression(x_beta[~np.isnan(z_beta)],z_beta[~np.isnan(z_beta)])
                    Bf.append(-m)

                Bf_all.append(Bf)

            return Bf_all


        def calcBT(toes):

            BT_all = []
            for ii in range(0,len(toes)):
                BT = -((toes[ii][1][:,2]-toes[ii][0][:,2])/(toes[ii][1][:,0]-toes[ii][0][:,0]))
                BT_all.append(BT)

            return BT_all

        def calcRatio(Bf,BT):

            BTBf_all = []
            for ii in range(0,len(Bf)):
                BTBf = BT[ii]/Bf[ii]
                BTBf_all.append(BTBf)

            return BTBf_all

            
        self.scarpRegions = idScarpRegions(regions_agg,regions_deg)
        self.T_scarps = extractScarpTransects(self.scarpRegions)
        if IDmethod is not 'manual':
            self.scarpToes = idScarpToesOnTransects(self.T_scarps,IDmethod)
        else:
            self.scarpToes = idScarpToesManually(file_pre,file_post,self.scarpRegions,self.T_scarps)
        self.Bf = calcBf(self.T_scarps,self.scarpToes)
        self.BT = calcBT(self.scarpToes)
        self.BTBf = calcRatio(self.Bf,self.BT)
        return self.BTBf





