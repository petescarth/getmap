"""
Flask app designed to work with timeseries groundcover data
to produce output rasters of ground cover anomalies and 
time series and spatial percentiles for grazing land management applications

Original script from JRSRP
Modified by peter.scarth@gmail.com to include custom anomalies
and serve via a Flask app.
Reguires timeseries ground cover mosaics in albers equal area to run.

"""

import os,sys, time, shutil, fnmatch, json, datetime

from osgeo import gdal,ogr,osr
import numpy as np
from scipy.stats import percentileofscore
from optparse import OptionParser
from rios import applier
from flask import Flask, request, current_app, make_response,send_from_directory
from functools import update_wrapper

# Seasonal image data location with JRSRP naming convention
IMGDIR = '/rdsi/public/data/landsat/seasonal_fractional_cover/ground_cover/aus/'

# Temp output directory
OUTPUT_DIR = '/mnt/getmap'

######################################################################################
# General Functions
#####################################################################################


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

  
    
def find(pattern, startdir=os.curdir):
    matches = []
    os.path.walk(startdir, findvisitor, (matches, pattern))
    matches.sort()
    return matches


def findvisitor((matches, pattern), thisdir, nameshere):
    for name in nameshere:
        if fnmatch.fnmatch (name, pattern):
            fullpath = os.path.join(thisdir, name)
            matches.append(fullpath)



def writeJson(obs):
    obs.shpfile = os.path.join(obs.outFolder,'clip.json')
    with open(obs.shpfile, 'w') as outfile:
        outfile.write(obs.json)

    return obs


def reprojectShp(obs):
    """
    reprojects shapefile into australian albers
    """
    
    obs.newShp = os.path.join(obs.outFolder,'clip3577.json')
    cmd = 'ogr2ogr -f GeoJSON -t_srs EPSG:3577 %s %s' % (obs.newShp,obs.shpfile)
    os.system(cmd)

    return obs



def getExtent(obs):
    """
    gets the extent  of the transformed shapefile
    """

    obs.DS = ogr.Open(str(obs.newShp))
    obs.layer = obs.DS.GetLayer()
    obs.srs = obs.layer.GetSpatialRef()

    obs.extent = obs.layer.GetExtent()
    (obs.ulx,obs.lrx,obs.lry,obs.uly) = obs.extent

    return obs




def getWindow(obs):
    """
    Extracts portion of raster defined by window
    """

    ## for our mask we subset a cover image

    obs.window = "lztmre_%s_window.img" % obs.site
    obs.window = os.path.join(obs.outFolder,obs.window)

    cmd = """
    gdal_translate -of HFA -b 1 -projwin %f %f %f %f  %s %s
     """ % ( obs.ulx, obs.uly, obs.lrx, obs.lry,
             obs.templateFile, obs.window)
    os.system(cmd)

    return obs





def getPaddockMask(obs):
    """
    creates a mask based on the reprojected shapefile for stats at paddock level.
    """

    obs.paddockMask =   "lztmre_%s_mask_paddock.img" % obs.site
    obs.paddockMask = os.path.join(obs.outFolder,  obs.paddockMask)

    shutil.copyfile(obs.window,obs.paddockMask)


    ## now we loop through the polygon layer and burn things in
    ## field/strata
    strataField = obs.strata

    ## just in case each polygon isn't unique in its strata
    strata = {}
    stratacounter = 1
    for feature in obs.layer:
        thisFeatVal = feature.GetField(strataField)
        #print thisFeatVal
        if not thisFeatVal in strata:
            strata[thisFeatVal] = stratacounter
            stratacounter = stratacounter + 1

    ## now make the mask
    cmd = """
    gdal_rasterize -i -burn 0 -l %s %s %s
    """ % (obs.layer.GetName(), obs.newShp, obs.paddockMask)
    os.system(cmd)

    for stratum in strata:
        cmd = """
    gdal_rasterize -burn %d -where "%s='%s'" -l %s %s %s
    """ % (strata[stratum], strataField, stratum, obs.layer.GetName(), obs.newShp, obs.paddockMask)
        os.system(cmd)


    return obs


def getPropertyMask(obs):
    """
    creates a mask based on the reprojected shapefile for stats at property level.
    """
    obs.propertyMask = "lztmre_%s_mask_wholeproperty.img"  % obs.site
    obs.propertyMask = os.path.join(obs.outFolder, obs.propertyMask)

    shutil.copyfile(obs.window,obs.propertyMask)

    cmd1 = """
    gdal_translate -of HFA \
     -b 1 \
     -projwin %f %f %f %f  %s %s
     """ % ( obs.extent[0], obs.extent[3], obs.extent[1], obs.extent[2],
             obs.templateFile, obs.propertyMask)
    os.system(cmd1)

    cmd2 = """
    gdal_rasterize -i -burn 0 -l %s %s %s
    """ % (obs.layer.GetName(), obs.newShp, obs.propertyMask)
    os.system(cmd2)

    ## now make the mask
    cmd3 = """
    gdal_rasterize -burn 1 -l %s %s %s
    """ % (obs.layer.GetName(), obs.newShp, obs.propertyMask)
    os.system(cmd3)


    return obs


######################################################################################
# RIOS FUNCTIONS
#####################################################################################


def getPercentiles(info, inputs, outputs, otherArgs):
    """
    Given a layer of groundcover, we allocate a rank
    or percentile per stratum.
    """
    (nBands, nRows, nCols) = inputs.sfc.shape
    outputs.outdata = np.zeros( (1,nRows, nCols), 'uint8')  ## hopefully fewer than 256 strata
    percentiles = np.ones( (nRows, nCols)) * 0
    bare = inputs.sfc[0]
    maskvals = np.unique(inputs.mask)

    for strata in maskvals:
        if strata > 0:
            validpixels = (bare > 99) & (bare < 201) & (inputs.mask[0]==strata)
            if validpixels.any():
                thisdata = 200 - bare[validpixels]
                pct = np.percentile(thisdata, list(np.arange(100)))
                pctatpoint = np.interp(thisdata, pct, np.arange(100))
                ## this can be between 0 and 100,
                percentiles[validpixels] = pctatpoint + 100
    ## now scale appropriately
    outputs.outdata[0] = percentiles.astype('uint8')



def percentilesRoutine(mask,sfcImage,outputImage,obs):

    """
    Set up files for rios function
    """
    infiles = applier.FilenameAssociations()
    infiles.mask = mask
    infiles.sfc = sfcImage
    outfiles = applier.FilenameAssociations()
    outfiles.outdata = outputImage
    controls = applier.ApplierControls()
    controls.setReferenceImage(mask)
    controls.setFootprintType(applier.INTERSECTION)
    ## and force it to not block the image up
    dS = gdal.Open(infiles.mask)
    controls.setWindowXsize(dS.RasterXSize)
    controls.setWindowYsize(dS.RasterYSize)

    otherArgs = applier.OtherInputs()
    dS = None
    applier.apply(getPercentiles, infiles, outfiles, otherArgs,
                   controls=controls)

    obs.percentileImages.append(outputImage)




def summarisePercent(info, inputs, outputs, otherArgs):
    """
    summarise layers
    """

    (nRows, nCols) = inputs.images[0].shape[1:]
    outimage = np.zeros( (1, nRows, nCols), 'uint8')
    nimages = len(inputs.images)
    dataStack = np.zeros( (nimages, nRows, nCols))
    for i in range(nimages):
        dataStack[i] = inputs.images[i]
    ## remove na's
    dataStack[dataStack==0] = np.nan

    medimage = np.nanmedian(dataStack, axis=0)
    medimage[np.isnan(medimage)] = 0
    outimage[0] = medimage.astype('uint8')
    outputs.psum = outimage


def summariseRoutine(inImages,output):
    """
    Set up files for summarise rios function
    """
    infiles = applier.FilenameAssociations()
    infiles.images = inImages
    outfiles = applier.FilenameAssociations()
    outfiles.psum = output

    otherArgs = applier.OtherInputs()
    applier.apply(summarisePercent, infiles, outfiles, otherArgs)


def summariseDiff(info, inputs, outputs, otherArgs):
    """
    summarise layers
    """

    (nRows, nCols) = inputs.images[0].shape[1:]
    outimage = np.zeros( (1, nRows, nCols), 'int16')
    nimages = len(inputs.images)
    dataStack = np.zeros( (nimages, nRows, nCols))
    dataStack = inputs.images[0][0]
    ## Compute the difference
    anomolyImage = 0.0 + inputs.compare[0] - inputs.images[0][0]
    ## Add in the nodata values
    anomolyImage[np.where(inputs.mask[0]==0)] = 32767
    anomolyImage[np.where(inputs.compare[0]==0)] = 32767
    anomolyImage[np.where(inputs.images[0][0]==0)] = 32767

    outimage[0] = anomolyImage.astype('int16')
    outputs.mean = outimage


def diffRoutine(mask,inImages,compareImage,output):
    """
    Set up files for summarise rios function
    """
    infiles = applier.FilenameAssociations()
    infiles.images = inImages
    infiles.mask = mask
    infiles.compare = compareImage
    outfiles = applier.FilenameAssociations()
    outfiles.mean = output
    otherArgs = applier.OtherInputs()
    controls = applier.ApplierControls()
    controls.setReferenceImage(mask)
    controls.setFootprintType(applier.INTERSECTION)
    controls.setStatsIgnore(32767)
    applier.apply(summariseDiff, infiles, outfiles, otherArgs, controls=controls)

def summariseMean(info, inputs, outputs, otherArgs):
    """
    summarise layers
    """

    (nRows, nCols) = inputs.images[0].shape[1:]
    outimage = np.zeros( (1, nRows, nCols), 'int16')
    nimages = len(inputs.images)
    dataStack = np.zeros( (nimages, nRows, nCols))
    for i in range(nimages):
        dataStack[i] = inputs.images[i][0]
    ## remove na's
    dataStack[dataStack==0] = np.nan
    ## Compute the Mean anomoly
    anomolyImage = inputs.compare[0] - np.nanmean(dataStack, axis=0)
    ## Add in the nodata values
    anomolyImage[np.isnan(anomolyImage)] = 32767
    anomolyImage[np.where(inputs.mask[0]==0)] = 32767
    anomolyImage[np.where(inputs.compare[0]==0)] = 32767

    outimage[0] = anomolyImage.astype('int16')
    outputs.mean = outimage


def meanRoutine(mask,inImages,compareImage,output):
    """
    Set up files for summarise rios function
    """
    infiles = applier.FilenameAssociations()
    infiles.images = inImages
    infiles.mask = mask
    infiles.compare = compareImage
    outfiles = applier.FilenameAssociations()
    outfiles.mean = output
    otherArgs = applier.OtherInputs()
    controls = applier.ApplierControls()
    controls.setReferenceImage(mask)
    controls.setFootprintType(applier.INTERSECTION)
    controls.setStatsIgnore(32767)
    applier.apply(summariseMean, infiles, outfiles, otherArgs, controls=controls)

def summariseRank(info, inputs, outputs, otherArgs):
    """
    summarise layers
    """

    (nRows, nCols) = inputs.images[0].shape[1:]
    outimage = np.zeros( (1, nRows, nCols), 'uint8')
    nimages = len(inputs.images)
    dataStack = np.zeros( (nimages, nRows, nCols))
    for i in range(nimages):
        dataStack[i] = inputs.images[i][0]
    ## remove na's
    dataStack[dataStack==0] = np.nan
    comparisonImage = inputs.compare[0].astype('float32')
    comparisonImage[comparisonImage==0] = np.nan
    ## Compute the rank of the comparison image 
    rankImage = np.zeros([nRows, nCols])
    for i in range(nRows):
        for j in range(nCols):
            rankImage[i,j] = (percentileofscore(dataStack[:,i,j],comparisonImage[i,j]))
    ## Add in the nodata values
    rankImage[np.isnan(rankImage)] = 255
    rankImage[np.where(inputs.mask[0]==0)] = 255
    outimage[0] = rankImage.astype('uint8')
    outputs.mean = outimage

def rankRoutine(mask,inImages,compareImage,output):
    """
    Set up files for summarise rios function
    """
    infiles = applier.FilenameAssociations()
    infiles.images = inImages
    infiles.mask = mask
    infiles.compare = compareImage
    outfiles = applier.FilenameAssociations()
    outfiles.mean = output
    otherArgs = applier.OtherInputs()
    controls = applier.ApplierControls()
    controls.setReferenceImage(mask)
    controls.setFootprintType(applier.INTERSECTION)
    controls.setStatsIgnore(255)
    applier.apply(summariseRank, infiles, outfiles, otherArgs, controls=controls)



######################################################################################
# Main Driver function
#####################################################################################



class Observation(object):
    def __init__(self,json,outFolder, strata,site,template):
        self.json = json
        self.outFolder = outFolder
        self.strata = strata
        self.site = site
        self.templateFile = template
        self.percentileImages = []

def createStats(site,strata,startDate,endDate,compareDate,atype,seasonal,jsonPoly):
    
    TEMPLATE = None
    dixRange = []
    dixImages = find('*dixa2*.vrt', IMGDIR)

    for dix in dixImages:
        if TEMPLATE==None:TEMPLATE=dix
        dixElem = os.path.basename(dix).split('_')
        dixStart = dixElem[2][1:7]
        dixEnd = dixElem[2][7:13]
        

        if endDate != None and (int(dixStart) >= int(startDate)  and  int(dixEnd) <= int(endDate)):
            dixRange.append(dix)

        elif endDate != None and (int(dixStart) >= int(startDate)  and  len(dixRange) < 1):
            dixRange.append(dix)

        elif endDate == None and int(dixStart) >= int(startDate):
            dixRange.append(dix)


        # Determine the comparison image
        if compareDate != None and (int(dixStart) <= int(compareDate)  and  int(dixEnd) >= int(compareDate)):
            compareImage = dix

    # Make the comparison image the last image
    if compareDate == None:
            compareImage = dixRange[-1]

    # Seasonal Selection
    if seasonal:
        seasonalRange=[]
        # Determine the starting month of the reference Image
        seasonMonth = os.path.basename(compareImage).split('_')[2][5:7]

        for dixImage in dixRange:
            # Check if the dix image is the correct season
            if seasonMonth == os.path.basename(dixImage).split('_')[2][5:7]:
                # Add it to the output list
                seasonalRange.append(dixImage)
        if len(seasonalRange) > 0:
            # Only replace the images with the seasonal ones if there is a season in the comparison
            dixRange = seasonalRange
            
    # Create the output folder. Timestamp may noy be unique...
    timestamp = str(int(time.time()))
    outFolder = os.path.join(OUTPUT_DIR,'s'+timestamp)
    os.makedirs(outFolder) if not os.path.exists(outFolder) else None

    # Buils the initial Object
    obs = Observation(jsonPoly,outFolder, strata,site,TEMPLATE)
    obs = writeJson(obs)
    obs = reprojectShp(obs)
    obs = getExtent(obs)
    obs = getWindow(obs)
    obs = getPropertyMask(obs)

    for analysis in atype.split(","):
        if analysis == 'property':
            mask = obs.propertyMask
            for dix in dixRange:
                outimg = os.path.basename(dix).replace('dix','percent_%s' % analysis).replace('aus',obs.site).replace('.vrt','.img')
                outimg = os.path.join(outFolder,outimg)
                percentilesRoutine(mask,dix,outimg,obs)
            summariseRoutine(obs.percentileImages, os.path.join(outFolder, 'lztmre_%s_%s.img' % (obs.site,analysis)))

        elif analysis == 'paddock':
            obs = getPaddockMask(obs)
            mask = obs.paddockMask
            for dix in dixRange:
                outimg = os.path.basename(dix).replace('dix','percent_%s' % analysis).replace('aus',obs.site).replace('.vrt','.img')
                outimg = os.path.join(outFolder,outimg)
                percentilesRoutine(mask,dix,outimg,obs)
            summariseRoutine(obs.percentileImages, os.path.join(outFolder, 'lztmre_%s_%s.img' % (obs.site,analysis)))

        elif analysis == 'mean':
            meanRoutine(obs.propertyMask,dixRange,compareImage,os.path.join(outFolder, 'lztmre_%s_%s.img' % (obs.site,analysis)))

        elif analysis == 'rank':
            rankRoutine(obs.propertyMask,dixRange,compareImage,os.path.join(outFolder, 'lztmre_%s_%s.img' % (obs.site,analysis)))

        elif analysis == 'diff':
            diffRoutine(obs.propertyMask,dixRange,compareImage,os.path.join(outFolder, 'lztmre_%s_psummarya2_%s.img' % (obs.site,analysis)))



    return obs


 

######################################################################################
# Crossdomain rewriter
#####################################################################################

def crossdomain(origin=None, methods=None, headers=None,
                max_age=21600, attach_to_all=True,
                automatic_options=True):
    if methods is not None:
        methods = ', '.join(sorted(x.upper() for x in methods))
    if headers is not None and not isinstance(headers, basestring):
        headers = ', '.join(x.upper() for x in headers)
    if not isinstance(origin, basestring):
        origin = ', '.join(origin)
    if isinstance(max_age, datetime.timedelta):
        max_age = max_age.total_seconds()

    def get_methods():
        if methods is not None:
            return methods

        options_resp = current_app.make_default_options_response()
        return options_resp.headers['allow']

    def decorator(f):
        def wrapped_function(*args, **kwargs):
            if automatic_options and request.method == 'OPTIONS':
                resp = current_app.make_default_options_response()
            else:
                resp = make_response(f(*args, **kwargs))
            if not attach_to_all and request.method != 'OPTIONS':
                return resp

            h = resp.headers

            h['Access-Control-Allow-Origin'] = origin
            h['Access-Control-Allow-Methods'] = get_methods()
            h['Access-Control-Max-Age'] = str(max_age)
            if headers is not None:
                h['Access-Control-Allow-Headers'] = headers
            return resp

        f.provide_automatic_options = False
        return update_wrapper(wrapped_function, f)
    return decorator



######################################################################################
# Flask App
#####################################################################################

dbApp = Flask(__name__)
@dbApp.route("/getmap", methods = ['POST','GET'])
@crossdomain(origin='*')
def flaskFun():
    # Get the JSON data sent from the form. Sensible defaults set
    site = request.args.get('site',default='Spyglass', type=str)
    strata = request.args.get('strata',default='Id', type=str)
    startDate = request.args.get('startDate',default='198501', type=str)
    endDate = request.args.get('endDate',default='209901', type=str)
    compareDate = request.args.get('compareDate',default=None, type=str)
    jsonPoly = request.args.get('jsonPoly',default='{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},"geometry":{"type":"Polygon","coordinates":[[[143.6,-20.8],[143.7,-20.8],[143.7,-20.7],[143.6,-20.7],[143.6,-20.8]]]}}]}', type=str)
    atype = request.args.get('atype',default='property', type=str) # 'property','paddock','mean','rank','diff'
    seasonal = request.args.get('seasonal',default=False, type=bool)   
      
    # TODO: Check inputs for bad values
      
    obs = createStats(site,strata,startDate,endDate,compareDate,atype,seasonal,jsonPoly)

    #copy the img files to the shape file folder
    for f in os.listdir(obs.outFolder):
        if f.endswith('_paddock.img'): filename = f
        if f.endswith('_property.img'): filename = f
        if f.endswith('_mean.img'):  filename = f
        if f.endswith('_rank.img'): filename = f
        if f.endswith('_diff.img'): filename = f

    # Cleanup
    # shutil.rmtree(OUTPUT_FOLDER)

    return send_from_directory(obs.outFolder, filename, as_attachment=True)



######################################################################################
# MAIN
#####################################################################################

if __name__ == "__main__":
    dbApp.run(host="0.0.0.0",port=int("8080"),debug=True)







