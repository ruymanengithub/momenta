#! /usr/bin/env python

# IMPORT STUFF
import numpy as num
from pdb import set_trace as stop
import pyfits
from time import time
import os
from momsource import momsource
# END IMPORT

isthere = os.path.exists

def getCLstats(vector):
    # IMPORT STUFF
    from numpy import mean,median,std
    from CommonTools.maths import percentiles
    # END IMPORT
    ok = num.where(vector != -99.0)
    vector = vector[ok]
    if len(vector >1):
        max = vector.max()
        min = vector.min()
        medianval = float(median(vector))
        meanval = float(mean(vector))
        sigma = float(std(vector))
        p25 = percentiles(vector,25)
        p75 = percentiles(vector,75)
    else:
        max = -99.0 ; min = -99.0 ; medianval = -99.0 ; meanval = -99.0
        sigma = -99.0 ; p25 = -99.0 ; p75 = -99.0
    return max,min,medianval,meanval,sigma,p25,p75

def getRstats(radii,fluxes):
    """Statistics on the radial distribution of the clumps"""
    # IMPORT STUFF
    from scipy.stats import median
    from CommonTools.maths import wmean
    # END IMPORT
    maxRcl = radii.max()
    minRcl = radii.min()
    FweightRcl = float(wmean(radii,fluxes))
    FmaxRcl = radii[fluxes.argmax()]
    medianRcl = float(median(radii))
    return maxRcl, minRcl, FweightRcl, FmaxRcl, medianRcl
    
def runSextractor(image1,image2,photfile,sexcat,segimage,\
    sexvalues,sexpath='.'):
    """Executes sextractor"""
    # IMPORT STUFF
    from clumps_lib2 import ClumpsPhotfile
    # END IMPORT
    
    cwd = os.getcwd() # current working directory
    if sexpath != '.' : os.chdir(sexpath)
    
    ClumpsPhotfile(photfile,segimage,sexvalues)
    
    os.system('sex %s %s -c %s -CATALOG_NAME %s' %\
    (image1,image2,photfile,sexcat))
    
    if sexpath != '.' : os.chdir(cwd)
    
def ClumpsPhotfile(photfile,segimagen,sexvalues):
    """Makes a sextractor photfile to be used by 'DetectClumps'."""
    # IMPORT STUFF
    import string
    from clumps_verbose import Lphotpars
    # END IMPORT
    
    if 'CHECKIMAGE_TYPE' in sexvalues:
       if not isinstance(sexvalues['CHECKIMAGE_TYPE'],list):
           sexvalues['CHECKIMAGE_TYPE'] = [sexvalues['CHECKIMAGE_TYPE']]
           sexvalues['CHECKIMAGE_NAME'] = [sexvalues['CHECKIMAGE_NAME']]
       sexvalues['CHECKIMAGE_TYPE'].append('SEGMENTATION')
       sexvalues['CHECKIMAGE_NAME'].append(segimagen) 
       # always with segmentation...
    else:
       sexvalues['CHECKIMAGE_TYPE'] = ['SEGMENTATION']
       sexvalues['CHECKIMAGE_NAME'] = segimagen
    
    strsexvalues = {}
    
    for key in sexvalues:
        if isinstance(sexvalues[key],list):
            item = tuple([string.strip('%s'%issue) for issue in sexvalues[key]])
            strsexvalues[key] = ('%s,'*len(item))[0:-1] % item
        else:
            strsexvalues[key] = '%s'% string.strip('%s'% sexvalues[key])
    
    asset = []
    for key in Lphotpars:
        if key in strsexvalues:
            if 'none' in string.lower(strsexvalues[key]):
                asset.append('#')
            else : asset.append('')
            asset.append(strsexvalues[key])
        else:
            asset.append('#') ; assent.append('')
    asset = tuple(asset)
    
    text = \
"""
# Default configuration file for SExtractor V1.2b14 - > 2.0
# EB 23/07/98
# (*) indicates parameters which can be omitted from this config file.
#-------------------------------- Catalog ------------------------------------
    
%sCATALOG_TYPE 	      %s	# "NONE","ASCII_HEAD","ASCII","FITS_1.0"
                                            # or "FITS_LDAC"
    
%sCATALOG_NAME       %s
%sPARAMETERS_NAME %s
#------------------------------- Extraction ----------------------------------
    
%sDETECT_TYPE           %s		# "CCD" or "PHOTO" (*)
%sFLAG_TYPE               %s
%sFLAG_IMAGE            %s
%sFITS_UNSIGNED       %s
%sDETECT_MINAREA	 %s		# minimum number of pixels above threshold
%sDETECT_THRESH	 %s	    # <sigmas> or <threshold>,<ZP> in     mag.arcsec-2
%sANALYSIS_THRESH	  %s      # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
%sTHRESH_TYPE          %s
%sFILTER		                  %s		# apply filter for detection ("Y" or "N")?
%sFILTER_NAME	          %s	    # name of the file containing the filter
    
%sDEBLEND_NTHRESH	%s		# Number of deblending sub-thresholds
%sDEBLEND_MINCONT	%s	    # Minimum contrast parameter for deblending
    
%sCLEAN		                    %s		# Clean spurious detections? (Y or N)?
%sCLEAN_PARAM	        %s		# Cleaning efficiency
    
%sMASK_TYPE	                %s	    # type of detection MASKing: can be one of
                                        # "NONE", "BLANK" or "CORRECT"
    
#------------------------------ Photometry -----------------------------------
    
%sPHOT_APERTURES      %s
%sPHOT_AUTOPARAMS	%s 	# MAG_AUTO parameters: <Kron_fact>,<min_radius>
%sPHOT_FLUXFRAC        %s
%sMAG_ZEROPOINT       %s
%sGAIN                           %s
%sSATUR_LEVEL	            %s		# level (in ADUs) at which arises saturation
%sMAG_GAMMA             %s         # gamma of emulsion (for photographic scans)
%sPIXEL_SCALE              %s		    # size of pixel in arcsec (0=use FITS WCS info).
    
#------------------------- Star/Galaxy Separation ----------------------------
    
%sSEEING_FWHM	       %s	# stellar FWHM in arcsec
%sSTARNNW_NAME	   %s	# Neural-Network_Weight table filename
    
#------------------------------ Background -----------------------------------
    
%sBACK_SIZE	              %s  # Background mesh: <size> or <width>,<height>
%sBACK_FILTERSIZE     %s  # Background filter: <size> or <width>,<height>
%sBACKPHOTO_TYPE    %s   # can be "GLOBAL" or "LOCAL" (*)
%sBACKPHOTO_THICK  %s  # thickness of the background LOCAL annulus (*)
%sBACK_TYPE                %s
%sBACK_VALUE             %s
    
#------------------------------ Check Image ----------------------------------
    
%sCHECKIMAGE_TYPE   %s    # can be one of "NONE", "BACKGROUND",
                                # "MINIBACKGROUND", "-BACKGROUND", "OBJECTS",
                                # "-OBJECTS", "SEGMENTATION", "APERTURES",
                                # or "FILTERED" (*)
%sCHECKIMAGE_NAME  %s
    
#--------------------- Memory (change with caution!) -------------------------
    
%sMEMORY_OBJSTACK	%s   		# number of objects in stack
%sMEMORY_PIXSTACK	%s	        # number of pixels in stack
%sMEMORY_BUFSIZE      %s    		# number of lines in buffer
    
#----------------------------- Miscellaneous ---------------------------------
    
%sVERBOSE_TYPE	        %s		    # can be "QUIET", "NORMAL" or "FULL" (*)
%sINTERP_MAXXLAG      %s
%sINTERP_MAXYLAG      %s
%sINTERP_TYPE              %s
    
    
""" % asset
    
    f = open(photfile,'w')
    f.write(text)
    f.close()
    
    return None


def ParseCatClumps(sexcat):
    """Function that parses a photometric catalog in search of clumps."""
    # IMPORT STUFF
    from Moments.momsource import momcat
    from clumps_verbose import Lclsexcols, Dclsexformats
    from copy import copy
    # END IMPORT
    
    if not isthere(sexcat) : return None
    
    father_cat_name = sexcat
    execpars = {'father_cat_name':father_cat_name,'sexpars':Lclsexcols}
    catcl = momcat(execpars)
    try: catcl.read_catsex(Dclsexformats)
    except RuntimeError: return None
    
    data = {}
    for key in Lclsexcols : data[key] = copy(catcl[key])
    
    return data
    
def profileofModel(self):
    """Retrieves a radial profile of the 2D model of the object provided by
    FlxAz"""
    # IMPORT STUFF
    from copy import copy
    # END IMPORT
    
    execpars = copy(self.execpars)
    dummy = momsource(execpars)

    # INPUTS
    dummy['STAMP'] = self['FlxAz'].copy()
    # image masking
    dummy['MASKOTHER'] = self['MASKOTHER'].copy()
    dummy['MASK'] = self['MASK'].copy()
    # geometrical parameters
    parlist = ['X_IMAGE','Y_IMAGE','MXMIN_IMAGE','MYMIN_IMAGE',
    'ELLIPTICITY','THETA_IMAGE','BOXY','BACKGROUND']
    for par in parlist : dummy[par] = self[par]
    
    if dummy.execpars['doFineRadial'][0]==0:
        dummy.radial_v1(dograph=False)
    else:
        dummy.radial_v3(dograph=False)
    
    self['M_MOD_RADIAL'] = copy(dummy['M_RADIAL'])
    del dummy
    return None
    
