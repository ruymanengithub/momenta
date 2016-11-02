#! /usr/bin/env python

# IMPORT STUFF
import os
from pdb import set_trace as stop
# END IMPORT

def runSextractor(image1,image2,photfile,sexcat,segimage,sexvalues,sexpath='.'):
    """Executes sextractor."""
    # IMPORT STUFF
    from sexclumps import doPhotfile
    # END IMPORT
    
    cwd = os.getcwd() # current working directory
    if sexpath != '.' : os.chdir(sexpath)
    
    doPhotfile(photfile,segimage,sexvalues)
    
    os.system('mysex %s %s -c %s -CATALOG_NAME %s' %\
    (image1,image2,photfile,sexcat))
    
    if sexpath != '.' : os.chdir(cwd)

def doPhotfile(photfile,segimage,sexvalues={}):
    """Writes a photfile (.sex) for running sextractor."""
    # IMPORT STUFF
    import string
    # END IMPORT
    
    sexpars = ['CATALOG_TYPE','CATALOG_NAME','PARAMETERS_NAME','DETECT_TYPE','FLAG_TYPE','FLAG_IMAGE','FITS_UNSIGNED','DETECT_MINAREA','DETECT_THRESH','ANALYSIS_THRESH','THRESH_TYPE','FILTER','FILTER_NAME','DEBLEND_NTHRESH','DEBLEND_MINCONT','CLEAN','CLEAN_PARAM','MASK_TYPE','PHOT_APERTURES','PHOT_AUTOPARAMS','PHOT_FLUXFRAC','MAG_ZEROPOINT','GAIN','SATUR_LEVEL','MAG_GAMMA','PIXEL_SCALE','SEEING_FWHM','STARNNW_NAME','BACK_SIZE','BACK_FILTERSIZE','BACKPHOTO_TYPE','BACKPHOTO_THICK','BACK_TYPE','BACK_VALUE','CHECKIMAGE_TYPE','CHECKIMAGE_NAME','MEMORY_OBJSTACK','MEMORY_PIXSTACK','MEMORY_BUFSIZE','VERBOSE_TYPE','INTERP_MAXXLAG','INTERP_MAXYLAG','INTERP_TYPE']
    
    if 'CHECKIMAGE_TYPE' in sexvalues:
       if not isinstance(sexvalues['CHECKIMAGE_TYPE'],list):
           sexvalues['CHECKIMAGE_TYPE'] = [sexvalues['CHECKIMAGE_TYPE']]
           sexvalues['CHECKIMAGE_NAME'] = [sexvalues['CHECKIMAGE_NAME']]
       sexvalues['CHECKIMAGE_TYPE'].append('SEGMENTATION')
       sexvalues['CHECKIMAGE_NAME'].append(segimage)
    else:
       sexvalues['CHECKIMAGE_TYPE'] = ['SEGMENTATION']
       sexvalues['CHECKIMAGE_NAME'] = segimage
    
    strsexvalues = {}
    
    for key in sexpars:
        if isinstance(sexvalues[key],list):
            item = tuple([string.strip('%s'%issue) for issue in sexvalues[key]])
            strsexvalues[key] = ('%s,'*len(item))[0:-1] % item
        else:
            strsexvalues[key] = '%s'% string.strip('%s'% sexvalues[key])
    
    asset = []
    for key in sexpars:
        if 'None' in strsexvalues[key] or 'none' in strsexvalues[key] :
            asset.append('#')
        else : asset.append('')
        asset.append(strsexvalues[key])
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
    
