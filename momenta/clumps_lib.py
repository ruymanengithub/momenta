#! /usr/bin/env python

"""Set of programs which find and make statistics on clumps in an image 
making use of Sextractor's segmentation tool"""

# IMPORT STUFF
import numpy as num
from pdb import set_trace as stop
import pyfits
from time import time
import os
from momsource import momsource
from copy import copy
# END IMPORT

isthere = os.path.exists

# DATA

from clumps_verbose import Lclsexcols,Dclsexformats,\
Lphotpars,Ddefaultphotpars

def Clumpscore(self):
    """Method which retrieves statistics on clumpy emission of an image. It 
    may be possible to substract first a model of the radial distribution of flux. 
    Also, fotometry using clumps detected in another band can be done. 
    This can be (as well) applied to a map of mass instead a map of intensity."""
    # IMPORT STUFF
    import string
    from Moments.easygui import multenterbox
    from trunc_lib import show_on_ds9
    # END IMPORT
    # INPUTS
    blanknegs = False
    clsubradial = self.execpars['clsubradial'][0] == 1
    cldetectother = self.execpars['cldetectother'][0] == 1
    object = self['name']
    sky = self['BACKGROUND']
    try: tag = self.execpars['tag'][0]
    except : tag = ''
    interactive = self.execpars['interactive'][0] == 1
    dopdf = 'pdf' in self
    maskoffinside = self.execpars['maskoffinside'][0] == 1
    doMEF = self.execpars['doMEF'][0] == 1
    clmaxradius = self.execpars['clmaxradius'][0] == 1
    petrofactor = self.execpars['petrofactor'][0]
    # END INPUTS (MORE BELOW!)
    
    done = False
    
    while not done:
        oriclimg = self['STAMP'].copy()
        try: mask = self['MASKOTHER'].copy()
        except AttributeError : mask = self['MASKOTHER']
        try: extmask = self['EXTMASK'].copy()
        except AttributeError : extmask = self['EXTMASK']
        oriclimg -= sky
        if mask is -1: mask = num.zeros(shape=climg.getshape(),type='Int8')
        if extmask is -1: 
            extmask = num.zeros(shape=climg.getshape(),type='Int8')
        if not maskoffinside:
            oriclimg[num.where(mask > 0)] = 0.
        else:
            pass
        oriclimg[num.where(extmask > 0)] = 0.
        
        if clsubradial:
            FlxAz = self.FlxAz() # creates 2D model of galaxy from radial profile
            if self.execpars['doRadialModel'] == 1:
                self.profileofModel() # makes a 1D profile of created 2D model
            
            climg = oriclimg - FlxAz
            if blanknegs:
                neg = num.where(climg<0.)
                climg[neg] = 0.
        else: climg = oriclimg.copy()
        
        if cldetectother:
            clcat_in = self['CL_cat_in']
            segimg_in = self['CL_seg_in']
            self.MeasureClumps(climg,clcat_in,segimg_in,mask=-1)
        else:
            self.DetectClumps(climg,mask=-1,cleanafter=False)
        
        self.AnalClumps()
        
        # PDF creation (opt)
        
        if dopdf or doMEF:
            
            # create figures : images and plots
            # It also produces a MEF... this is new...
            
            figlist = []
            imglist = []
            oriclimgn = 'CL_%s%s_IMG1.eps' % (object,tag)
            azimgn = 'CL_%s%s_IMG2.eps' % (object,tag)
            climgn = 'CL_%s%s_IMG3.eps' % (object,tag)
            clsegn = 'CL_%s%s_IMG4.eps' % (object,tag)
            
            # Input image
            if dopdf:
                self.ClDisplay(oriclimg,oriclimgn)
                figlist.append(oriclimgn)
            imglist.append(oriclimg)
            # Segmented clumps
            oclseg = self['CL_segdata'].copy()
            clseg = oclseg * 0.
            clseg[num.where(oclseg > 0)] = 1
            if dopdf:
                self.ClDisplay(clseg.astype('Float32'),clsegn)
                figlist.append(clsegn)
            imglist.append(oclseg.astype('Float32'))
            
            # Radial model (opt)
            if clsubradial:
                if dopdf:
                    self.ClDisplay(FlxAz,azimgn)
                    figlist.append(azimgn)
                imglist.append(FlxAz)
            else:
                imglist.append((oriclimg * 0.).copy())    
            if dopdf:
                self.ClDisplay(climg,climgn)
                figlist.append(climgn)
            imglist.append(climg)
            
            if dopdf:
                plotlist = []
                plintensn = 'CL_%s%s_PLOT1.eps' % (object,tag)
                plfluxn = 'CL_%s%s_PLOT2.eps' % (object,tag)
                plareasn = 'CL_%s%s_PLOT3.eps' % (object,tag)
                plotlist = [plintensn,plfluxn,plareasn]
                self.ClRadialPlot('intens',plintensn)
                self.ClRadialPlot('clfluxes',plfluxn)
                self.ClRadialPlot('clareas',plareasn)
                for plotn in plotlist:
                    if not isthere(plotn): plotlist[plotlist.index(plotn)] = None
            
            if dopdf: self.ClLatex(figlist,plotlist)
            if doMEF:
                MEF = 'CL_%s%s_MEF.fits' % (object,tag)
                try: 
                    pyfits.writeto(MEF,num.array(imglist))
                    self['MEF'] = MEF
                except :
                    if isthere(MEF) : os.system('rm %s' % MEF)
                    for imef in range(len(imglist)):
                        pyfits.writeto('%s_%i' % (MEF,imef+1),num.array(imglist[imef]))
                    self['MEF'] = '%s_#' % MEF
                del imglist
        
        if interactive:
            # show images on ds9
            tempN = 'DS9_%s%s'%(object,tag) + '_%s.fits'
            
            imagesds9 = [oriclimg,clseg.astype('Float32')]
            namesds9 = [tempN%'ORI',tempN%'SEG']
            if clsubradial: 
                imagesds9 = imagesds9 + [FlxAz,climg]
                namesds9 = namesds9 + [tempN%'AZI',tempN%'CLU']
            for imname in namesds9:
                if isthere(imname) : os.system('rm %s' % imname)
            myds9 = show_on_ds9(imagesds9,namesds9)
            
            # show pdf
            if 'pdf' in self: os.system('xpdf %s &' % self['pdf'])
            
            try: CLTHRESH = copy(self['CLTHRESH'])
            except KeyError : CLTHRESH = copy(self.execpars['CLTHRESH'])
            try: CLMINAREA = copy(self['CLMINAREA'])
            except KeyError : CLMINAREA = copy(self['CLMINAREA'])
            
            GUInames = ['CLTHRESH','CLMINAREA','clsubradial']
            GUIvalues = [CLTHRESH,CLMINAREA,clsubradial]
            
            message = 'Wanna change values?' ; GUItitle = 'CLUMPS'
            updated = multenterbox(message,GUItitle,GUInames,GUIvalues)
            if updated == None : updated = GUIvalues
            for i in range(len(GUInames)): print GUInames[i],updated[i]
            print '\n'
            try:
                self['CLTHRESH'] = float(updated[0])
                self['CLMINAREA'] = float(updated[1])
                clsubradial = int(updated[2]) == 1
            except : pass
            
            done = string.lower(raw_input('Done with the object? yes = y '))\
            == 'y'
            
            for imname in namesds9:
                if isthere(imname) : os.system('rm %s' % imname)
            
            isopen = myds9.isOpen()
            if isopen : myds9.xpaset('exit')
            del myds9
            
        else: done = True
    
    return None

def MeasureClumps(self,climage,clcat_in,segimg_in,mask=-1):
    """Measure Clumps on an image, given a segmentation map and a 
    photometric catalog. Intended to be used when the detection of the clumps 
    is done on another image."""
    # IMPORT STUFF
    from Moments.momsource import momcat
    from clumps_lib2 import ParseCatClumps
    from copy import copy
    # END IMPORT
    # INPUTS
    object = self['name']
    try: tag = self.execpars['tag'][0]
    except KeyError: tag = ''
    clmaxradius = self.execpars['clmaxradius'][0] == 1
    clglobal = self.execpars['clglobal'][0] == 1
    # END INPUTS
    
    # Naming files
    
    clcat_out = 'CL_%s%s.cat' % (object,tag)
    segimg_out = 'CL_%s%s_seg.fits' % (object,tag)
    
    # Masking
    
    if mask is not -1 : 
        climage[num.where(mask > 0.)] = 0.
    
    radialmask = self.getradialmask()
    areaobject = float(len(num.where(self['SEXMASK']==0)[0]))
    ellip = self['ELLIPTICITY']
    maxradius = (areaobject / (num.pi * (1.-ellip)))**0.5
    petrofactor = self.execpars['petrofactor'][0]
    if petrofactor > 0.: maxradius *= petrofactor
    
    if clmaxradius:
        climage[num.where(radialmask>maxradius)] = 0.
    
    # 'NUMBER','X_IMAGE','Y_IMAGE','A_IMAGE', 
    # 'B_IMAGE', 'THETA_IMAGE','ELLIPTICITY', 'FLUX_ISO','ISOAREA_IMAGE'
    
    # Read out input catalog and prepare output catalog
    
    execp = {'father_cat_name':clcat_in,'sexpars':Lclsexcols}
    try:
        data_in = momcat(execp)
        data_in.read_catsex(Dclsexformats)
    except:
        self['CL_cat'] = clcat_out
        self['CL_catdata'] = {}
        self['CL_seg'] = 'None'
        self['CL_segdata'] = climage * 0.
        return None
    
    data_out = copy(data_in) #
    data_out.execpars['outfile_name'] = clcat_out
    data_out.execpars['to_output'] = Lclsexcols
    
    # Loading segmentation map and dumping to output sementation map
    
    segclumps = pyfits.getdata(segimg_in).astype('Int32')
    pyfits.writeto(segimg_out,segclumps.astype('Int32'))
    
    # Updating 'FLUX_ISO'
    
    FLUX_ISO = []
    
    for i in range(len(data_in['NUMBER'])):
        NUMBER = data_in['NUMBER'][i]
        ix = num.where(segclumps == NUMBER)
        FLUX_ISO.append(num.sum(climage[ix]))
    
    FLUX_ISO = num.array(FLUX_ISO)
    FLUX_ISO[num.where(FLUX_ISO<=0.)] = -99.0
    
    data_out['FLUX_ISO'] = FLUX_ISO.copy()
    
    # Dumping output catalog
    comments = {}
    data_out.dumpcat(Lclsexcols,comments,Dclsexformats)
    
    self['CL_cat'] = clcat_out
    self['CL_catdata'] = ParseCatClumps(clcat_out)
    self['CL_seg'] = segimg_out
    self['CL_segdata'] = segclumps.copy()
    
    if clglobal:
        self.getclglobal(radialmask,maxradius,segclumps,climage)
    
    return None
    
def DetectClumps(self,climage,mask=-1,cleanafter=False):
    """Finds clumps in an image making use of Sextractor's segmentation
    tool. Key parameters are the THRESHOLD and MINAREA for detection."""
    # IMPORT STUFF
    from clumps_lib2 import runSextractor, ParseCatClumps
    from numpy import random_array
    # END IMPORT
    # INPUTS
    try:
        threshold = self['CLTHRESH'] # counts
        if threshold == 0.: raise KeyError
    except KeyError : threshold = self.execpars['CLTHRESH'][0]
    try:
        minarea = self['CLMINAREA'] # pixels
        if minarea == 0 : raise KeyError
    except KeyError: minarea = self.execpars['CLMINAREA'][0]
    object = self['name']
    try: tag = self.execpars['tag'][0]
    except KeyError : tag = ''
    clpath = os.getenv('CL_DIR')
    clmaxradius = self.execpars['clmaxradius'][0] == 1
    clglobal = self.execpars['clglobal'][0] == 1
    # END INPUTS
    
    if 'CLTHRESH' in self : self['CLTHRESH'] = threshold
    if 'CLMINAREA' in self : self['CLMINAREA'] = minarea
    
    # Masking
    
    if mask is not -1 : climage[num.where(mask > 0)] = 0.
    
    radialmask = self.getradialmask()
    areaobject = float(len(num.where(self['SEXMASK']==0)[0]))
    ellip = self['ELLIPTICITY']
    maxradius = (areaobject / (num.pi * (1.-ellip)))**0.5
    petrofactor = self.execpars['petrofactor'][0]
    if petrofactor > 0.: maxradius *= petrofactor
    
    if clmaxradius:
        climage[num.where(radialmask>maxradius)] = 0.

    climagen = 'tmpDetectClumps_%i.fits' % random_array.randint(1,10000)
    pyfits.writeto(climagen,climage.astype('Float32'))
    
    segimagen = 'CL_%s%s_seg.fits' % (object,tag)
    photfile = 'CL_%s%s_phot.txt' % (object,tag)
    sexcat = 'CL_%s%s.cat' % (object,tag)
    sexvalues = {}
    
    # Sextractor parameters : sexvalues
    
    # no filtering, deblending true, background computed or substracted
    # in photometry.
    
    sexvalues = {'CATALOG_TYPE':'ASCII_HEAD',
    'PARAMETERS_NAME':'%s/clsexpars.txt' % clpath,'DETECT_TYPE':'CCD',
    'DETECT_MINAREA':minarea,'DETECT_THRESH':threshold,
    'ANALYSIS_THRESH':threshold,'THRESH_TYPE':'ABSOLUTE',
    'FILTER':'N','DEBLEND_NTHRESH':2,'DEBLEND_MINCONT':0.5,
    'CLEAN':'N','CLEAN_PARAM':1.0,'MASK_TYPE': 'None',
    'PHOT_AUTOPARAMS':[2.5,3.5],'MAG_ZEROPOINT':25,'GAIN':1.0,
    'BACK_TYPE': 'MANUAL','BACK_VALUE':[0.0,0.0],
    'BACK_SIZE':'NONE','BACK_FILTERSIZE':'NONE',
    'BACKPHOTO_TYPE': 'GLOBAL',
    'MEMORY_OBJSTACK' : 50000,'MEMORY_PIXSTACK':10000000,
    'MEMORY_BUFSIZE':32000,
    'VERBOSE_TYPE': 'NORMAL','INTERP_MAXXLAG':2,
    'INTERP_MAXYLAG':2,'INTERP_TYPE':'ALL'}
    
    for key in Lphotpars : 
        if key not in sexvalues and 'CHECKIMAGE' not in key: 
            sexvalues[key] = 'None'
    
    # Run sextractor
    
    runSextractor(climagen,climagen,photfile,sexcat,segimagen,sexvalues,
    sexpath='.')
    
    # Did it go right?
    
    if (not isthere(segimagen)) and (not isthere(sexcat)):
        self['CL_cat'] = None
        return None
    
    # Load segmentation
    
    self['CL_seg'] = segimagen
    segclumps = pyfits.getdata(segimagen).astype('Int32')
    self['CL_segdata'] = segclumps.copy()

    # Read and parse output catalog
    # save from catalog: NUMBER, ISOAREA_IMAGE, THETA_IMAGE, 
    # ELLIPTICITY, X_IMAGE, Y_IMAGE, FLUX_ISO, A_IMAGE, B_IMAGE
    
    dataclumps = ParseCatClumps(sexcat)
    
    if clglobal :
        self.getclglobal(radialmask,maxradius,segclumps,climage)
    
    self['CL_cat'] = sexcat
    self['CL_catdata'] = dataclumps
    
    os.system('rm %s %s' % (climagen,photfile))
    if cleanafter: os.system('rm %s %s' % (segimagen,sexcat))
    
    return None

    
def AnalClumps(self):
    """Analyzes Clumps found in an image. """
    # IMPORT STUFF
    from clumps_lib2 import getCLstats,getRstats
    # END IMPORT
    # INPUTS
    # END INPUTS
    # TO MEASURE:
    #   - number of distinct clumps
    #   - fluxes (max, min, median, mean,sigma,p25,p75), 
    #   - areas (max, min, median, mean,sigma,p25,p75)
    #   - Total flux in clumps.
    #   - Total area of clumps.
    #   - Radial profiles of clumps fluxes, and areas : save RADII
    
    data = self['CL_catdata']
    
    if not isinstance(data,dict) or len(data.keys()) == 0:
        Ncl = 0
##        toupdate= ['maxAcl','maxFcl','minAcl','minFcl',
##        'medianAcl','medianFcl','meanAcl','meanFcl',
##        'sigmaAcl','sigmaFcl','p25Acl','p25Fcl','p75Acl',
##        'p75Fcl','maxRcl','minRcl','FweightRcl','FmaxRcl',
##        'medianRcl','phiAcl','phiFcl']
##        for key in toupdate: self[key] = -99.0
    else:
        NUMBER = data['NUMBER']
        X_IMAGE = data['X_IMAGE'] ; Y_IMAGE = data['Y_IMAGE']
        XCENTER = self['X_IMAGE'] - self['MXMIN_IMAGE']
        YCENTER = self['Y_IMAGE'] - self['MYMIN_IMAGE']
        X = X_IMAGE - XCENTER
        Y = Y_IMAGE - YCENTER
        RADII = (X**2. + Y**2.)**0.5
        self['CL_catdata']['RADII'] = RADII
        
        FLUXES = data['FLUX_ISO']
        AREAS = data['ISOAREA_IMAGE']
        
        Fcl = num.sum(FLUXES)
        Acl = num.sum(AREAS)
        Ncl = len(NUMBER)
    
    if Ncl >= 2:
        
        maxFcl,minFcl,medianFcl,meanFcl,\
        sigmaFcl,p25Fcl,p75Fcl = getCLstats(FLUXES)
        
        maxAcl,minAcl,medianAcl,meanAcl,\
        sigmaAcl,p25Acl,p75Acl = getCLstats(AREAS)
        
        maxRcl, minRcl, FweightRcl, FmaxRcl, medianRcl = getRstats(RADII,FLUXES)
        
    if Ncl == 1:
        maxFcl = minFcl = medianFcl = meanFcl = p25Fcl = p75Fcl = FLUXES[0]
        sigmaFcl = 0.
        maxAcl = minAcl = medianAcl = meanAcl = p25Acl = p75Acl= AREAS[0]
        sigmaAcl = 0.
        maxRcl = minRcl = FweightRcl = FmaxRcl = medianRcl = RADII[0]
     
    if Ncl == 0:
        Fcl = -99.0 ; Acl = -99.0 ; Ncl = 0
        maxFcl = minFcl = medianFcl = meanFcl = \
        sigmaFcl = p25Fcl = p75Fcl = -99.0
        maxAcl = minAcl = medianAcl = meanAcl = \
        sigmaAcl = p25Acl = p75Acl = -99.0
        maxRcl = minRcl = FweightRcl = FmaxRcl = medianRcl = -99.0 
    
    self['maxAcl'] = maxAcl ; self['maxFcl'] = maxFcl ;
    self['minAcl'] = minAcl ; self['minFcl'] = minFcl ;
    self['medianAcl'] = medianAcl ; self['medianFcl'] = medianFcl ;
    self['meanAcl'] = meanAcl ; self['meanFcl'] = meanFcl ;
    self['sigmaAcl'] = sigmaAcl ; self['sigmaFcl'] = sigmaFcl ;
    self['p25Acl'] = p25Acl ; self['p25Fcl'] = p25Fcl ;
    self['p75Acl'] = p75Acl ; self['p75Fcl'] = p75Fcl ;
    self['maxRcl'] = maxRcl ;self['minRcl'] = minRcl
    self['FweightRcl'] = FweightRcl ; self['FmaxRcl'] = FmaxRcl
    self['medianRcl'] = medianRcl
    self['Fcl'] = Fcl ; self['Acl'] = Acl
    self['Ncl'] = Ncl
    
    # Filling factors (area and flux)
    
    if Acl > 0.:
        try: phiAcl = max(0., Acl / self['M_NPIX']) # M_NPIX may be -99.0 ?
        except KeyError : phiAcl = -99.0
    else : phiAcl = -99.0
    
    if Fcl > 0.:
        try: phiFcl = max(0.,Fcl / self['M_FLUX']) # M_FLUX may be -99.0 ?
        except KeyError : phiFcl = -99.0
    else: phiFcl = -99.0
    
    self['phiAcl'] = phiAcl ; self['phiFcl'] = phiFcl
    
    return None
    

def FlxAz(self):
    """Returns an azimuthally averaged model of a galaxy's Flux Distribution.
    It uses an Intensity Profile and geometric parameter (center,ellipticity,
    thetha), and a mask."""
    # IMPORT STUFF
    from Moments.ellipse import dist_superellipse
    from Moments.radial import RadialModel
    # END IMPORT
    
    # INPUTS
    img = self['STAMP'].copy()
    mask = self['MASKOTHER'].copy()
    radial = self['M_RADIAL'].copy()
    xcenter = self['X_IMAGE']  ; ycenter = self['Y_IMAGE']
    xcenter = xcenter - self['MXMIN_IMAGE']
    ycenter = ycenter - self['MYMIN_IMAGE']
    ellip = self['ELLIPTICITY'] 
    pa = self['THETA_IMAGE'] - 90. # Astronomical position angle.
    c = self['BOXY']
    areaobject = float(len(num.where(self['SEXMASK']==0)[0]))
    maxradius = (areaobject / (num.pi * (1.-ellip)))**0.5
    maxradius *= self.execpars['petrofactor'][0] 
    # END INPUTS

    radii = radial['radii']
    intens = radial['intens']
    goodr = num.where(radii <= maxradius)
    radii = radii[goodr] ; intens = intens[goodr]
    
    n = img.shape
    dims = n[-1::-1]
    q = 1. - ellip
    icenter = (xcenter,ycenter)
    FlxAz = RadialModel(radii,intens,dims,icenter,q,pa,c=0.)
    
    # FlxAz[num.where(mask)] = 0.
    negative = num.where(FlxAz<0.)
    FlxAz[negative] = 0.
    
    self['FlxAz'] = FlxAz.copy()
    
    return FlxAz    

def getradialmask(self):
    """Retrieves a radial mask for an object."""
    # IMPORT STUFF
    from ellipse import dist_superellipse
    # END IMPORT
    # INPUTS
    xcenter = self['X_IMAGE']  ; ycenter = self['Y_IMAGE']
    xcenter = xcenter - self['MXMIN_IMAGE']
    ycenter = ycenter - self['MYMIN_IMAGE']
    ellip = self['ELLIPTICITY'] 
    pa = self['THETA_IMAGE'] - 90. # Astronomical position angle.
    c = self['BOXY']
    n = self['STAMP'].shape
    # END INPUTS
    q = 1. - ellip
    ecenter = (ycenter,xcenter)
    return dist_superellipse(n,ecenter,q=q,pos_ang=pa,c=c)
    
def getclglobal(self,radialmask,maxradius,segclumps,climage):
        inix = num.where((radialmask <= maxradius) & (segclumps > 0.) & \
        (climage != 0.))
        outix = num.where((radialmask <= maxradius) & (segclumps == 0.) &\
        (climage != 0.))
        allix = num.where((radialmask <= maxradius) & (climage != 0.))
        try:
            FL_IN = num.sum(climage[inix])
            N_IN = len(inix[0])
        except: FL_IN = -99.0 ; N_IN = 0
        try: 
            FL_OUT = num.sum(climage[outix])
            N_OUT = len(outix[0])
        except : FL_OUT = -99.0 ; N_OUT = 0
        try: 
            FL_ALL = num.sum(climage[allix])
            N_ALL = len(allix[0])
        except : FL_ALL = -99.0 ; N_ALL = 0
        
        self['FL_IN'] = FL_IN ; self['FL_OUT'] = FL_OUT ; self['FL_ALL'] = FL_ALL
        self['N_IN'] = N_IN ; self['N_OUT'] = N_OUT ; self['N_ALL'] = N_ALL
        magzero = self.execpars['magzero'][0]
        for appendix in ['IN','OUT','ALL']:
            FL = self['FL_%s' % appendix]
            if FL>0.: self['MG_%s' % appendix] =  magzero - 2.5 * num.log10(FL)
            else: self['MG_%s' % appendix] = -99.0
        scale = self.execpars['scale'][0]
        for appendix in ['IN','OUT','ALL']:
            NP = self['N_%s' % appendix]
            if NP >0 : self['AR_%s' % appendix] =  NP * scale**2.
            else : self['AR_%s' % appendix] = -99.0
            
        return None
        
class Clumpssource(momsource):
    # IMPORT STUFF
    from clumps_lib import Clumpscore, MeasureClumps,DetectClumps,\
    AnalClumps,FlxAz,getradialmask,getclglobal
    from clumps_lib2 import profileofModel
    from clumps_lib3 import ClLatex,ClRadialPlot,ClDisplay
    # END IMPORT

