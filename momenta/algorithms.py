#! /usr/bin/env/ python

"""Set of functions that operate on a image + mask, executing some algorithm.

Functions: gini, ellipse, moments (1st and 2nd), axis_symmetry, 
asymmetry (A), clumpiness (S), concentration (C), peak, M20.
"""
# IMPORT STUFF
import numpy as num
from pdb import set_trace as stop
import sys
from pyraf import iraf
import pyfits
from time import time
from Moments.graphs import Convert, doStamp
# END IMPORT

def getsky(self):
    """Gets basic info about background."""
    sky_image = self['STAMP'].copy()
    skymask = self['SKYMASK'].copy()
    sky_image[num.where(skymask != 0)] = -10000.
    sky_median = get_stat(sky_image,'midpt',minimum=-9999.)
    sky_sigma = get_stat(sky_image,'stddev',minimum=-9999.)
    self['SKY_MEDIAN'] = sky_median
    self['SKY_SIGMA'] = sky_sigma
        
def ellipse(self):
    """Returns A, B, THETA elongation and ellipticity of distribution of pixels."""
    # IMPORT STUFF
    from algorithms import moments
    # END IMPORT
    
    sky = self['BACKGROUND']
    img = self['STAMP'].copy() - sky
    
    # MASKING
    try : mask = self['MASK'].copy()
    except AttributeError: mask = self['MASK']
    active = num.where(mask == 0)
    
    x0 = self['MXMIN_IMAGE']
    y0 = self['MYMIN_IMAGE']
    
    x2, y2, xy = moments(img,x0,y0,order=2,active=active)
    
    p1 = (x2 + y2) / 2.
    p2 = num.sqrt(((x2 - y2)/2.)**2. + xy**2.)
    A = num.sqrt(p1 + p2)
    B = num.sqrt(num.abs(p1 - p2)) # ?? sometimes p1-p2 is negative...
    
    # solve for theta
    # tan theta_0 = 2 * xy / (x2 - y2)
    disc = (xy)/(num.abs(xy))
    try:
        theta_1 = (num.arctan(2*xy/(x2 - y2)))/2. # radians
        # sign of theta must be the same as that of xy
        
        if theta_1 * disc > 0. : theta = theta_1
        elif theta_1 * disc <=0.:
            theta = theta_1 + disc * num.pi/2.
        
    except ZeroDivisionError:
        theta = 0.
    # convert to decimal degrees
    theta = theta * 180. / (num.pi)
    
    ellip = 1. - B/A
    elongation = A / B
    
    # self['M_ELLIP'] = {}
    self['M_A'] = A
    self['M_B'] = B
    self['M_THETA'] = theta
    self['M_ELONG'] = elongation
    self['M_ELLIP'] = ellip
    
    return None
    
    
def moments(img,x0,y0,order,active=-1):
    """Returns moment of given order.
    
    Coordinates refered to stamp!!"""
    
    # IMPORT STUFF
    # END IMPORT
    
    # Allows 'external' masking while invoking
    if active == -1 : 
        active = num.where(img == img)
     
    xabs = active[1].copy() + x0
    yabs = active[0].copy() + y0
    
    sumx = num.sum(img[active] * xabs)
    sumy = num.sum(img[active]* yabs)
    deno = num.sum(img[active])
    x = sumx / deno 
    y = sumy / deno 
    
    if order == 2:
        sumx2 = num.sum(img[active] * xabs**2.)
        sumy2 = num.sum(img[active]* yabs**2.)
        sumxy = num.sum(img[active]* xabs * yabs)
        x2 = sumx2 / deno - x**2
        y2 = sumy2 / deno - y**2
        xy = sumxy / deno - x * y
    
    if order == 1: moments = (x,y)
    if order == 2: moments = (x2,y2,xy)
    
    return moments
    
def snr(self):
    """Updates in self the average Signal to Noise ratio per pixel."""
    # IMPORT STUFF
    # END IMPORT 
    
    sigma_sky = self.execpars['sigma_sky'][0]
    image = self['STAMP'].copy()
    
    Dosnr = True
    
    # ALTERNATIVE SHOT
    if int(sigma_sky) is -1:
        sky_img = image.copy()
        sky_mask = self['SKYMASK']
        sky_img[num.where(sky_mask != 0)] = -100.
        sigma_sky = get_stat(sky_img,'stddev',minimum=-99)
        if sigma_sky == None:
            snr = 0.
            self['SNR'] = snr
            Dosnr = False
        del sky_img, sky_mask
    
    if Dosnr:
    
        sky = self['BACKGROUND']
        try : mask = self['MASK'].copy()
        except AttributeError : mask = self['MASK']
        
        if mask is -1 : mask = num.zeros(shape=image.getshape(),\
        type='Int8')    
        active = num.where(mask == 0)
        x = image[active] - sky
        
        snr = x / sigma_sky # num.sqrt(sigma_sky**2.+x)
        
        try: snr = snr.mean()
        except ZeroDivisionError : snr = 0.
        print 'snr = %5.3f\n' % snr
        self['SNR'] = snr

def increase_window(image,bias=0):
    """Increases the size of a STAMP for the correct operation of 
    asymmetry algorithms."""
    # IMPORT STUFF
    # END IMPORT
    oldshape = image.shape
    valued = num.where(image != 0)
    boxx = max(valued[1]) - min(valued[1]) + 1 + 10
    boxy = max(valued[0]) - min(valued[0]) + 1 + 10
    buffer_lx = min(valued[1])
    buffer_rx = oldshape[1] - (max(valued[1])+1)
    buffer_by = min(valued[0])
    buffer_ty = oldshape[0] - (max(valued[0])+1)
    
    incbuffer_lx = 0 ; incbuffer_rx = 0
    incbuffer_by = 0 ; incbuffer_ty = 0
    if boxx > buffer_lx : incbuffer_lx = boxx - buffer_lx
    if boxx > buffer_rx : incbuffer_rx = boxx - buffer_rx
    if boxy > buffer_by : incbuffer_by = boxy - buffer_by
    if boxy > buffer_ty : incbuffer_ty = boxy - buffer_ty
    
    increased = num.zeros(shape=(oldshape[0]+incbuffer_by+incbuffer_ty,\
    oldshape[1]+incbuffer_lx+incbuffer_rx),type='Float32') + bias
    increased[incbuffer_by:incbuffer_by+oldshape[0],\
    incbuffer_lx:incbuffer_lx+oldshape[1]] = image.copy()
    return increased,incbuffer_lx,incbuffer_by
    

def peak(self,dograph=True):
   """Updates in self the peak center: xpeak, ypeak."""
   # IMPORT STUFF
   from numpy.convolve import boxcar
   # END IMPORT
   
   # INPUTS
   img = self['STAMP'].copy()
   try : mask = self['MASK'].copy()
   except AttributeError : mask = self['MASK']
   box2petro = self.execpars['box2petro'][0]
   x0 = self['MXMIN_IMAGE']
   y0 = self['MYMIN_IMAGE']
   if self.execpars['useExPetro'][0] == 1:
        r_petro = self['EX_R_PETRO']
   else: r_petro = self['R_PETRO']
   sky = self['BACKGROUND']
   # END INPUTS
   
   # the image should be smoothed first
   # In accordance to Lotz et al. 2004, a boxcar filter of side = 0.25 * rp is used.
   
   img = img - sky
   
   boxwidth = int(num.around(box2petro * r_petro))
   if boxwidth % 2 == 0 : boxwidth += 1 # always odd scale.
   if boxwidth < 3 : boxwidth = 3  # Minimum scale for smoothing!
   
   def maskit(stamp,mask):
        masked = num.zeros(shape=stamp.shape,type=stamp.type())
        active = num.where(mask == 0)
        masked[active] = stamp[active]
        return masked
   
   img = maskit(img,mask)
   simg = boxcar(img,(boxwidth,boxwidth))
   peak = num.where(simg == simg.max())
   if len(peak[0]) > 1 :
        dummie = num.ones(shape = simg.getshape(),type = 'Float32')
        xpeak, ypeak = moments(dummie,x0,y0,order=1,active=peak)
   else: xpeak, ypeak = peak[1][0] + x0, peak[0][0] + y0
   
   self['M_XPEAK'] = xpeak
   self['M_YPEAK'] = ypeak
   
def getmoments(self):
    """Returns moments of order 1 and 2."""
    # IMPORT STUFF
    from algorithms import moments
    # END IMPORT
    
    img = self['STAMP'].copy()
    sky = self['BACKGROUND']
    # MASKING.
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    try: img -= sky
    except : stop()
    active = num.where(mask == 0)
    
    x0 = self['MXMIN_IMAGE']
    y0 = self['MYMIN_IMAGE']
    
    x,y = moments(img,x0,y0,order=1,active=active)
    x2, y2, xy = moments(img,x0,y0,order=2,active=active)
    
    print 'X2 = %7.2f\n' % x2
    print 'Y2 = %7.2f\n' % y2
    
    self['M_X'] = x
    self['M_Y'] = y
    self['M_X2'] = x2
    self['M_Y2'] = y2
    self['M_XY'] = xy
    
def makeskystamp(self,sketchimage,bads=0):
    """Makes a fake stamp of the sky."""
    # IMPORT STUFF
    from random import shuffle
    # from algorithms_II import get_stat
    # END IMPORT
    
    # INPUTS
    sky_image = self['STAMP'].copy()
    skymask = self['SKYMASK'].copy()
    # END INPUTS
    
    sky_sigma= self.execpars['sigma_sky'][0]
    sky_image[num.where(skymask != 0)] = -100.
    
    if sky_sigma is -1:
        sky_sigma = get_stat(sky_image,'stddev',minimum=-99)
        if sky_sigma == None : sky_sigma = 1000000
    
    sky_median = get_stat(sky_image,'midpt',minimum=-99)
    
    if sky_median == None : sky_median = 0.
    
    notlowvals = sky_image > sky_median - 3. * sky_sigma
    nothighvals = sky_image < sky_median + 3. * sky_sigma
    
    sky_image[num.where(notlowvals ^ nothighvals)] = -1000.
    
    sky_vals = sky_image[num.where(sky_image != -1000.)]
    skytemp = sketchimage.copy()
    tofill = num.where(skytemp != 0)
    skytemp[tofill] = 1000.
    
    len_sky_vals = len(sky_vals)
    len_skytemp = len(tofill[0])
    
    if len_skytemp <= len_sky_vals :
        fake_sky = sky_vals.copy()
        shuffle(fake_sky)
        skytemp[tofill] = fake_sky[0:len_skytemp].copy()
    else :
        fake_sky = num.zeros((len_skytemp/len_sky_vals+1)*len_sky_vals)
        for i in range(len_skytemp/len_sky_vals+1) :
            fake_sky[i*len_sky_vals:(i+1)*len_sky_vals] = sky_vals.copy()
        shuffle(fake_sky)
        skytemp[tofill] = fake_sky[0:len_skytemp].copy()
    
    self['SKYSTAMP'] = skytemp.copy()
    return skytemp.copy()

def get_stat(image,param,minimum=-99.):
    """Gets standard deviation of an image using iraf."""
    # IMPORT STUFF
    from os import access, F_OK, system
    from pyraf import iraf
    import pyraf
    import pyfits
    from time import time
    # END IMPORT
    
    #       image - the image name
    #       npix - the number of pixels used to do the statistics
    #       mean - the mean of the pixel distribution
    #       midpt - estimate of the median of the pixel distribution
    #       mode - the mode of the pixel distribution
    #       stddev - the standard deviation of the pixel distribution
    #       skew - the skew of the pixel distribution
    #       kurtosis - the kurtosis of the pixel distribution
    #       min - the minimum pixel value
    #       max - the maximum pixel value
    
    tmpstatfits = 'tmpstat_%f.fits' % time()
    tmpstattxt = 'tmpstat_%f.txt' % time()
    isthere_fits = access(tmpstatfits,F_OK)
    isthere_txt = access(tmpstattxt,F_OK)
    if isthere_fits: system('rm %s' % tmpstatfits)
    if isthere_txt: system('rm %s' % tmpstattxt)
    
    pyfits.writeto(tmpstatfits,image.astype('Float32'))
    
    # if minimum < 0. : corr = 0.9
    # if minimum >= 0 : corr = 1.1
    corr = 1.
    iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr()
    value = 0.
    
    try:
        iraf.imstatistics(images=tmpstatfits,fields=param,\
        lower=corr*minimum,upper='INDEF',nclip=3,lsigma=3,\
        usigma=3,binwidt=0.1,format='no',\
        cache='no',mode='al',Stdout=tmpstattxt)
    except OSError : stop()
    except pyraf.irafglobals.IrafError : 
        value = None
    iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr()
    
    if value != None:
        f = open(tmpstattxt,'r')
        data = f.readlines()
        f.close()
        try: value = float(data[0][0:-1])
        except ValueError : value = None
    
    system('rm %s %s' % (tmpstatfits,tmpstattxt))
    
    return value

