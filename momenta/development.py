#! /usr/bin/env python

from pdb import set_trace as stop
import numpy as num
from flags import addflag, isflagon, allflags


def Basics(self,dograph=False):
    """Returns number of pixels of object, total flux, and therefore
    average intensity."""
    # INPUTS
    sky = self['BACKGROUND']
    img = self['STAMP'].copy() - sky
    # MASKING.
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')
    active = num.where(mask == 0)
    # END INPUTS

    nactive = len(num.where(active)[0])
    totflux = num.sum(img[active])
    avint = totflux / nactive
    self['M_NPIX'] = nactive ; self['M_FLUX'] = totflux
    self['M_AVINT'] = avint
    
def M2(self,dograph=False):
    """Returns 2nd order moment of an object, in pixels^2, relative to first
    order moment center."""
    
    # INPUTS
    xcenter = self['M_X'] - self['MXMIN_IMAGE']
    ycenter = self['M_Y'] - self['MYMIN_IMAGE']
    sky = self['BACKGROUND']
    img = self['STAMP'].copy() - sky
    # END INPUTS
    
    # MASKING.
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')
    active = num.where(mask == 0)
    i = img[active]
    x = active[1] ; y = active[0]
    r = num.sqrt((x-xcenter)**2.+(y-ycenter)**2.)
    M2 = num.sum(i * r**2.) / num.sum(i)
    self['M_M2'] = M2
    
    print 'M2**0.5 = %.2f' % (num.abs(M2))**0.5
    
def Radii(self,dograph=False):
    """Returns radii which enclose 20%, 50% and 80% of 
    total flux, in pixels. Reuse code in 'concent'."""
    # IMPORT STUFF
    import scipy.interpolate as interp
    # from pylab import plot,show
    # END IMPORT
    
    # INPUTS
    rprof = self['M_RADIAL']
    # END INPUTS
    
    if isflagon(self['flags'],allflags['NORADIAL']):         
	self['M_R20'] = -99.0 ; self['M_R50'] = -99.0 ; self['M_R80'] = -99.0	
	return None
    
    cumulflx = rprof['cumulflx']
    radii = rprof['radii']
    normal = cumulflx/cumulflx[-1] * 100.
    
    def lazy(normal,radii,proc,addflag):
        to_solve_f = normal - proc
        to_solve_s = interp.splrep(radii,to_solve_f,s=0,k=1)
        try: roots = interp.sproot(to_solve_s)
        except : roots = num.array([])
        radius = -1
        flag = 0L
        if len(roots) == 0:  flag = addflag(flag,1L) # no solution
        else : 
            radius = max(roots)
            if len(roots) > 1: flag = addflag(flag,2L) # more than one radius.
        
        return radius, flag
    
    r_20, flag_20 = lazy(normal,radii,20,addflag)
    r_50, flag_50 = lazy(normal,radii,50,addflag)
    r_80, flag_80 = lazy(normal,radii,80,addflag)
    
    if flag_20 != 0L: r_20 = -99.
    if flag_50 != 0L: r_50 = -99.
    if flag_80 != 0L: r_80 = -99.
    
    self['M_R20'] = r_20 ; self['M_R50'] = r_50 ; self['M_R80'] = r_80
    
    print 'r20, r50, r80 = %.2f, %.2f, %.2f' % (r_20,r_50,r_80)


def Radii2(self,dograph=False):
    """Returns radii which enclose 20%, 50% and 80% of 
    total flux, in pixels. Reuse code in 'concent'."""
    # IMPORT STUFF
    import scipy.interpolate as interp
    # from pylab import plot,show
    from momlib import findzeros
    # END IMPORT
    
    # INPUTS
    rprof = self['M_RADIAL']
    # END INPUTS
    
    if isflagon(self['flags'],allflags['NORADIAL']): 
        self['M_R20'] = -99.0 ; self['M_R50'] = -99.0 ; self['M_R80'] = -99.0	
	return None
    
    cumulflx = rprof['cumulflx']
    radii = rprof['radii']
    normal = cumulflx/cumulflx[-1] * 100.
    
    def lazy(normal,radii,proc,addflag):
        to_solve_f = normal - proc
	roots = findzeros(radii,to_solve_f)
        radius = -1
        flag = 0L
        if len(roots) == 0:  flag = addflag(flag,1L) # no solution
        else : 
            radius = max(roots)
            if len(roots) > 1: flag = addflag(flag,2L) # more than one radius.
        
        return radius, flag
    
    r_20, flag_20 = lazy(normal,radii,20,addflag)
    r_50, flag_50 = lazy(normal,radii,50,addflag)
    r_80, flag_80 = lazy(normal,radii,80,addflag)
    
    if flag_20 != 0L: r_20 = -99.
    if flag_50 != 0L: r_50 = -99.
    if flag_80 != 0L: r_80 = -99.
    
    self['M_R20'] = r_20 ; self['M_R50'] = r_50 ; self['M_R80'] = r_80
    
    print 'r20, r50, r80 = %.2f, %.2f, %.2f' % (r_20,r_50,r_80)


def ApFlx(self,dograph=False):
    """Returns Aperture Flux on circular/eliptical apertures, on 
    given radii. Reuses radial profile."""
    # IMPORT STUFF
    import scipy.interpolate as interp
    # END IMPORT
    # INPUTS
    rprof = self['M_RADIAL']
    apflxradius = self['APFLXRADIUS']
    # END INPUTS
    
    
    if isflagon(self['flags'],allflags['NORADIAL']) or (apflxradius < 0.): 
        self['M_APFLX'] = -99.0
        return None
    
    cumulflx = rprof['cumulflx']
    radii = rprof['radii']
    
    if apflxradius < radii.min() or apflxradius > radii.max() :
        apflx =  -99.0
    else :
        try:
            spline = interp.splrep(radii,cumulflx,s=0,k=1)
            apflx = interp.splev(apflxradius,spline)
        except : apflx = -99.0
    
    print '\nAPFLX = %.3e\n' % apflx 
    self['M_APFLX'] = apflx

    
def AC(self,dograph=False):
    """Returns 'Angular Contrast', a measure of how evenly distributed in 
    P.A. is the flux of an object relative to a certain center. Definition: Given a 
    certain center, the flux of the object in XX deg radial sectors is integrated.  
    'AC' is the maximum difference in contained flux amongst oposed sectors, 
    in absolute value, divided by the flux in all sectors. The center is given by 
    the first order moment of the intensities. For an even angular distribution, 
    AC = 0., and for an object with all flux locked in a single sector, AC = 1."""
    # IMPORT STUFF
    # END IMPORT 
    
    # INPUTS
    sky = self['BACKGROUND']
    image = self['STAMP'].copy() - sky
    xcenter = self['M_X'] - self['MXMIN_IMAGE']
    ycenter = self['M_Y'] - self['MYMIN_IMAGE']
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')
    sector8 = 45.
    sector4 = 90.
    # END INPUTS
    
    image[num.where(mask != 0)] = 0.
    image[int(ycenter),int(xcenter)] = 0.
    
    # Create array with PAs
    
    shape = image.shape
    trshape = (shape[1],shape[0])
    x = num.arange(shape[1]) * \
    num.ones(shape,type='Float32')
    y = num.transpose(num.arange(shape[0])*\
    num.ones(trshape,type='Float32'))
    
    x -= xcenter
    y -= ycenter
    x[num.where(x==0)] = 99999.
    tan = y / x
    angles = num.arctan(tan) * 180. / num.pi
    xneg = num.where(x<0.)
    xposyneg = num.where((x>0.) & (y<0.))
    angles[xneg] += 180.
    angles[xposyneg] += 360.
    angles[num.where( (x  == 99999.) & (y>=0.))] = 90.
    angles[num.where( (x  == 99999.) & (y<0))] = 270.
    
    # Add and substract intensities in sectors
    
    totflux = num.sum(num.sum(image))
    
    def getAC(angles,image,sector):
        
        right = 0.
        left = right + sector
        insector = []
        while right < 180.:
            indexes1 = num.where((angles>right) & \
            (angles <=left))
            indexes2 = num.where((angles>(right+180.)) & \
            (angles <=(left+180.)))
            acc1 = num.sum(image[indexes1])
            acc2 = num.sum(image[indexes2])
            insector.append(float(num.abs(acc1-acc2)))
            
            right += sector
            left = right + sector
    
        insector = num.array(insector)
        insector /= totflux
        AC = insector.max() * 100.
        rightmax = sector * insector.argmax()
        ACMEAN = insector.mean() * 100.
        
        return AC,rightmax,ACMEAN
    
    AC8,right8,AC8M = getAC(angles,image,sector8)
    AC4,right4,AC4M = getAC(angles,image,sector4)
    
    self['M_AC8'] = AC8
    self['M_AC8M'] = AC8M
    self['M_AC4'] = AC4
    self['M_AC4M'] = AC4M
    
    print 'AC8 = %.2f o/o, AC4 = %.2f o/o' %\
    (AC8,AC4)
    
    if dograph:
        
        masksector = (1-mask) & ((angles>right8) & \
        (angles<=(right8+sector8)) | ((angles>(right8+180.)) &\
        (angles<=(right8+sector8+180.))))
        
        self.AC_graph(image,masksector)
    
def FindPeaksII(self,dograph=False):
    """Returns number of 'significative peaks' in an image."""
    # IMPORT STUFF
    import numpy as num
    from pdb import set_trace as stop
    from numpy.nd_image import shift
    from flags import allflags,addflag
    # END IMPORT
    
    # INPUTS
    sky = self['BACKGROUND']
    image = self['STAMP'].copy() - sky
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    # ibin = 0.0001 # !!! f(zeropoint, magbin, pixel,sigma_sky?)
    # 'quantum' of intensity
    radius = max(3,2**0.5) # Pixel radius to discard unreliable
    # secondary peaks. ~ FHWM
    nbins = 10.
    sigma_sky = self.execpars['sigma_sky'][0]
    nsigma = 3.
    # END INPUTS
    
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')
    
    image[num.where(mask != 0)] = 0.
    
    DoPeaks = True
    
    if sigma_sky is -1:
        sky_image = self['STAMP'].copy()
        skymask = self['SKYMASK'].copy()
        sky_image[num.where(skymask != 0)] = -100.
        sigma_sky = get_stat(sky_image,'stddev',minimum=-99)
        if sigma_sky == None : 
            DoPeaks = False
            self['flags'] = addflag(self['flags'],allflags['NOPEAKS'])
        del sky_image,skymask,sky_sigma
    
    if DoPeaks:
     
        maxval = image.max()
        ibin = maxval / nbins
        # More than 3 levels of intensity
        if (maxval-sigma_sky*nsigma)/ibin > 3 :
            nlevels = int(num.ceil(maxval/ibin))
            ilevels = num.arange(nbins) * ibin
            detect = image * 0.
            # Quantization of values
            for ilevel in ilevels: detect[num.where((image>=ilevel) &  \
            (image<ilevel+ibin))] = ilevel
            
            # Find Local Maxima in "quantizated" "detect" image.
            maxima = num.ones(shape=image.shape,type='Bool')
            for i in range(-1,2,1):
                for j in range(-1,2,1):
                    if i!=0 or j!=0:
                        tmpshift = detect.copy() * 0. 
                        shift(detect,(i,j),output = tmpshift)
                        maxima = maxima & (detect >= tmpshift) & \
                        (detect >= nsigma * sigma_sky)
            
            # Solve for adjacent maxima
            maxcoo = num.where(maxima)
            maxinten = image[maxcoo]
            intenorder = maxinten.argsort()
            maxinten = maxinten[intenorder]
            maxcoo = (maxcoo[0][intenorder],maxcoo[1][intenorder])
            nmaxini = len(intenorder)
            
            distances = num.zeros(shape=(nmaxini,nmaxini),type='Float32')
            
            x = maxcoo[1]
            y = maxcoo[0]
            for i in range(nmaxini-1):
                x0 = x[i] ; y0 = y[i]
                distances[i,i+1:] = ((x0-x[i+1:])**2.+(y0-y[i+1:])**2.)**0.5
            
            # Trim Tree-of-Maxima (up-down) 
            
            survivors = [nmaxini-1]
            abduced = []
            for i in range(nmaxini-2,-1,-1):
                if i not in abduced:
                    close = num.where((distances[i,:]>0) & (distances[i,:]<radius))
                    if len(close[0]>0) : abduced += close[0].tolist()
                    survivors.append(i)
                else : pass
            
            # Count Maxima
            NPEAKS = len(survivors)
            abduced = (num.array(abduced),)
            abducedcoo= (maxcoo[0][abduced],maxcoo[1][abduced])
            maxima[abducedcoo] = 0
        else:
            nlevels = 0
            NPEAKS = None
            detect  = None
            maxima = None
            
        self['FP_LOC'] = maxima
        self['FP_DETECT'] = detect
        self['M_NPEAKS'] = NPEAKS
        
        if dograph and maxima != None:
            self.FindPeaks_graph()
        
    else:
        self['FP_LOC'] = None
        self['FP_DETECT'] = None
        self['M_NPEAKS'] = -99
    
    return None

def FindPeaks(self,norm=-1,dograph=False):
    """Returns number of 'significative peaks' in an image."""
    # IMPORT STUFF
    import numpy as num
    from pdb import set_trace as stop
    from numpy.nd_image import shift
    from numpy.nd_image.filters import uniform_filter
    import pyfits
    import os
    #from Moments.algorithms import get_stat
    # END IMPORT
    
    # INPUTS
    sky = self['BACKGROUND']
    image = self['STAMP'].copy() - sky
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    sigma_sky = self.execpars['sigma_sky'][0]
    # skystamp = self['SKYSTAMP'].copy()
    # END INPUTS
    
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')    
    
    filtered3 = num.zeros(shape=image.shape,type='Float32')
    # filtered5 = num.zeros(shape=image.shape,type='Float32')
    # filtered7 = num.zeros(shape=image.shape,type='Float32')
    filtered9 = num.zeros(shape=image.shape,type='Float32')
    filtered15 = num.zeros(shape=image.shape,type='Float32')
    
    uniform_filter(image,(3,3),output=filtered3,\
    mode='constant',cval=0.)
    #uniform_filter(image,(5,5),output=filtered5,\
    #mode='constant',cval=0.)
    #uniform_filter(image,(7,7),output=filtered7,\
    #mode='constant',cval=0.)
    uniform_filter(image,(9,9),output=filtered9,\
    mode='constant',cval=0.)
    uniform_filter(image,(15,15),output=filtered15,\
    mode='constant',cval=0.)
    
    # detect = 49. * filtered7 - (49.*(225.*filtered15-49.*filtered7)/176.)
    
    deno = filtered9 - filtered15
    nonzero = num.where(deno > sigma_sky)
    numerator = filtered3 - filtered15
    detect = filtered3 * 0.
    detect[nonzero] = numerator[nonzero] / deno[nonzero]
    detect[num.where(mask !=0)] = 0.
    
    maxima = num.ones(shape=image.shape,type='Bool')
    
    for i in range(-2,3,1):
        for j in range(-2,3,1):
            if i!=0 or j!=0:
                tmpshift = image.copy() * 0. 
                shift(detect,(i,j),output = tmpshift)
                maxima = maxima & (detect > tmpshift) & (tmpshift != 0.)
    
    #relevance = detect > sigma_sky * 7.9 * 10.
    mindetect = 3.
    relevance = detect > mindetect
    relevance2 = image > 3. * sigma_sky
    Fmaxima = maxima & relevance & relevance2
    
    self['FP_LOC'] = Fmaxima
    self['FP_DETECT'] = detect
    self['M_NPEAKS'] = len(num.where(Fmaxima)[0])
    
    if dograph and maxima != None:
        self.FindPeaks_graph()
    
    return None

def FFactor(self,dograph=False):
    """Returns FFactor, some sort of Filling Factor... on tests."""
    # IMPORT STUFF
    import numpy as num
    from flags import addflag, allflags
    # END IMPORT
    toexec = self.execpars['toexec']
    if 'BASICS' in toexec and 'M20' in toexec:
        Icut = self['M_I20']
        Iavg = self['M_AVINT']
        if Iavg != -99. and Icut != -99.:
            FFactor = num.log10(Icut / Iavg)
        else:
            FFactor = -99.
            self['flags'] = addflag(self['flags'],allflags['NOFFACTOR'])
    else:
        # INPUTS
        sky = self['BACKGROUND']
        img = self['STAMP'].copy() - sky
        try: mask = self['MASK'].copy()
        except AttributeError : mask = self['MASK'].copy()
        # END INPUTS
        if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')
        
        active = num.where(mask == 0)
        
        ord_int = num.sort(img[active])
        
        ord_cumul = num.zeros(shape=ord_int.shape,type='Float32')
        ord_cumul[-1] = ord_int[-1]
        for i in range(len(ord_int)-2,-1,-1):
            ord_cumul[i] = ord_int[i] + ord_cumul[i+1]
        
        ord_cumul = ord_cumul / ord_cumul[0]
        ord_cumul = ord_cumul - 0.2 # brightest pixels which sum up 20% of flux
        cut = num.where(ord_cumul > 0)[0].max()
        Icut = ord_int[cut]
        
        Iavg = ord_int.mean()
        
        if Icut / Iavg > 0.:
            FFactor = num.log10(Icut / Iavg)
        else :
            FFactor = -99.
            self['flags'] = addflag(self['flags'],allflags['NOFFACTOR'])
    #if dograph:
    #    self.FFactor_graph()
    
    self['M_FF'] = FFactor
    print 'Filling Factor = %.2f' % FFactor
    
    return None
