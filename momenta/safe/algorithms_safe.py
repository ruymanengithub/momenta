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
# END IMPORT

def gini_core(x):
    """Computes gini parameter."""
    n = len(x)
    x = num.sort(x)
    # G = (<X> n(n-1) )^-1 * SUM(i->n) (2i - n - 1) Xi
    mean_x = x.mean()
    summatory = sum(( 2.* num.arange(1.,n+1.) - n - 1.) * x )
    G = summatory / (mean_x * n * (n-1.))
    return G
    
def gini(self):
    """Returns the gini parameter of a distribution of pixels.
    Algorithm should be checked.
    """
    # IMPORT STUFF
    import sys
    from algorithms import gini_core
    # END IMPORT
    # INPUTS
    sky = self['BACKGROUND']
    image = self['STAMP'].copy() - sky
    try: mask = self['MASK'].copy()
    except AttributeError: mask = self['MASK']
    # END INPUTS
    
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')    
    active = num.where((mask == 0) & (image > 0.))
    x = image[active]
    # if x.min() < 0. : x -= x.min() # !!!
    # x = x[num.where(x >=0)]
    x = num.abs(x)
    
    G = gini_core(x)
    
    if G > 1. or G < 0. : 
        print 'G = %8.4f, stopping...' % G
        sys.exit()
    
    print 'G = %4.2f' % G
    self['M_GINI'] = G
    
def M20(self):
    """Returns M20 as defined by Lotz et al., 2004... or quasi!"""
    # IMPORT STUFF
    import sys
    from flags import addflag,allflags
    # END IMPORT
    # INPUTS
    sky = self['BACKGROUND']
    img = self['STAMP'].copy() - sky
    try: mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK'].copy()
    x0 = self['MXMIN_IMAGE']
    y0 = self['MYMIN_IMAGE']
    # END INPUTS

    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')

    active = num.where((mask == 0) & (img > 0))
    #x2, y2, xy = moments(img,x0,y0,order=2,active=active)
    
    ord_int = num.sort(img[active])
    ord_cumul = num.zeros(shape=ord_int.shape,type='Float32')
    ord_cumul[-1] = ord_int[-1]
    for i in range(len(ord_int)-2,-1,-1):
        ord_cumul[i] = ord_int[i] + ord_cumul[i+1]
    
    ord_cumul = ord_cumul / ord_cumul[0]
    ord_cumul_20 = ord_cumul - 0.2
    ord_cumul_90 = ord_cumul - 0.9
    M20 = -99
    try: 
        def get_active(img,ord_cumul,ord_int,mask):
            #trimmer = num.where(ord_cumul == ord_cumul[num.where(ord_cumul>0)].min())[0][0]
            trimmer = num.where(ord_cumul > 0)[0].max()
            cut = ord_int[trimmer]
            active = num.where((mask == 0) & (img>cut))
            
            return active,cut
            
        active20,cut20 = get_active(img,ord_cumul_20,ord_int,mask)
        active90,cut90 = get_active(img,ord_cumul_90,ord_int,mask)    
        print 'Cut20 at %f, Cut90 at %f, background included' % (cut20+sky,cut90+sky)
        
        x2_20, y2_20, xy_20 = moments(img,x0,y0,order=2,active=active20)
        x2_90, y2_90, xy_90 = moments(img,x0,y0,order=2,active=active90)

        #M2 = num.sqrt(num.abs(x2)+num.abs(y2))
        M2_20 = num.sqrt(num.abs(x2_20)+num.abs(y2_20))
        M2_90 = num.sqrt(num.abs(x2_90)+num.abs(y2_90))
        
        M20 = num.log10(M2_20/M2_90)
        
    except ZeroDivisionError:
        self['flags'] = addflag(self['flags'],allflags['NOM20'])
    
    print 'M20 = %4.2f' % M20
    self['M20'] = M20
    
    
def ellipse(self):
    """Returns A, B, THETA elongation and ellipticity of distribution of pixels."""
    # IMPORT STUFF
    from algorithms import moments
    # END IMPORT
    
    sky = self['BACKGROUND']
    img = self['STAMP'].copy() - sky
    
    # MASKING.
    try : mask = self['MASK'].copy()
    except AttributeError: mask = self['MASK']
    active = num.where(mask == 0)
    
    x0 = self['MXMIN_IMAGE']
    y0 = self['MYMIN_IMAGE']
    
    x2, y2, xy = moments(img,x0,y0,order=2,active=active)
    
    p1 = (x2 + y2) / 2.
    p2 = num.sqrt(((x2 - y2)/2.)**2. + xy**2.)
    A = num.sqrt(p1 + p2)
    B = num.sqrt(p1 - p2)
    
    # solve for theta
    # tan theta_0 = 2 * xy / (x2 - y2)
    
    try:
        theta = (num.arctan(2*xy/(x2 - y2)))/2. # radians
        # sign of theta must be the same as that of xy
        theta = theta * (xy)/(num.abs(xy))
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
    # ALTERNATIVE SHOT
    if int(sigma_sky) is -1:
        sky_img = image.copy()
        sky_mask = self['SKYMASK']
        sky_img[num.where(sky_mask != 0)] = -100.
        sigma_sky = get_stat(sky_img,'stddev',minimum=-99)
        del sky_img, sky_mask
    
    sky = self['BACKGROUND']
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')    
    active = num.where(mask == 0)
    x = image[active] - sky
    x = x[num.where(x>0.)] # NO NEGATIVES
    
    snr = x / num.sqrt(sigma_sky**2.+x)
    
    try: snr = snr.mean()
    except ZeroDivisionError : snr = 0.
    print 'snr = %5.3f\n' % snr
    self['SNR'] = snr
    
def axis_asymmetry(self,axis):
    """Updates in self the degree of axis symmetry respect to minor/major semi-major axis.
    
    Ax = sum | I(i,j) - I_180(i,j)| / (2 * sum| I(i,j) |)  - Ax_180 
    
    where I is the galaxy's image and I_180 is the image bended about the major or minor axis of the galaxy, and
    B_180 is the average asymmetry of the background. Ax is summed over all pixels within 1.5 rp of the galaxy's
    center. The central pixel is fixed. Which one to use??"""
    # IMPORT STUFF
    from flags import addflag
    import pyfits
    # END IMPORT
    
    # INPUTS
    image = self['STAMP'].copy()
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    center = (self['Y_IMAGE']-self['MYMIN_IMAGE'],self['X_IMAGE']-self['MXMIN_IMAGE'])
    pa = self['THETA_IMAGE']  # RELATIVE TO X AXIS
    sky = self['BACKGROUND']
    doAxisAsym_sky = self.execpars['doAxisAsym_sky'][0] == 1
    doAxisAsym_iter = self.execpars['doAxisAsym_iter'][0] == 1
    mode = 'iraf' # iraf, numpy
    # END INPUTS
    
    if axis == 'minor' or axis=='Minor': 
        pa = 90. + pa
        label = '_MIN'
    else: label = '_MAJ'
    
    y, x = center
    
    masked = num.zeros(shape = image.shape,type=image.type())
    active = num.where(mask == 0)
    masked[active] = image[active] - sky
    del active
    
    masked,xshift,yshift = increase_window(masked)
    x += xshift ; y += yshift
    center = (y,x)
    
    if doAxisAsym_sky:
        skystamp = self.makeskystamp(masked,bads=0)
        active = num.where(skystamp != 0)
        skystamp[active] = skystamp[active] - sky
    del active
    
    axsky = 0.0
    bunch, Ax_flags = axisasym_solver(masked,center,pa,doAxisAsym_iter,axsky,norm=-1,mode=mode)
    Ax_center = bunch[0]
    Aximg = bunch[1]
    
    self['M_AXS'+label] = Aximg
    self['flags'] = addflag(self['flags'],Ax_flags)
    
    if doAxisAsym_sky:
        ignore = False
        norm = num.abs(masked)
        norm = num.sum(num.sum(norm))
        norm = norm * 2.
        bunchsky,ignore = axisasym_solver(skystamp,Ax_center,pa,ignore,axsky,norm=norm,mode=mode)
        Axsky = bunchsky[1]
        self['M_AXS%s_SKY'%label] = Axsky
    else :
        Axsky = -99.
    
    print 'Ax_iraf%s = %6.4f    Ax_iraf%s_sky = %6.4f' % (label,Aximg,label,Axsky)
    print 'Ax_shift = %4.2f pixels\n' % num.sqrt((center[1]-Ax_center[1])**2.+(center[0]-Ax_center[0])**2.)

def axisasym_solver(image,center,pa,iterative,axsky,norm=-1,mode='iraf'):
    """Returns minimum Axis Asymmetry of given image, and point with respect to which, Axis Asymmetry is 
    such a minimum."""
    # IMPORT STUFF
    from scipy.optimize import fmin
    from flags import addflag, allflags
    # END IMPORT STUFF
    
    flags = 0L
   
    if iterative :
        bunch = fmin(axisasym_func,center,args=(image,pa,axsky,norm,mode),full_output=1,maxiter = 1000.,disp=0)
        # xopt = bunch[0] # miminum A center
        # fopt = bunch[1] # A
        # iter = bunch[2] # number of iterations
        # funcalls = bunch[3] # number of function calls
    else :
        Ax = axisasym_func(center,image,pa,axsky,norm,mode)
        bunch = [[center[0],center[1]],Ax,0,0,0]
    warnings = bunch[4] # warnings: 1= max number of func. evals, 2=max number of iterations
    if warnings == 2 : flags = addflag(flags,allflags['MAXITER_AXAS'])
    
    return bunch, flags
    
def axisasym_func(center,image,pa,zero=0,norm=-1,mode='iraf'):
    """Computes degree of axis symmetry for generic image and axis of symmetry. Uses either iraf or numpy."""
    # IMPORT STUFF
    from pyraf import iraf
    import pyfits
    import os
    from algorithms_II import bend180numpy
    # END IMPORT    
    
    x = center[1] ; y = center[0]
    
    Ax_flags = 0L
    
    if mode == 'iraf':
        image1, image2 = bend180iraf(image,x,y,pa)
    elif mode == 'numpy': 
        image1, image2 = bend180numpy(image,x,y,pa)
    
    # nonzero = num.where(True ^ ((image1 != 0) ^ (image2 != 0)))
    nonzero = num.where(image1 == image1)
    up = num.abs(image1[nonzero] - image2[nonzero])
    up = num.sum(up)
    # from pylab import imshow,show
    # print 'showing Axis symetry image... %i pixels left out' % len(blank[0])
    # imshow(up,origin='lower') ; show()
    
    if norm == -1:
        norm = num.abs(image1)
        norm = num.sum(num.sum(norm))
        norm = norm * 2.
    
    Ax = up / down - zero
    
    return Ax # , Ax_flags
  
def bend180iraf(image,x,y,pa):
    """Bends an image over an axis, using iraf."""
    # IMPORT STUFF
    import os
    import pyfits
    from pyraf import iraf
    from time import time
    # END IMPORT
    
    tmpimage1 = 'tmpbend180_1_%f.fits' % time()
    tmpimage2 = 'tmpbend180_2_%f.fits' % time()
    
    if os.path.exists(tmpimage1) : os.system('rm %s' % tmpimage1)
    if os.path.exists(tmpimage2) : os.system('rm %s' % tmpimage2)
    
    pyfits.writeto(tmpimage1,image)
    
    iraf.load('images',doprint=0)
    iraf.images.imgeom()
    
    xshift = ((image.shape[1]-1)/2.) - x #
    yshift = ((image.shape[0]-1)/2.) - y #
    
    # print 'xshift = %f, yshift = %f' % (xshift,yshift)
    iraf.images.imgeom.imshift(input=tmpimage1,output=tmpimage1,xshift = xshift,yshift = yshift,interp='linear',boundar='constant',constan=0.)
    iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() 
    iraf.images.imgeom.imshift.unlearn()
    
    iraf.images.imgeom.rotate(input=tmpimage1,output=tmpimage1,rotation=-pa,xin='INDEF',yin='INDEF',\
    xout='INDEF',yout='INDEF',interpo='linear',boundar='constant',constan=0,verbose='no')
    
    iraf.images.imgeom.imlintran(input=tmpimage1,output=tmpimage2,xrotatio=0,yrotatio=180,xmag=1,\
    ymag=1,xin='INDEF',yin='INDEF',xout='INDEF',yout='INDEF',ncols=0.,nlines=0.,interpo='linear',\
    boundar='constant',constan=0.,fluxcon='yes',nxblock=512,nyblock=512,verbose='no')
    
    image1 = pyfits.getdata(tmpimage1)
    image2 = pyfits.getdata(tmpimage2)
    
    os.system('rm %s %s' % (tmpimage1,tmpimage2))
    
    return image1, image2
    
def asymmetry(self):
    """Updates in self the 'Asymmetry' (A) parameter, according to the definition given in Lotz. et al., 2004. 
    
    A = sum | I(i,j) - I_180(i,j)| / (2 * sum| I(i,j) |)  - B_180 
    
    where I is the galaxy's image and I_180 is the image rotated by 180deg about the galaxy's central pixel, and
    B_180 is the average asymmetry of the background. A is summed over all pixels within 1.5 rp of the galaxy's
    center. The central pixel is determined by minimizing A.
  
    """
    # IMPORT STUFF
    from flags import addflag,isflagon,allflags
    import pyfits
    import os
    # from pylab import imshow,show
    # END IMPORT 
    
    # INPUTS
    stamp = self['STAMP'].copy()
    xcenter = self['X_IMAGE'] - self['MXMIN_IMAGE']
    ycenter = self['Y_IMAGE'] - self['MYMIN_IMAGE']
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    sky = self['BACKGROUND']
    iterative = self.execpars['Asym_iterative'][0] == 1
    doAsym_sky = self.execpars['doAsym_sky'][0] == 1
    mode = 'iraf' # iraf, numpy
    # END INPUTS
    
    # To speed up...
    if num.any(num.array(stamp.shape) > 500) : mode = 'iraf'
    
    masked = num.zeros(shape=stamp.shape,type=stamp.type()) 
    active = num.where(mask == 0)
    masked[active] = stamp[active] - sky
    masked,xshift,yshift = increase_window(masked)
    xcenter += xshift
    ycenter += yshift
    center = (ycenter,xcenter)
    del active
    
    if doAsym_sky:
        skystamp = self.makeskystamp(masked,bads=0)
        active = num.where(skystamp != 0)
        skystamp[active] = skystamp[active] - sky
        del active
    
    B_180 = 0.
    
    bunch, A_flags = asym_solver(masked,center,B_180,iterative,mode=mode)
    A_center = bunch[0]
    A = bunch[1]
    
    if isflagon(A_flags,allflags['MAXITER_AS']):
        bunch, dump = asym_solver(masked,center,B_180,iterative=False,mode=mode)
        A_center = bunch[0]
        A = bunch[1]
    
    if doAsym_sky:
        norm = num.abs(masked)
        norm = num.sum(num.sum(norm))
        norm = norm * 2.
        bunch_sky,ignore = asym_solver(skystamp,A_center,B_180,iterative=False,norm=norm,mode=mode)
        A_sky = bunch_sky[1]
        self['M_AS_SKY'] = A_sky
    else : A_sky = -99
  
    self['M_AS_X'] = A_center[1] - xshift + self['MXMIN_IMAGE']
    self['M_AS_Y'] = A_center[0] - yshift + self['MYMIN_IMAGE']
    self['M_AS'] = A
    
    print '\nAsym = %4.2f    Asym_sky = %4.2f' % (A,A_sky)
    print 'A_shift = %4.2f pixels\n' % num.sqrt((center[1]-A_center[1])**2.+(center[0]-A_center[0])**2.)
    
    self['flags'] = addflag(self['flags'],A_flags)
    
def asym_solver(image,center,asky,iterative,norm=-1,mode='iraf'):
    """Returns minimum 180deg Asymmetry of given image, and point 
    with respect to which, Asymmetry is such a minimum."""    
    # IMPORT STUFF
    from scipy.optimize import fmin
    from flags import addflag, allflags
    # END IMPORT STUFF
    
    flags = 0L
    maxiter = 50
    
    if iterative :
        bunch = fmin(asym_func,center,args=(image,asky,norm,mode),full_output=1,maxiter = maxiter,disp=0)
        # xopt = bunch[0] # miminum A center
        # fopt = bunch[1] # A
        # iter = bunch[2] # number of iterations
        # funcalls = bunch[3] # number of function calls
    else :
        A = asym_func(center,image,asky,norm,mode)
        bunch = [[center[0],center[1]],A,0,0,0]
    warnings = bunch[4] # warnings: 1= max number of func. evals, 2=max number of iterations
    if warnings == 2 : 
        flags = addflag(flags,allflags['MAXITER_AS'])
        print 'MAXIMUM NUMBER OF ITERATIONS (%i) reached in ASYMMETRY' % maxiter
    
    return bunch, flags
    
def asym_func(center,image,zero,norm=-1,mode='iraf'):    
    """Returns asymmetry value given masked image and center."""
    # IMPORT STUFF
    import os
    from algorithms_II import rot180numpy
    #END IMPORT
    
    if mode == 'iraf':
        image1,image2 = rot180iraf(image,center)
    elif mode == 'numpy':
        image1, image2 = rot180numpy(image,center)
    
    # nonzero = num.where(True ^ ((image1 != 0) ^ (image2 != 0)))
    nonzero = num.where(image1 == image1)
    up = num.abs(image1[nonzero] - image2[nonzero])
    up = num.sum(up)
    if norm == -1:
        down = num.abs(image1)
        down = num.sum(num.sum(down))
        norm = down * 2.
    
    A = up / norm - zero
    
    return A
 
def rot180iraf(image,center):
    """Rotates an image 180 degrees respect to a center, using iraf."""    
    # IMPORT STUFF
    import pyfits
    from pyraf import iraf
    import os
    from time import time
    # END IMPORT
    
    #print time()
    
    tmpimage1 = 'tmp180_1_%f.fits' % time()
    tmpimage2 = 'tmp180_2_%f.fits' % time()
    
    if os.path.exists(tmpimage1) : os.system('rm %s' % tmpimage1)
    if os.path.exists(tmpimage2) : os.system('rm %s' % tmpimage2)
    
    pyfits.writeto(tmpimage1,image)
     
    iraf.load('images',doprint=0)
    iraf.images.imgeom()
    
    x = center[1] ; y = center[0]
    xshift = ((image.shape[1]-1)/2.) - x # 
    yshift = ((image.shape[0]-1)/2.) - y # 
    
    iraf.images.imgeom.imshift(input=tmpimage1,output=tmpimage1,xshift = xshift,yshift = yshift,\
    interp='linear',boundar='constant',constan=0.)
    iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() 
    
    iraf.images.imgeom.rotate(input=tmpimage1,output=tmpimage2,rotation=180,xin='INDEF',yin='INDEF',\
    xout='INDEF',yout='INDEF',interpo='linear',boundar='constant',constan=0,verbose='no')
    
    image1 = pyfits.getdata(tmpimage1)
    image2 = pyfits.getdata(tmpimage2)
    os.system('rm %s %s' % (tmpimage1,tmpimage2))
    
    return image1,image2
 
def increase_window(image,bias=0):
    """Increases the size of a STAMP for the correct operation of asymmetry algorithms."""
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
    increased[incbuffer_by:incbuffer_by+oldshape[0],incbuffer_lx:incbuffer_lx+oldshape[1]] = \
    image.copy()
    return increased,incbuffer_lx,incbuffer_by
    
def clumpy(self):
    """Updates in self the 'Clumpiness' parameter (S) according to the definition given in Lotz et al., 2004.
    
    S = sum | I(i,j) - Is(i,j)| / sum | I(i,j) | - Bs
    
    where Is is the galaxy's image smoothed by a boxcar of width 0.25 rp (rp = petrosian radius), and Bs is 
    the average smoothness of the background. Like A, S is summed over the pixels within 1.5rp of the galaxy's
    center. Bs is average 'clumpiness' of sky.
    """
    # IMPORT STUFF
    from numpy.convolve import boxcar
    from flags import addflag
    import pyfits
    # END IMPORT
    
    # INPUTS
    sky = self['BACKGROUND']
    stamp = self['STAMP'].copy()
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    if self.execpars['useExPetro'][0] == 1:
        r_petro = self['EX_R_PETRO']
    else: r_petro = self['R_PETRO']
    box2petro = self.execpars['box2petro'][0]
    doClumpy_sky = self.execpars['doClumpy_sky'][0] == 1
    # END INPUTS
    
    masked = num.zeros(shape=stamp.shape,type=stamp.type())
    active = num.where(mask == 0)
    masked[active] = stamp[active] - sky
    if doClumpy_sky : 
        skystamp = self.makeskystamp(masked,bads=0)
        skystamp[active] = skystamp[active] - sky
    del active
    
    boxwidth = int(num.around(box2petro * r_petro))
    if boxwidth % 2 == 0 : boxwidth += 1
    if boxwidth < 3 : boxwidth = 3

    smoothed = boxcar(masked,(boxwidth,boxwidth),mode='reflect') 

    def g(difimage,image,num):
        up = num.abs(difimage)
        down = num.abs(image)
        up = num.sum(num.sum(up))
        down = num.sum(num.sum(down))
        return up/down
    
    S = g(masked-smoothed,masked,num)
    
    S_flags = 0L
    
    self['M_S'] = S
    self['flags'] = addflag(self['flags'],S_flags)
    if doClumpy_sky :
        S_sky = g(skystamp,masked,num)
        self['M_S_SKY'] = S_sky
    else : S_sky = -99
    
    print 'S = %4.2f   S_sky = %4.2f' % (S,S_sky)
    
    
def concent(self,verbose=False):
    """Updates in self the 'Concentration', acording to the definition given in Bershady et al., 2000 
    (extracted from Lotz et al., 2004)
    
    C = 5 log (r80/r20)
    where r? is the radius of the circular aperture that contains ? % of the total flux of the object. The total flux
    is that contained within 1.5rp of the galaxy's center, as in Lotz et al., 2004 and Conselice 2003.
    How is the center determined? Lotz et al. use that given by the determination of the 'A' parameter.
    """
    # IMPORT STUFF
    import scipy.interpolate as interp
    from flags import addflag, isflagon, allflags
    # from pylab import plot,show
    # END IMPORT
    
    # INPUTS
    rprof = self['M_RADIAL']
    r_big_proc = self.execpars['r_big_proc'][0]
    r_small_proc = self.execpars['r_small_proc'][0]
    # END INPUTS
    
    if not isflagon(self['flags'],allflags['NORADIAL']):
        cumulflx = rprof['cumulflx']
        radii = rprof['radii']
        normal = cumulflx/cumulflx[-1] * 100.
        print 'Min radius = %f pixels' % radii[0]
        # plot(radii,normal) ; show()
        
        def lazy(normal,radii,proc,addflag):
            # IMPORT STUFF
            # END IMPORT
            to_solve_f = normal - proc
            to_solve_s = interp.splrep(radii,to_solve_f,s=0)
            roots = interp.sproot(to_solve_s)
            radius = -1
            flag = 0L
            if len(roots) == 0:  flag = addflag(flag,1L) # no solution
            else : 
                radius = max(roots)
                if len(roots) > 1: flag = addflag(flag,2L) # more than one radius.
            
            return radius, flag
        
        r_small, flag_small = lazy(normal,radii,r_small_proc,addflag)
        r_big, flag_big = lazy(normal,radii,r_big_proc,addflag)
        
        cflag = 0L
        if isflagon(flag_small,1L) or isflagon(flag_big,1L) : cflag = addflag(cflag,allflags['NOCONC'])
        if isflagon(flag_small,2L) or isflagon(flag_big,2L):
            cflag = addflag(cflag,allflags['NOCONC'])
            cflag = addflag(cflag,allflags['MANYCONC'])
        C = -99.0
        if not isflagon(cflag,allflags['NOCONC']) : \
        C = 5. * num.log10(r_big / r_small)
        
        print 'C = %5.3f\n' % C
        self['M_C'] = C
        if not isflagon(cflag,allflags['NOCONC']) and\
        not isflagon(cflag,allflags['MANYCONC']) : 
            self['M_RBIG'] = r_big
            self['M_RSMALL'] = r_small
        else:
            self['M_RBIG'] = -99.0
            self['M_RSMALL'] = -99.0
        self['flags'] = addflag(self['flags'],cflag)
        
    else : 
        self['M_C'] = -99.0
        self['flags'] = addflag(self['flags'],allflags['NOCONC'])
    
   
    if verbose : return C,r_small,r_big
    
    
def peak(self):
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
   if boxwidth % 2 == 0 : boxwidth += 1
   if boxwidth < 3 : boxwidth = 3  
   
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
    img -= sky
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
    
    return skytemp

def get_stat(image,param,minimum=-99):
    """Gets standard deviation of an image using iraf."""
    # IMPORT STUFF
    from pyraf import iraf
    import pyfits
    import pyraf
    from os import access, F_OK, system
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
    
    pyfits.writeto(tmpstatfits,image)
    
    if minimum < 0. : corr = 0.9
    if minimum >= 0 : corr = 1.1
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
        except ValueError : value = data[0][0:-1]
    
    system('rm %s %s' % (tmpstatfits,tmpstattxt))
    
    return value

def FindPeaks(self):
    """Returns number of 'significative peaks' in an image, and radial 
    distribution."""
    # IMPORT STUFF
    from numpy.nd_image import shift
    from numpy.nd_image.filters import uniform_filter
    from algorithms import gini_core
    # END IMPORT
    
    # INPUTS
    sky = self['BACKGROUND']
    image = self['STAMP'].copy() - sky
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    sigma_sky = self.execpars['sigma_sky'][0]
    xcenter = self['M_X'] ; ycenter = self['M_Y']
    try:
        Rbig = self['M_RBIG'] ; Rsmall = self['M_RSMALL']
    except KeyError :
        Rbig = -99.0 ; Rsmall = -99.0
    # END INPUTS
    
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')    
    image[num.where(mask != 0)] = 0.
    
    filtered3 = num.zeros(shape=image.shape,type='Float32')
    filtered7 = num.zeros(shape=image.shape,type='Float32')
    filtered15 = num.zeros(shape=image.shape,type='Float32')
    
    uniform_filter(image,(3,3),output=filtered3,\
    mode='constant',cval=0)
    uniform_filter(image,(7,7),output=filtered7,\
    mode='constant',cval=0)
    uniform_filter(image,(15,15),output=filtered15,\
    mode='constant',cval=0)
    
    detect = 49. * filtered7 - (49.*(225.*filtered15-49.*filtered7)/176.)
    
    maxima = num.ones(shape=image.shape,type='Bool')
    
    for i in range(-2,3,1):
        for j in range(-2,3,1):
            if i==0 and j==0:
                pass
            else:
                tmpshift1 = detect.copy() * 0.
                shift(detect,(i,j),output = tmpshift1)
                maxima = maxima & (detect > tmpshift1)
    relevance = filtered3 / sigma_sky > 1.
    relevance2 = detect > 3. * 13.13 * sigma_sky #detsigma
    maxima = maxima & relevance & relevance2
    xMAXIMA = num.where(maxima)[1]
    yMAXIMA = num.where(maxima)[0]
    radMAXIMA = num.sqrt((xcenter-xMAXIMA)**2.+(ycenter-yMAXIMA)**2.)
    
    iMAXIMA = detect[num.where(maxima)]
    
    GINI_PEAK = gini_core(iMAXIMA)
    
    if Rbig > 0. and Rsmall > 0.:
        self['M_NPEAKS_S'] = len(num.where(radMAXIMA<= Rsmall)[0])
        self['M_NPEAKS_L'] = len(num.where(radMAXIMA <= Rbig)[0])
    else:
        self['M_NPEAKS_S'] = -99.0
        self['M_NPEAKS_L'] = -99.0
    
    self['M_NPEAKS'] = len(num.where(maxima)[0])
    self['M_GINI_PEAK'] = GINI_PEAK
    
    return None
    
