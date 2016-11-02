"""Set of functions that operate on a image + mask, executing some algorithm.

Functions: gini, ellipse, moments (several order moments)
axis_symmetry, asymmetry (A), clumpiness (S), concentration, peak.
"""

import numpy as num
from pdb import set_trace as stop

def gini(self):
    """Returns the gini parameter of a distribution of pixels.
    Algorithm should be checked.
    """
    # IMPORT STUFF
    import sys
    # END IMPORT
    
    sky = self['BACKGROUND']
    image = self['STAMP'] - sky
    mask = self['MASK']
    
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')    
    active = num.where(mask == 0)
    x = image[active]
    # if x.min() < 0. : x -= x.min() # !!!
    x = x[num.where(x >=0)]
    n = len(x)
    x = num.sort(x)
    
    
    # G = (<X> n(n-1) )^-1 * SUM(i->n) (2i - n - 1) Xi
    mean_x = x.mean()
    summatory = sum(( 2.* num.arange(1.,n+1.) - n - 1.) * x )
    G = summatory / (mean_x * n * (n-1.))
    if G > 1. or G < 0. : 
        print 'G = %8.4f, stopping...' % G
        sys.exit()

    self['M_GINI'] = G
    
    
def ellipse(self):
    """Returns A, B, THETA elongation and ellipticity of distribution of pixels."""
    # IMPORT STUFF
    from algorithms import moments
    # END IMPORT
    
    img = self['STAMP']
    sky = self['BACKGROUND']
    # MASKING.
    mask = self['MASK']
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
    
def axis_symmetry(self,axis):
    """Updates in self the degree of axis symmetry respect to minor/major semi-major axis.
    
    Ax = sum | I(i,j) - I_180(i,j)| / sum| I(i,j) |  - Ax_180 
    
    where I is the galaxy's image and I_180 is the image bended about the major or minor axis of the galaxy, and
    B_180 is the average asymmetry of the background. Ax is summed over all pixels within 1.5 rp of the galaxy's
    center. The central pixel is fixed. Which one to use??"""
    # IMPORT STUFF
    from flags import addflag
    from algorithms import axissym
    # END IMPORT
    
    # INPUTS
    image = self['STAMP']
    mask = self['MASK']
    # mask_sky = self['MASK_SKY']
    center = (self['Y_IMAGE']-self['MYMIN_IMAGE'],self['X_IMAGE']-self['MXMIN_IMAGE'])
    #center_sky = (self['YCENTER_SKY']-self['MYMIN_IMAGE'],self['XCENTER_SKY']-self['MXMIN_IMAGE'])
    pa = self['THETA_IMAGE']  # RELATIVE TO X AXIS
    sky = self['BACKGROUND']
    # END INPUTS
    
    if axis == 'minor' or axis=='Minor': 
        pa = 90. + pa
        label = '_MIN'
    else: label = '_MAJ'
    
    y, x = center
    # xsky,ysky = center_sky
    
    masked = num.zeros(shape = image.shape,type=image.type())
    active = num.where(mask == 0)
    masked[active] = image[active] - sky
    del active
    
    # masked_sky = num.zeros(shape = image.shape,type=image.type()) + sky
    # active = num.where(mask_sky == 0)
    # masked_sky[active] = image[active]
    # del active
    
    # Axsky, ignore = axissym(masked_sky,xsky,ysky,pa,const =0)
    Axsky = 0.0
    Aximg, Ax_flags = axissym_iraf(masked,x,y,pa,Axsky)
    Aximgnoiraf,ignore = axissym(masked,x,y,pa,Axsky)
    print 'Ax_iraf = %4.2f, Ax_noiraf = %4.2f' % (Aximg,Aximgnoiraf)
    
    self['M_AXS'+label] = Aximg
    self['flags'] = addflag(self['flags'],Ax_flags)
    
    
def axissym_iraf(image,x,y,pa,zero=0):
    """Computes degree of axis symmetry for generic image and axis of symmetry. Uses iraf."""
    # IMPORT STUFF
    from pyraf import iraf
    import pyfits
    import os
    # END IMPORT    
    
    Ax_flags = 0L
    
    tmpimage1 = 'tmpbend180_1.fits'
    tmpimage2 = 'tmpbend180_2.fits'
    
    pyfits.writeto(tmpimage1,image)
    
    iraf.load('images',doprint=0)
    iraf.images.imgeom()
    
    xshift = (image.shape[1]/2. + 0.5) - (x+1)
    yshift = (image.shape[0]/2. + 0.5) - (y+1)
    
   
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
    
    blank1 = image1 == 0 ; blank2 = image2 == 0
    blank = num.where( blank1 | blank2)
    image1[blank] = 0. ; image2[blank] = 0.
    
    up = num.abs(image1 - image2)
##    from pylab import imshow,show
##    imshow(up,origin='lower') ; show()
##    imshow(image1,origin='lower') ; show()
    
    up = num.sum(num.sum(up))
    down = num.abs(image1)
    down = num.sum(num.sum(down))
    
    Ax = up / down - zero
    
    return Ax, Ax_flags
    
    
def axissym(image,x,y,pa,zero=0):
    """Computes degree of axis symmetry for generic image and axis of symmetry."""
    # IMPORT STUFF
    # END IMPORT
    
    Ax_flags = 0L
    
    bimage = bend180(x,y,pa,image)
    
    # dummie = num.ones(shape=image.shape,type='Int8')
    # dummie[num.where(image != 0)] = 0
    # dummie[num.where(bimage != 0)] = 0
    # active = num.where(dummie == 0)
    
    # up = num.abs(image[active] - simage[active])
    # up = num.sum(num.sum(up))
    # down = num.abs(image[active])
    up = num.abs(image - bimage)
    up = num.sum(num.sum(up))
    down = num.abs(image)
    
    down = num.sum(num.sum(down))
    
    Ax = up / down - zero
    
    return Ax, Ax_flags
    
def bend180(x,y,pa,image):
    """Returns a rotation of 180deg of given image around given line (defined by x,y and pa). 
    Returned image has same dimensions as input image."""
    # IMPORT STUFF
    # END IMPORT
    
    radeg = 180. / num.pi
    ang = pa / radeg
    
    if abs(pa) != 0 and abs(pa) != 180 and abs(pa) != 90. and abs(pa) != 270:
        m = num.tan(ang) # slope of axis 
        mp = - 1./m
        b = y - m * x
        indexes = num.where(image == image)
        xp = indexes[1].astype('Float32').copy()
        yp = indexes[0].astype('Float32').copy()
        bp = yp - mp * xp
        xcross = (bp - b) / (m - mp)
        ycross = m * xcross + b
        incx = 2.*(xcross - xp)
        incy = 2.*(ycross - yp)
        xp += incx ; yp += incy
        indexesp = (yp,xp)
    elif abs(pa) == 90. or abs(pa) == 270.:
        indexes = num.where(image == image)
        xp = indexes[1].copy()
        xcross = x
        incx = 2.*(xcross - xp)
        xp += incx
        indexesp = (indexes[0],xp)
    elif abs(pa) == 0. or abs(pa) == 180.:
        indexes = num.where(image == image)
        yp = indexes[0].copy()
        ycross = y
        incy = 2.*(ycross - yp)
        yp += incy
        indexesp = (yp,indexes[1])

    def intindex(indexes,num):
        return (num.around(indexes[0]).astype('Int32').copy(),\
        num.around(indexes[1]).astype('Int32').copy())
    
    indexes = intindex(indexes,num)
    indexesp = intindex(indexesp,num)
    
    maxx = max(indexes[1]) -1
    maxy = max(indexes[0]) -1
    
    dummie = num.zeros(shape=(len(indexes[0]),),type='Int8')
    dummie[num.where(indexesp[0] < 0 )] = 1
    dummie[num.where(indexesp[0] > maxy )] = 1
    dummie[num.where(indexesp[1] < 0 )] = 1
    dummie[num.where(indexesp[1] > maxx )] = 1
    
    def trimindex(indexes,mask,num):
        goods = num.where(mask==0)
        return (indexes[0][goods],indexes[1][goods])
    
    indexes = trimindex(indexes,dummie,num)
    indexesp = trimindex(indexesp,dummie,num)
    bimage = num.zeros(shape=image.shape,type='Float32') 
    bimage[indexesp] = image[indexes]
    # from pylab import show, imshow # DEVEL
    # stop()                                                # DEVEL
    
    return bimage

def asymmetry(self):
    """Updates in self the 'Asymmetry' (A) parameter, according to the definition given in Lotz. et al., 2004. 
    
    A = sum | I(i,j) - I_180(i,j)| / sum| I(i,j) |  - B_180 
    
    where I is the galaxy's image and I_180 is the image rotated by 180deg about the galaxy's central pixel, and
    B_180 is the average asymmetry of the background. A is summed over all pixels within 1.5 rp of the galaxy's
    center. The central pixel is determined by minimizing A.
  
    """
    # IMPORT STUFF
    from flags import addflag
    # END IMPORT 
    
    # INPUTS
    stamp = self['STAMP']
    xcenter = self['X_IMAGE'] - self['MXMIN_IMAGE']
    ycenter = self['Y_IMAGE'] - self['MYMIN_IMAGE']
    center = (ycenter,xcenter)
    mask = self['MASK']
    # mask_sky = self['MASK_SKY']
    # ycentersky = self['XCENTER_SKY'] - self['MXMIN_IMAGE']
    # xcentersky = self['YCENTER_SKY'] - self['MYMIN_IMAGE']
    # centersky = (ycentersky,xcentersky)
    sky = self['BACKGROUND']
    iterative = self.execpars['Asym_iterative'][0] == 1
    # END INPUTS
    
    masked = num.zeros(shape=stamp.shape,type=stamp.type()) 
    active = num.where(mask == 0)
    masked[active] = stamp[active] - sky
    del active
    
    # asky = 0
    # B_180 = asym_func(masked_sky,centersky,asky,sky) #  RIGHT??
    B_180 = 0.
    
    asky = B_180
    bunch, A_flags = asym_solver(masked,center,asky,iterative)
    A_center = bunch[0]
    A = bunch[1]
    
    A_center[0] += self['MYMIN_IMAGE']
    A_center[1] += self['MXMIN_IMAGE']
    self['M_AS_X'] = A_center[1]
    self['M_AS_Y'] = A_center[0]
    self['M_AS'] = A
    
    self['flags'] = addflag(self['flags'],A_flags)
    
def asym_solver(image,center,asky,iterative):
    """Returns minimum 180deg Asymmetry of given image, and point with respect to which, Asymmetry is 
    such a minimum."""    
    # IMPORT STUFF
    from scipy.optimize import fmin
    from flags import addflag, allflags
    # END IMPORT STUFF
    
    flags = 0L
    
    if iterative :
        bunch = fmin(asym_func_iraf,center,args=(image,asky),full_output=1,maxiter = 1000.,disp=0)
        # xopt = bunch[0] # miminum A center
        # fopt = bunch[1] # A
        # iter = bunch[2] # number of iterations
        # funcalls = bunch[3] # number of function calls
    else :
        A = asym_func_iraf(center,image,asky)
        bunch = [[center[0],center[1]],A,0,0,0]
    warnings = bunch[4] # warnings: 1= max number of func. evals, 2=max number of iterations
    if warnings == 2 : flags = addflag(flags,allflags['MAXITER_AS'])
    
    return bunch, flags
    
def asym_func_iraf(center,image,zero):    
    """Returns asymmetry value given masked image and center. Uses Iraf."""
    # IMPORT STUFF
    from pyraf import iraf
    import pyfits
    import os
    #END IMPORT
    
    tmpimage1 = 'tmp180_1.fits'
    tmpimage2 = 'tmp180_2.fits'
    
    pyfits.writeto(tmpimage1,image)
    
    iraf.load('images',doprint=0)
    iraf.images.imgeom()
    
    x = center[1] ; y = center[0]
    xshift = (image.shape[1]/2. + 0.5) - (x+1)
    yshift = (image.shape[0]/2. + 0.5) - (y+1)
    
    iraf.images.imgeom.imshift(input=tmpimage1,output=tmpimage1,xshift = xshift,yshift = yshift,\
    interp='linear',boundar='constant',constan=0.)
    iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() 
    
    iraf.images.imgeom.rotate(input=tmpimage1,output=tmpimage2,rotation=180,xin='INDEF',yin='INDEF',\
    xout='INDEF',yout='INDEF',interpo='linear',boundar='constant',constan=0,verbose='no')
    
    image1 = pyfits.getdata(tmpimage1)
    image2 = pyfits.getdata(tmpimage2)
    
    os.system('rm %s %s' % (tmpimage1,tmpimage2))
    
    blank1 = image1 == 0 ; blank2 = image2 == 0
    blank = num.where( blank1 | blank2)
    image1[blank] = 0. ; image2[blank] = 0.
    
    from pylab import imshow,show
    imshow(num.abs(image1-image2),origin='lower') ; show()
    
    up = num.abs(image1 - image2)
    up = num.sum(num.sum(up))
    down = num.abs(image1)
    down = num.sum(num.sum(down))
    
    Ax = up / down - zero
    
    # print 'Partial Ax = %4.2f' % Ax
    
    return Ax
    
##def asym_func(center,image,zero):
##    """Returns asymmetry value given masked image and center."""
##    # IMPORT STUFF
##    #END IMPORT
##    
##    x = center[1] ; y = center[0]
##    simage = rot180iraf(x,y,image)
##    
##    blank1 = image == 0 ; blank2 = simage == 0.
##    blank = num.where(blank1 | blank2)
##    image[blank] = 0. ; simage[blank] = 0.
##    
##    up = num.abs(image - simage)
##    up = num.sum(num.sum(up))
##    down = num.abs(image)
##    down = num.sum(num.sum(down))
##    
##    A = up / down - zero
##   
####    from pylab import imshow,show
####    imshow(num.abs(image - simage),origin='lower'); show()
####    print 'Asym_partial = %4.2f' % A
####    stop()
##    
##    return A
##    
##def rot180iraf(x,y,image):
##    """Returns a rotation of 180deg of given image around given center. 
##    Returned image has same dimensions as input image. Uses iraf."""
##    """Returns a rotation of 180deg of given image around given center. 
##    Returned image has same dimensions as input image."""
##    # IMPORT STUFF
##    from pdb import set_trace as stop
##    from pyraf import iraf
##    import pyfits
##    import os
##    # END IMPORT
##    tmpimage= 'tmp180.fits'
##    pyfits.writeto(tmpimage,image)
##    iraf.load('images',doprint=0)
##    iraf.images.imgeom()
##    iraf.images.imgeom.rotate(input=tmpimage,output=tmpimage,rotation=180,xin=x,yin=y,interpo='linear',\
##    boundar = 'constant',constan=0.,verbose='no')
##    simage = pyfits.getdata(tmpimage)
##    os.system('rm %s' % tmpimage)
##    
##    return simage
    
##def rot180(x,y,image):
##    """Returns a rotation of 180deg of given image around given center. 
##    Returned image has same dimensions as input image."""
##    # IMPORT STUFF
##    from algorithms import trimindex,intindex
##    import sys
##    # END IMPORT
##    print 'for some reason, this task does not work properly'
##    sys.exit()
##    
##    indexes = num.where(image == image)
##    xp = indexes[1].astype('Float32').copy()
##    yp = indexes[0].astype('Float32').copy()
##    incx = 2.*(x - xp)
##    incy = 2.*(y - yp)
##    xp += incx ; yp += incy
##    rotated = (yp,xp)
##    
##    indexes = intindex(indexes,num)
##    rotated = intindex(rotated,num)
##    maxx = max(indexes[1]) 
##    maxy = max(indexes[0])
##    
##    dummie = num.zeros(shape=(len(indexes[0]),),type='Int8')
##    lesserthan0_0 = rotated[0] < 0
##    biggermaxy_0 = rotated[0] > maxy
##    lesserthan0_1 = rotated[1] < 0
##    biggermaxx_1 = rotated[1] > maxx
##    dummie[num.where(lesserthan0_0 | biggermaxy_0 | lesserthan0_1 | biggermaxx_1)] = 1
##    
##    
##    rotated = trimindex(rotated,dummie,num)
##    indexes = trimindex(indexes,dummie,num)
##    simage = num.zeros(shape=image.shape,type='Float32') 
##    simage[indexes] = image[rotated]
##    # from pylab import show, imshow # DEVEL
##    # stop()                                                # DEVEL
##    
##    return simage
## 
##def trimindex(indexes,mask,num):
##    goods = num.where(mask==0)
##    return (indexes[0][goods],indexes[1][goods])
##
##def intindex(indexes,num):
##    return (num.around(indexes[0]).astype('Int32').copy(),\
##    num.around(indexes[1]).astype('Int32').copy())

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
    # END IMPORT
    
    # INPUTS
    sky = self['BACKGROUND']
    stamp = self['STAMP'] 
    mask = self['MASK']
    # mask_sky = self['MASK_SKY']
    r_petro = self['R_PETRO']
    box2petro = self.execpars['box2petro'][0]
    # END INPUTS
    
    stamp -= sky
    
    def maskit(stamp,mask):
        masked = num.zeros(shape=stamp.shape,type=stamp.type())
        active = num.where(mask == 0)
        masked[active] = stamp[active]
        return masked

    masked = maskit(stamp,mask)
    # masked_sky = maskit(stamp,mask_sky)
    boxwidth = int(num.around(box2petro * r_petro))
    if boxwidth % 2 == 0 : boxwidth += 1
    if boxwidth < 3 : boxwidth = 3

    smoothed = boxcar(masked,(boxwidth,boxwidth),mode='reflect') 
    # skysmoothed = boxcar(masked_sky,(boxwidth,boxwidth),mode='reflect') 
    
    def g(image,simage,num):
        up = num.abs(image - simage)
        down = num.abs(image)
        up = num.sum(num.sum(up))
        down = num.sum(num.sum(down))
        return up/down
    
    # Bs = g(masked_sky,skysmoothed)
    Bs = 0
    S = g(masked,smoothed,num) - Bs
    
    S_flags = 0L
    
    self['M_S'] = S
    self['flags'] = addflag(self['flags'],S_flags)
    
    # return S, S_flags
    
    
def concent(self):
    """Updates in self the 'Concentration', acording to the definition given in Bershady et al., 2000 
    (extracted from Lotz et al., 2004)
    
    C = 5 log (r80/r20)
    where r? is the radius of the circular aperture that contains ? % of the total flux of the object. The total flux
    is that contained within 1.5rp of the galaxy's center, as in Lotz et al., 2004 and Conselice 2003.
    How is the center determined? Lotz et al. use that given by the determination of the 'A' parameter.
    """
    # IMPORT STUFF
    import scipy.interpolate as interp
    from radial import radial
    from flags import addflag, isflagon, allflags
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
        
        def lazy(normal,radii,proc,addflag):
            # IMPORT STUFF
            from pdb import set_trace as stop
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
        if isflagon(flag_small,2L) or isflagon(flag_big,2L) : 
            cflag = addflag(cflag,allflags['NOCONC'])
            cflag = addflag(cflag,allflags['MANYCONC'])
        C = -99.0
        if not isflagon(cflag,allflags['NOCONC']) : C = 5. * num.log10(r_big / r_small)
        
        self['M_C'] = C
        self['flags'] = addflag(self['flags'],cflag)
        
    else : 
        self['M_C'] = -99.0
        self['flags'] = addflag(self['flags'],allflags['NOCONC'])
    
def peak(self):
   """Updates in self the peak center: xpeak, ypeak."""
   # IMPORT STUFF
   from numpy.convolve import boxcar
   # END IMPORT
   
   # INPUTS
   img = self['STAMP']
   mask = self['MASK']
   box2petro = self.execpars['box2petro'][0]
   x0 = self['MXMIN_IMAGE']
   y0 = self['MYMIN_IMAGE']
   r_petro = self['R_PETRO']
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
    
    img = self['STAMP']
    sky = self['BACKGROUND']
    # MASKING.
    mask = self['MASK']
    active = num.where(mask == 0)
    
    img -= sky
    
    x0 = self['MXMIN_IMAGE']
    y0 = self['MYMIN_IMAGE']
    
    x,y = moments(img,x0,y0,order=1,active=active)
    x2, y2, xy = moments(img,x0,y0,order=2,active=active)
    
    self['M_X'] = x
    self['M_Y'] = y
    self['M_X2'] = x2
    self['M_Y2'] = y2
    self['M_XY'] = xy
    
    
