#! /usr/bin/env python

"""Axis Asymmetry."""

# IMPORT STUFF
import numpy as num
from pdb import set_trace as stop
import sys
from pyraf import iraf
import pyfits
from time import time
from Moments.graphs import Convert, doStamp
# END IMPORT

def axis_asymmetry(self,axis,dograph=False):
    """Updates in self the degree of axis symmetry respect to 
    minor/major semi-major axis.
    
    Ax = sum | I(i,j) - I_180(i,j)| / (2 * sum| I(i,j) |)  - Ax_180 
    
    where I is the galaxy's image and I_180 is the image bended about 
    the major or minor axis of the galaxy, and
    B_180 is the average asymmetry of the background. Ax is summed 
    over all pixels within 1.5 rp of the galaxy's
    center. The central pixel is fixed. Which one to use??"""
    # IMPORT STUFF
    from flags import addflag
    import pyfits
    # END IMPORT
    
    # INPUTS
    image = self['STAMP'].copy()
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    center = (self['Y_IMAGE']-self['MYMIN_IMAGE'],self['X_IMAGE']-\
    self['MXMIN_IMAGE'])
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
        if 'SKYSTAMP' in self: skystamp = self['SKYSTAMP'].copy()
        else: skystamp = self.makeskystamp(masked,bads=0)
        active = num.where(skystamp != 0)
        skystamp[active] = skystamp[active] - sky
    del active
    
    axsky = 0.0
    bunch, Ax_flags = axisasym_solver(masked,center,pa,doAxisAsym_iter,\
    axsky,norm=-1,mode=mode)
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
    
    if dograph:
        self.axis_asymmetry_graph(label,masked,pa,Ax_center)
    else : pass

def axisasym_solver(image,center,pa,iterative,axsky,norm=-1,\
    mode='iraf'):
    """Returns minimum Axis Asymmetry of given image, and point 
    with respect to which, Axis Asymmetry is such a minimum."""
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
    """Computes degree of axis symmetry for generic image and axis 
    of symmetry. Uses either iraf or numpy."""
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
    
    pyfits.writeto(tmpimage1,image.astype('Float32'))
    
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
    
    image1 = pyfits.getdata(tmpimage1).astype('Float32')
    image2 = pyfits.getdata(tmpimage2).astype('Float32')
    
    os.system('rm %s %s' % (tmpimage1,tmpimage2))
    
    return image1, image2
    
