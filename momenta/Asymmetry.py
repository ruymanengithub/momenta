#! /usr/bin/env python

"""Asymmetry"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as num
import pyfits
import os
# END IMPORT

def asymmetry(self,dograph=False):
    """Updates in self the 'Asymmetry' (A) parameter, according to 
    the definition given in Lotz. et al., 2004.
    
    A = sum | I(i,j) - I_180(i,j)| / (2 * sum| I(i,j) |)  - B_180 
    
    where I is the galaxy's image and I_180 is the image rotated by 
    180deg about the galaxy's central pixel, and
    B_180 is the average asymmetry of the background. A is summed over all
    pixels within 1.5 rp of the galaxy's center (or segmentation mask). The
    central pixel is determined by minimizing A.
  
    """
    # IMPORT STUFF
    from flags import addflag,isflagon,allflags
    from algorithms import increase_window
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
    rotmode = 'iraf' # iraf, numpy
    minimode = self.execpars['minModeAsym'][0] # fmin, con00
    try: Asymstep = self.execpars['Asymstep'][0]
    except KeyError : Asymstep = 1.
    # END INPUTS
    
    # To speed up... 
    # if num.any(num.array(stamp.shape) > 500) : rotmode = 'iraf'
    
    masked = num.zeros(shape=stamp.shape,type=stamp.type()) 
    active = num.where(mask == 0)
    masked[active] = stamp[active] - sky
    masked,xshift,yshift = increase_window(masked)
    xcenter += xshift
    ycenter += yshift
    center = (ycenter,xcenter)
    del active
    
    if doAsym_sky:
        if 'SKYSTAMP' in self: skystamp = self['SKYSTAMP'].copy()
        else: skystamp = self.makeskystamp(masked,bads=0)
        active = num.where(skystamp != 0)
        skystamp[active] = skystamp[active] - sky
        skystamp,ignore,ignore = increase_window(skystamp)
        stop('check size(systamp) == size(mask)')
        del active
    
    B_180 = 0.
    
    bunch, A_flags = asym_solver(masked,center,B_180,iterative,\
    rotmode=rotmode,minimode=minimode)
    A_center = bunch[0]
    A = bunch[1]
    
    if isflagon(A_flags,allflags['MAXITER_AS']):
        bunch, dump = asym_solver(masked,center,B_180,iterative=False,\
        minimode=minimode)
        A_center = bunch[0]
        A = bunch[1]
     
    if doAsym_sky:
        norm = num.abs(masked)
        norm = num.sum(num.sum(norm))
        norm = norm * 2.
        bunch_sky,ignore = asym_solver(skystamp,A_center,0.,\
        iterative=False,norm=norm,minimode=minimode)
        A_sky = bunch_sky[1]
        self['M_AS_SKY'] = A_sky
    else : A_sky = -99
    
    self['M_AS_X'] = A_center[1] - xshift + self['MXMIN_IMAGE']
    self['M_AS_Y'] = A_center[0] - yshift + self['MYMIN_IMAGE']
    self['M_AS'] = A
    
    print '\nAsym = %4.2f    Asym_sky = %4.2f' % (A,A_sky)
    print 'A_shift = %4.2f pixels\n' % \
    num.sqrt((center[1]-A_center[1])**2.+(center[0]-A_center[0])**2.)
    
    self['flags'] = addflag(self['flags'],A_flags)
    
    if dograph:
        self.asymmetry_graph(masked,A_center)
    else : pass
    
def asym_solver(image,center,asky,iterative,norm=-1,rotmode='iraf',\
    minimode='fmin',step=1.):
    """Returns minimum 180deg Asymmetry of given image, and point 
    with respect to which, Asymmetry is such a minimum."""    
    # IMPORT STUFF
    from scipy.optimize import fmin
    # from Asymmetry import asym_mini_con00
    from flags import addflag, allflags
    # END IMPORT STUFF
    
    flags = 0L
    maxiter = 50
    
    if iterative :
        args = (image,asky,norm,rotmode)
        if minimode == 'fmin':
            bunch = fmin(asym_func,center,args,\
            full_output=1,maxiter=maxiter,disp=0)
        elif minimode == 'con00':
            bunch = asym_mini_con00(asym_func,center,args,\
            maxiter=maxiter,step=step) # step in pixels
        # xopt = bunch[0] # miminum A center
        # fopt = bunch[1] # A
        # iter = bunch[2] # number of iterations
        # funcalls = bunch[3] # number of function calls
    else:
        A = asym_func(center,image,asky,norm,rotmode)
        bunch = [[center[0],center[1]],A,0,0,0]
    
    warnings = bunch[4] 
    # warnings: 1= max number of func. evals, 2=max number of iterations
    if warnings == 2 : 
        flags = addflag(flags,allflags['MAXITER_AS'])
        print 'MAXIMUM NUMBER OF ITERATIONS (%i) reached in ASYMMETRY' %\
        maxiter
    
    return bunch, flags
    
def asym_func(center,image,zero,norm=-1,rotmode='iraf'):    
    """Returns asymmetry value given masked image and center."""
    # IMPORT STUFF
    from algorithms_II import rot180numpy
    #END IMPORT
    
    if rotmode == 'iraf':
        image1,image2 = rot180iraf(image,center)
    elif rotmode == 'numpy':
        image1, image2 = rot180numpy(image,center)
    
    # nonzero = num.where(True ^ ((image1 != 0) ^ (image2 != 0)))
    nonzero = num.where(image1 == image1)
    up = num.sum(num.abs(image1[nonzero] - image2[nonzero]))
    
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
    
    pyfits.writeto(tmpimage1,image.astype('Float32'))
    
    iraf.load('images',doprint=0)
    iraf.images.imgeom()
    
    x = center[1] ; y = center[0]
    xshift = ((image.shape[1]-1)/2.) - x # 
    yshift = ((image.shape[0]-1)/2.) - y # 
    
    iraf.images.imgeom.imshift(input=tmpimage1,output=tmpimage1,\
    xshift = xshift,yshift = yshift,interp='linear',\
    boundar='constant',constan=0.)
    iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; iraf.flpr() ; 
    iraf.flpr() ; iraf.flpr() ; iraf.flpr() 
    
    iraf.images.imgeom.rotate(input=tmpimage1,output=tmpimage2,\
    rotation=180,xin='INDEF',yin='INDEF',\
    xout='INDEF',yout='INDEF',interpo='linear',boundar='constant',\
    constan=0,verbose='no')
    
    image1 = pyfits.getdata(tmpimage1).astype('Float32')
    image2 = pyfits.getdata(tmpimage2).astype('Float32')
    os.system('rm %s %s' % (tmpimage1,tmpimage2))
    
    return image1,image2

def asym_mini_con00(asym_func,inicenter,args,maxiter=20,step=1.):
    """Finds center for which Asymmetry has a minimum value.
    Algorithm inspired by that defined in Conselice et al. 2000"""
    # IMPORT STUFF
    from copy import copy
    # END IMPORT
    
    # INPUTS
    image = copy(args[0]) ; asky = copy(args[1])
    norm = copy(args[2]) ; mode = copy(args[3])
    idims = image.shape
    # END INPUTS
    
    nsteps = 3                                  # always odd dimensions!
    cindx = (nsteps/2,nsteps/2)      # index of central value
    
    center = copy(inicenter)           # initial value
    Amatrix = num.zeros(shape=(nsteps,nsteps),type='Float32') 
    # matrix of "A" values
    
    Acenter0 = copy(inicenter)                                             # default value
    Amin0 = asym_func(Acenter0,image,asky,norm,mode) # default value
    
    go = True
    niter = 1             # number of iterations

    for niter in range(1,maxiter+1):
        
        # Compute A in a nstep x nstep grid
        for xs in range(-(nsteps/2),nsteps/2+1):
            for ys in range(-(nsteps/2),nsteps/2+1):
                tcenter = (center[0]+ys*step,center[1]+xs*step)
                # tcenter in bounds?
                inbounds = ((tcenter[0]>=0) and (tcenter[0]<idims[0]) and\
                (tcenter[1]>=0) and (tcenter[1]<idims[1]))
                
                if inbounds:
                    Amatrix[ys+cindx[0],xs+cindx[1]] = \
                    asym_func(tcenter,image,asky,norm,mode)
                else:
                    #print 'not in bounds'
                    break
        
        # Which center in grid gives a minimum value of A?
        Amin = Amatrix.min()
        Amin_indx = num.where(Amatrix == Amin)
        manysols = len(Amin_indx[0]) > 1               # many solutions?
        #if manysols : print 'manysols'
        Acenter = center[0]+(Amin_indx[0][0]-cindx[0])*step,\
        center[1]+(Amin_indx[1][0]-cindx[1])*step
        
        offset = num.sqrt((Acenter[0]-inicenter[0])**2.+\
        (Acenter[1]-inicenter[1])**2.)
        
        #print 'minimum Asymmetry, offset (iter): %.3f, %.3f, %i, %i' % \
        #(Amin, offset,Acenter[1],Acenter[0])
        #print Amatrix
        
        if Amin_indx[0][0] == cindx[0] and Amin_indx[1][0] == cindx[1]:
            break
        else:
            center = copy(Acenter)
            Amatrix[:,:] = 0.
        
        if not inbounds or manysols:
            Amin = copy(Amin0)
            Acenter = copy(Acenter0)
            niter = maxiter # any error will be considered as a
            #maximum iteration error
            break
    
    warning = 0
    if niter == maxiter : warning = 2
    solution = [Acenter,Amin,niter,niter*nsteps**2.,warning]
    
    return solution
