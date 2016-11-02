#! /usr/bin/env python

import numpy as num
from pdb import set_trace as stop

def clumpy(self,dograph=False):
    """Updates in self the 'Clumpiness' parameter (S) according to 
    the definition given in Lotz et al., 2004.
    
    S = sum | I(i,j) - Is(i,j)| / sum | I(i,j) | - Bs
    
    where Is is the galaxy's image smoothed by a boxcar of width 0.25 rp 
    (rp = petrosian radius), and Bs is 
    the average smoothness of the background. Like A, S is summed over 
    the pixels within 1.5rp of the galaxy's
    center (or in segmentation area). Bs is average 'clumpiness' of sky.
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
    try: NullCenter_S = self.execpars['NullCenter_S'][0] == 1
    except KeyError: NullCenter_S = False
    # END INPUTS
    
    masked = num.zeros(shape=stamp.shape,type=stamp.type())
    active = num.where(mask == 0)
    masked[active] = stamp[active] - sky
    if doClumpy_sky : 
        if 'SKYSTAMP' in self: skystamp = self['SKYSTAMP'].copy()
        else: skystamp = self.makeskystamp(masked,bads=0)
        skystamp[active] = skystamp[active] - sky
    del active
    
    boxwidth = int(num.around(box2petro * r_petro))
    if boxwidth % 2 == 0 : boxwidth += 1
    if boxwidth < 3 : boxwidth = 3
    
    if NullCenter_S :
        xcenter = int(num.around(self['X_IMAGE'] - self['MXMIN_IMAGE']))
        ycenter= int(num.around(self['Y_IMAGE'] - self['MYMIN_IMAGE']))
        sw = boxwidth/2
        masked[ycenter-sw:ycenter+sw:,xcenter-sw:xcenter+sw] = 0.
    
    smoothed = boxcar(masked,(boxwidth,boxwidth),mode='constant',\
    cval=0.0)

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
    
    if dograph:
        img = masked-smoothed
        self.clumpy_graph(img)
    else : pass

