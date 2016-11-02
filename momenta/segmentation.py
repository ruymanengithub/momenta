#! /usr/bin/env python

def FindClumps(self,dograph=False):
    """Finds clumps in an image."""
    # IMPORT STUFF
    from numpy.nd_image.filters import uniform_filter
    import numpy as num
    from pdb import set_trace as stop
    from flags import allflags,addflag
    # END IMPORT
    
    # INPUTS
    sky = self['BACKGROUND']
    image = self['STAMP'].copy() - sky
    try : mask = self['MASK'].copy()
    except AttributeError: mask = self['MASK']
    fluxratio = 0.5
    box = 5
    MAXITER = 1000.
    # END INPUTS
    
    if mask is -1 : mask = num.zeros(shape=image.shape,type='Int8')
    image[num.where(mask != 0)] = 0.
    
    nactive = len(num.where(mask==0)[0])
    
    Norm = image.sum()
    
    success = True
    if Norm <= 0. : success = False
    
    sbox = int(box / 2)
    aggregated = num.zeros(shape=image.shape,type='Float32')
    
    clumps = []
    xcoo = [] ; ycoo = []
    
    # GoOn = True
    niter = 1
    for niter in range(1,MAXITER+1):
        
        Npix = num.ones(shape=image.shape,type='Float32')
        Npix[num.where(image == 0.)] = 0
        uniform_filter(Npix,(box,box),output=Npix,mode='constant',cval=0.)
        Npix *= box ** 2.
        
        uniform_filter(image,(box,box),output=aggregated,\
        mode='constant',cval=0.)
        aggregated *= Npix
        maxval = aggregated.max()
        coomax = num.where(aggregated == maxval)
        xmax = coomax[1][0]
        ymax = coomax[0][0]
        xcoo.append(xmax) ; ycoo.append(ymax)
        fluxrel = maxval/Norm
        clumps.append(fluxrel)
     
        image[ymax-sbox:ymax+sbox+1,xmax-sbox:xmax+sbox+1] = 0.
        
        accumulated = num.sum(clumps)
        # print 'acc=%f'%accumulated
        
        if accumulated >= fluxratio : break
        
    if niter == MAXITER : success = False
    
    if success:
        
        NClumps = len(clumps)
        MaxClump = clumps[0]
        MinClump = clumps[-1]
        
        argmaxima = (num.array(ycoo),num.array(xcoo))
        maxima = num.zeros(shape=image.shape,type='Int8')
        maxima[argmaxima] = 1
        
        self['CL_LOC'] = maxima.copy()
        self['M_NUM_CL'] = NClumps
        self['M_MAX_CL'] = MaxClump ; self['M_MIN_CL'] = MinClump
        self['M_ACC_CL'] = accumulated
        self['M_FAR_CL'] = NClumps * box**2. / nactive
        print '%i Clumps, %.2f o/o to %.2f o/o in flux' % \
        (NClumps,MaxClump*100.,MinClump*100.)
        
        if dograph:
            self.FindClumps_graph()
    else:
        self['flags'] = addflag(self['flags'],allflags['NOCLUMPS'])
        self['M_NUM_CL'] = -99
        self['M_MAX_CL'] = -99. ; self['M_MIN_CL'] = -99.
        self['M_ACC_CL'] = -99. ; self['M_FAR_CL'] = -99.
        print 'No Clumps'
        
    return None
    
