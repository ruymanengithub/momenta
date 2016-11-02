#! /usr/bin/env python

def radial_v1(self,MaxBadRatio=0.25,ExtrapolBadPix=False,IgnoreBadValues=False,dograph=False):
    """Updates in self the radial profile of source."""
    # makes curve of integrated flux in elliptical appertures
    # the radial coordinate should be "corrected" to give the integer number 
    # of pixels enclosed by the appertures. Dubious...!!
    # It should be possible to mask out given pixels.
    # The value of the integrated flux should be corrected from the masking,
    # extrapolating the values ? I'd say no... but it's an option.
    # The intensity profile is built derivating the growth-curve.
    # The growth and intensity curves, the radial coordinate, and the number 
    # of pixels and the number of dismissed pixels in the each are given on 
    # return. 
    # If less than 25% of the pixels in an apperture or more are flagged out, 
    # the apperture is not valid.
    # the initial list of appertures is fixed. The final list of appertures depends 
    # on where is the center of the appertures and the flagged pixels.
    
    # DIFFERENCES WITH RESPECT TO VERSION 0:
    # npix is now the number of pixels inside the elliptical aperture, including  
    # the flagged pixels.
    
    # IMPORT STUFF
    from ellipse import effective_radius, area_superellip, dist_superellipse
    import numpy as num
    from flags import allflags, addflag
    from scipy import interpolate as interpol
    from momlib import mergemasks
    from pdb import set_trace as stop
    from time import time
    from deriv import deriv
    # END IMPORT STUFF
    
    # INPUTS
    stamp = self['STAMP'].copy()
    # image masking
    mask = self['MASKOTHER'].copy()
    xcenter = self['X_IMAGE']  ; ycenter = self['Y_IMAGE']
    xcenter = xcenter - self['MXMIN_IMAGE']
    ycenter = ycenter - self['MYMIN_IMAGE']
    ellip = self['ELLIPTICITY'] 
    pa = self['THETA_IMAGE'] - 90. # Astronomical position angle.
    c = self['BOXY']
    sky = self['BACKGROUND']
    # END INPUTS
    
    # initialization of radial distances.
    
    # minimum radius : 1 pix
    # maximum radius : depends on the area given in stamp.
    # logarithmic increase of factor 0.1
    
    n = stamp.shape
    q = 1. - ellip
    center = (ycenter,xcenter)
    radial_mask = dist_superellipse(n,center,q=q,pos_ang=pa,c=c)
    
    # Geometric progression of radius
    LogRadial = self.execpars['doLogRadial'][0] == 1
    FineRadial = self.execpars['doFineRadial'][0] == 1
    
    if LogRadial and not FineRadial:
        logminradius = 0.
        logmaxradius = num.log10(max(radial_mask.flat))
        logstep = 0.03
        fini_radial = num.arange(logminradius,logmaxradius+logstep,\
        logstep,type='Float32')
        fini_radial = 10.**fini_radial
        ini_radial = num.zeros(shape=(len(fini_radial)+1,),type='Float32')
        ini_radial[1:] = fini_radial.copy()
    elif not LogRadial and not FineRadial:
        # ALTERNATIVE SHOT: Mixture of linear and geometric progression of 
        # radius.
        logminradius = 0.
        logmaxradius = min(num.log10(max(radial_mask.flat)),1.)
        logstep = 0.03
        fini_radial = num.arange(logminradius,logmaxradius+logstep,\
        logstep,type='Float32')
        fini_radial = 10.**fini_radial
        if max(radial_mask.flat) > 10. : 
            lini_radial = num.arange(10.,max(radial_mask.flat),type='Float32')
            if len(lini_radial) > 1 : lini_radial = lini_radial[1:]
        ini_radial = num.zeros(shape=(len(fini_radial)+len(lini_radial)+1,),\
        type='Float32')
        ini_radial[1:len(fini_radial)+1] = fini_radial.copy()
        ini_radial[len(fini_radial)+1:] = lini_radial.copy()
    elif FineRadial:
        maxradius = max(radial_mask.flat)
        ini_radial = num.arange(float(maxradius))
    
    radiusflag = num.ones(shape=ini_radial.shape,type='Int8')
    
    cumulflx = num.array(shape=ini_radial.shape,type='Float32')
    radii = num.array(shape=ini_radial.shape,type='Float32')
    npix = num.array(shape=ini_radial.shape,type='Int32')
    npixout = num.array(shape=ini_radial.shape,type='Int32')
    cumulflx[:] = -1.0 ; radii[:] = -1.0
    npix[:] = -1 ; npixout[:] = -1
    
    index = 0
    
    for radius in ini_radial:
        englobed = num.where(radial_mask <= radius)
        nenglobed = len(englobed[0])
        dummie = mask + (radial_mask > radius).astype('Int8') * 2
        active = num.where(dummie == 0)
        boundary = (0 in active[0]) or (stamp.shape[0]-1 in active[0]) \
        or (0 in active[1]) or (stamp.shape[1]-1 in active[1])
        if boundary :
            radiusflag[index:] = 0
            break
        else:
            nactive = len(active[0])
            null = num.where(dummie == 1) ; nnull = len(null[0])
            del dummie
            
            effradius = effective_radius(nenglobed,q,c=0)
            # effradius = radius
            
            cumul = num.sum(stamp[active]) - sky * nactive
            try : corrfactor = float(nenglobed)/float(nactive)
            except ZeroDivisionError : corrfactor = 1.
            if ExtrapolBadPix : cumul = cumul * corrfactor
            #print cumul/nenglobed,1./corrfactor
            if (effradius in radii) or (nnull > MaxBadRatio * nenglobed) or \
            (nenglobed in npix) or (cumul in cumulflx): 
                radiusflag[index] = 0
                # "nenglobed in npix" and "cumul in cumulflx" avoids 
                # duplication of values which make splrep crash
            else : 
                radii[index] = effradius
                cumulflx[index] = cumul
                npix[index] = nenglobed
                npixout[index] = nnull
            index += 1
    
    radiusflags = 0L
    
    goodradii = num.where(radiusflag != 0)[0]
    radii = radii[goodradii]
    cumulflx = cumulflx[goodradii]
    npix = npix[goodradii]
    npixout = npixout[goodradii]
    
    if not IgnoreBadValues:
        
        # CHECK FOR NEGATIVE GROWTH
        
        maxcumul = cumulflx.max()
        maxcumulindex = cumulflx.argmax()
        if maxcumulindex < len(cumulflx)-1:
            postmax = (num.arange(maxcumulindex+1,len(cumulflx)),)
            cumpostmax = cumulflx[postmax]
            below90 = num.where(cumpostmax < maxcumul * 0.9)
            if len(below90[0]) != 0:
               lastgood = postmax[0][below90][0]
               radii = radii[0:lastgood]
               cumulflx = cumulflx[0:lastgood]
               npix = npix[0:lastgood]
               npixout = npixout[0:lastgood]
               print 'Growth curve going negative!'
               radiusflags = addflag(radiusflags,allflags['NEGGROWTH'])
            else : pass
        else : pass
        
        if len(radii) < 8 : radiusflags = addflag(radiusflags,allflags['BADRADIAL'])
        if len(radii) < 4 : 
            radiusflags = addflag(radiusflags,allflags['NORADIAL'])
            radii = None ; cumulflx = None ; intens2 = None ; 
            npix = None ; npixout = None ; 
        else:
            # spline = interpol.splrep(npix,cumulflx,s=0,k=3)
            intens2 = deriv(npix,cumulflx)
    else:
        radiusflags = addflag(radiusflags,allflags['NONCHECKEDRADIAL'])
        if len(npix) > 3:
            intens2 = deriv(npix,cumulflx)
        else : 
            radiusflags = addflag(radiusflags,allflags['NORADIAL'])
            radii=None ; cumulflx=None
            intens2=None ; npix=None ; npixout=None
    
    larousse = {'radii':radii,'cumulflx':cumulflx,'intens':intens2,'npix':npix,'npixout':npixout,'radiusflags':radiusflags}
    
    self['flags'] = addflag(self['flags'],radiusflags)
    self['M_RADIAL'] = larousse
    
    if dograph: self.radial_graph()
    
    return None
    
##def radial(self):
##    """Updates in self the radial profile of source."""
##    # make curve of integrated flux in elliptical appertures
##    # the radial coordinate should be "corrected" to give the integer number of pixels enclosed by the
##    # appertures.
##    # It should be possible to mask out given pixels.
##    # Should the value of the integrated flux be corrected from the masking, extrapolating the values ?
##    # The intensity profile is built derivating the B-spline interpolation of the growth-curve.
##    # the growth and intensity curves, the radial coordinate, and the number of pixels and the number
##    # of dismissed pixels in the each are given on return. 
##    # If less than 50% of the pixels in an apperture are flagged out, the apperture is not valid.
##    # the initial list of appertures is fixed. The final list of appertures depends on where is the center of the 
##    # appertures and the flagged pixels.
##   
##    # IMPORT STUFF
##    from ellipse import effective_radius
##    from ellipse import area_superellip
##    from ellipse import dist_superellipse
##    import numpy as num
##    from flags import allflags, addflag
##    from scipy import interpolate as interpol
##    from momlib import mergemasks
##    from pdb import set_trace as stop
##    from time import time
##    # END IMPORT STUFF
##    
##    # t1 = time()
##    
##    stamp = self['STAMP'].copy()
##    
##    # image masking
##    # mask = self['MASK'] # WRONG! the segmentation map of the object is included in the mask!
##    mask = self['MASKOTHER'].copy() # Right. Only other objects are flagged out.
##    
##    xcenter = self['X_IMAGE'] ; ycenter = self['Y_IMAGE']
##    xcenter = xcenter - self['MXMIN_IMAGE']
##    ycenter = ycenter - self['MYMIN_IMAGE']
##    ellip = self['ELLIPTICITY'] 
##    pa = self['THETA_IMAGE'] - 90. # Astronomical position angle.
##    c = 0
##    sky = self['BACKGROUND']
##   
##    # initialization of radial distances.
##   
##    # minimum radius : 1 pix
##    # maximum radius : depends on the area given in stamp.
##    # logarithmic increase of factor 0.1
##    # for each radius > 5, the area (in pixels) englobed should be computed and compared against
##    # theoretical area. When they differ by >50%, that radius and bigger are dropped.
##    
##    n = stamp.shape
##    q = 1. - ellip
##    center = (ycenter,xcenter)
##    radial_mask = dist_superellipse(n,center,q=q,pos_ang=pa,c=c)
##    
##    
##    logminradius = 0.
##    logmaxradius = num.log10(max(radial_mask.flat))
##    logstep = 0.1
##    ini_radial = num.arange(logminradius,logmaxradius+logstep,logstep,type='Float32')
##    ini_radial = 10.**ini_radial
##    
##    index = 0
##    for index in range(len(ini_radial)):
##        actarea = num.where(radial_mask <= ini_radial[index])
##        criterio = (0 in actarea[0]) or (stamp.shape[0]-1 in actarea[0]) \
##           or (0 in actarea[1]) or (stamp.shape[1]-1 in actarea[1])
##        if criterio : ini_radial[index] = -1
##    del index
##    
##    ini_radial = ini_radial[num.where(ini_radial != -1)]
##    
##    radiusflag = num.array(shape=ini_radial.shape,type='Int8')
##    radiusflag[:] = 1
##    
##    cumulflx = num.array(shape=ini_radial.shape,type='Float32')
##    radii = num.array(shape=ini_radial.shape,type='Float32')
##    npix = num.array(shape=ini_radial.shape,type='Int32')
##    npixout = num.array(shape=ini_radial.shape,type='Int32')
##    cumulflx[:] = -1.0
##    radii[:] = -1.0
##    npix[:] = 0
##    npixout[:] = 0
##    
##    index = 0
##    
##    for radius in ini_radial:
##        englobed = num.where(radial_mask <= radius)
##        nenglobed = len(englobed[0])
##        dummie = num.ones(shape=radial_mask.shape,type='Int8')+1 # all 2
##        dummie[num.where(radial_mask <= radius)] = 0
##        dummie[num.where(mask != 0)] += 1
##        active = num.where(dummie == 0)
##        nactive = len(active[0])
##        null = num.where(dummie == 1) 
##        nnull = len(null[0])
##        
##        effradius = effective_radius(nenglobed,q,c=0)
##        cumul = num.sum(stamp[active]) - sky * nactive
##        
##        if (effradius in radii) or (nnull > 0.25 * nenglobed) or (nactive in npix) : 
##            radiusflag[index] = 0
##            # "nactive in npix" avoids duplication of values which make splrep crash
##        else : 
##             radii[index] = effradius
##             cumulflx[index] = cumul
##             npix[index] = nactive 
##             # or better use nenglobed?? That would be the same as assuming that flagged pixels have only sky
##             npixout[index] = nnull
##        
##        index += 1
##    
##    radiusflags = 0L
##    
##    goodradii = num.where(radiusflag !=0)[0]
##    if len(goodradii) < 8 : radiusflags = addflag(radiusflags,allflags['BADRADIAL'])
##    if len(goodradii) < 4 : 
##        radiusflags = addflag(radiusflags,allflags['NORADIAL'])
##        larousse = {'radii':None,'cumulflx':None,'intens':None,'npix':None,'npixout':None,'radiusflags':radiusflags}
##    else:
##        radii = radii[goodradii]
##        cumulflx = cumulflx[goodradii]
##        npix = npix[goodradii]
##        npixout = npixout[goodradii]
##        
##        spline = interpol.splrep(npix,cumulflx,s=0,k=3)
##        intens = interpol.splev(npix,spline,der=1)  # TAKE CARE WITH THIS!!!
##        from deriv import deriv
##        intens2 = deriv(npix,cumulflx) # I PREFER THIS!
##        
##        larousse = {'radii':radii,'cumulflx':cumulflx,'intens':intens2,'npix':npix,'npixout':npixout,'radiusflags':radiusflags}
##    
##    self['flags'] = addflag(self['flags'],radiusflags)
##    self['M_RADIAL'] = larousse
##    
##    # t2 = time()
##    # lapsus = t2 - t1
##    # print '%s seconds in finding radial profile for object %i' % (lapsus,self['NUMBER'])
    
def SaveRadial(self,filename):
    """Saves a Radial Profile to a text file."""
    # IMPORT STUFF
    import time as t
    import string
    from pdb import set_trace as stop
    # END IMPORT
    
    t_str = t.ctime(t.time())
    
    radial = self['M_RADIAL']
    # radial = {'radii':Vector,'cumulflx':Vector,'intens':Vector,'npix':Vector,\
    #'npixout':Vector,'radiusflags':Scalar}
    radiusflags = radial['radiusflags']
    header = \
    """
    # Radial profile saved on %s
    # radiusflags = %i
    # radii cumulflx intens npix npixout\n""" % (t_str,radiusflags)
    
    f = open(filename,'w')
    
    f.write(header)
    
    doradial = True
    try : nradial = len(radial['radii'])
    except TypeError: doradial = False
    
    if doradial:
        for i in range(nradial):
            # 'radii','cumulflx','intens','npix','npixout'
            radii = string.rjust('%9.4f' % radial['radii'][i],8)
            cumulflx = string.rjust('%11.5f' % radial['cumulflx'][i],15)
            intens = string.rjust('%11.5f' % radial['intens'][i],15)
            npix = string.rjust('%9i' % radial['npix'][i],10)
            npixout = string.rjust('%9i' % radial['npixout'][i],10)
            line = '%s %s %s %s %s\n' % (radii,cumulflx,intens,npix,npixout)
            f.write(line)
    
    f.close()
    
    self['radial_file'] = filename
    
    return None

def radial_v2(self,MaxBadRatio=0.25,IgnoreBadValues = False,dograph=False):
    """Updates in self the radial profile of source."""
    # IMPORT STUFF
    from ellipse import effective_radius, area_superellip,dist_superellipse
    import numpy as num
    from flags import allflags, addflag
    from scipy import interpolate as interpol
    from momlib import mergemasks
    from pdb import set_trace as stop
    from time import time
    from deriv import deriv
    #from pylab import imshow,show,plot
    # END IMPORT STUFF
 
    t1 = time()
    
    # INPUTS
    stamp = self['STAMP'].copy()
    # image masking
    # mask = self['MASK'] # WRONG! the segmentation map of the object is included in the mask!
    mask = self['MASKOTHER'].copy()
    xcenter = self['X_IMAGE']  ; ycenter = self['Y_IMAGE']
    xcenter = xcenter - self['MXMIN_IMAGE']
    ycenter = ycenter - self['MYMIN_IMAGE']
    ellip = self['ELLIPTICITY'] 
    pa = self['THETA_IMAGE'] - 90. # Astronomical position angle.
    c = self['BOXY']
    sky = self['BACKGROUND']
    # END INPUTS
    
    # initialization of radial distances.
    
    # minimum radius : 0 pix
    # maximum radius : depends on the area given in stamp.
    # logarithmic increase of factor 0.1
    # for each radius > 5, the area (in pixels) englobed should be 
    # computed and compared against theoretical area. 
    # When they differ by >50%, that radius and bigger are dropped.
    
    n = stamp.shape
    q = 1. - ellip
    center = (ycenter,xcenter)
    radial_mask = dist_superellipse(n,center,q=q,pos_ang=pa,c=c)
    
    # Geometric progression of radius
    if self.execpars['doLogRadial'][0] == 1:
        logminradius = 0.
        logmaxradius = num.log10(max(radial_mask.flat))
        logstep = 0.05
        fini_radial = num.arange(logminradius,logmaxradius+logstep,logstep,type='Float32')
        fini_radial = 10.**fini_radial
        ini_radial = num.zeros(shape=(len(fini_radial)+1,),type='Float32')
        ini_radial[1:] = fini_radial.copy()
    else:
        # ALTERNATIVE SHOT: Mixture of linear and geometric progression of radius.
        logminradius = 0.
        logmaxradius = min(num.log10(max(radial_mask.flat)),1)
        logstep = 0.1
        fini_radial = num.arange(logminradius,logmaxradius+logstep,logstep,type='Float32')
        fini_radial = 10.**fini_radial
        if max(radial_mask.flat) > 10. : 
            lini_radial = num.arange(10.,max(radial_mask.flat),type='Float32')
            if len(lini_radial) > 1 : lini_radial = lini_radial[1:]
        ini_radial = num.zeros(shape=(len(fini_radial)+len(lini_radial)+1,),type='Float32')
        ini_radial[1:len(fini_radial)+1] = fini_radial.copy()
        ini_radial[len(fini_radial)+1:] = lini_radial.copy()
    
    radiusflag = num.ones(shape=ini_radial.shape,type='Int8')
    
    cumulflx = num.array(shape=ini_radial.shape,type='Float32')
    radii = num.array(shape=ini_radial.shape,type='Float32')
    npix = num.array(shape=ini_radial.shape,type='Int32')
    npixout = num.array(shape=ini_radial.shape,type='Int32')
    cumulflx[:] = -1.0 ; radii[:] = -1.0
    npix[:] = -1 ; npixout[:] = -1
    
    index = 0
    
    for radius in ini_radial:
        englobed = num.where(radial_mask <= radius)
        nenglobed = len(englobed[0])
        dummie = mask + (radial_mask > radius).astype('Int8') * 2
        active = num.where(dummie == 0)
        boundary = (0 in active[0]) or (stamp.shape[0]-1 in active[0]) \
           or (0 in active[1]) or (stamp.shape[1]-1 in active[1])
        if boundary : 
            radiusflag[index:] = 0
            break
        else:
            nactive = len(active[0])
            null = num.where(dummie == 1) ; nnull = len(null[0])
            #imshow(dummie,origin='lower') ; show()
            #plot(radii,cumulflx,'ro') ; show()
            del dummie
           
            effradius = effective_radius(nenglobed,q,c=0)
            # effradius = radius
            
            cumul = num.sum(stamp[active]) - sky * nactive
            
            if (effradius in radii) or (nnull > MaxBadRatio * nenglobed) or (nenglobed in npix) or (cumul in cumulflx): 
                radiusflag[index] = 0
                # "nenglobed in npix" and "cumul in cumulflx" avoids 
                # duplication of values which make splrep crash
            else : 
                radii[index] = effradius
                cumulflx[index] = cumul
                npix[index] = nenglobed
                npixout[index] = nnull
            index += 1
    
    radiusflags = 0L
    
    goodradii = num.where(radiusflag !=0)[0]
    radii = radii[goodradii]
    cumulflx = cumulflx[goodradii]
    npix = npix[goodradii]
    npixout = npixout[goodradii]
    
    if not IgnoreBadValues:
    
        # CHECK FOR NEGATIVE GROWTH
        
        maxcumul = cumulflx.max()
        maxcumulindex = cumulflx.argmax()
        if maxcumulindex < len(cumulflx)-1:
            postmax = (num.arange(maxcumulindex+1,len(cumulflx)),)
            cumpostmax = cumulflx[postmax]
            below90 = num.where(cumpostmax < maxcumul * 0.9)
            if len(below90[0]) != 0:
               lastgood = postmax[0][below90][0]
               radii = radii[0:lastgood]
               cumulflx = cumulflx[0:lastgood]
               npix = npix[0:lastgood]
               npixout = npixout[0:lastgood]
               print 'Growth curve going negative!'
               radiusflags = addflag(radiusflags,allflags['NEGGROWTH'])
            else : pass
        else : pass
    
        dointens = False
    
        if len(radii) < 8 : radiusflags = addflag(radiusflags,allflags['BADRADIAL'])
        if len(radii) < 4 : 
            radiusflags = addflag(radiusflags,allflags['NORADIAL'])
            larousse = {'radii':None,'cumulflx':None,'intens':None,'npix':None,'npixout':None,'radiusflags':radiusflags}
        else: doIntens = True
    
    if IgnoreBadValues or doIntens:
        npixin = npix - npixout
        incpixin = num.zeros(shape=cumulflx.shape,type='Float32')
        flxinc = num.zeros(shape=cumulflx.shape,type='Float32')
        
        incpixin[0] = 1
        incpixin[1:] = npixin[1:] - npixin[0:-1]
        x = int(num.around(num.array([xcenter]))[0])
        y = int(num.around(num.array([ycenter]))[0])
        flxinc[0] = stamp[y,x] - sky
        flxinc[1:] = cumulflx[1:] - cumulflx[0:-1]
        intens = flxinc / incpixin.astype('Float32')
        nintens = intens.copy()
        nintens[1:-1] = (intens[0:-2] + intens[1:-1])/2.
        
        intens = nintens.copy()
        
        #from pylab import plot,show
        #plot(radii,intens) ; show()
        #plot(radii,cumulflx) ; show()
        
    larousse = {'radii':radii,'cumulflx':cumulflx,'intens':intens,'npix':npix,'npixout':npixout,'radiusflags':radiusflags}
    
    if IgnoreBadValues : 
        radiusflags = addflag(radiusflags,allflags['NONCHECKEDRADIAL'])
    
    self['flags'] = addflag(self['flags'],radiusflags)
    self['M_RADIAL'] = larousse
    
    t2 = time()
    lapsus = t2 - t1
    print '%s seconds in finding radial profile for object %i' % (lapsus,self['NUMBER'])

    if dograph: self.radial_graph()
    
    return None


def radial_v3(self,MaxBadRatio=0.25,ExtrapolBadPix=False,IgnoreBadValues=False,dograph=False):
    """Updates in self the radial profile of source."""
    # This version computes intensity profile as difference in accumulated
    # intensities... I'll try to retrieve the same outputs version 1 does, plus
    # errors in intensities and acccumulated flux.
    # IMPORT STUFF
    from ellipse import effective_radius, area_superellip, dist_superellipse
    import numpy as num
    from flags import allflags, addflag
    from scipy import interpolate as interpol
    from momlib import mergemasks
    from pdb import set_trace as stop
    from time import time
    from deriv import deriv
    # END IMPORT STUFF
    
    # INPUTS
    stamp = self['STAMP'].copy()
    # image masking
    mask = self['MASKOTHER'].copy()
    # geometrical parameters
    xcenter = self['X_IMAGE']  ; ycenter = self['Y_IMAGE']
    xcenter = xcenter - self['MXMIN_IMAGE']
    ycenter = ycenter - self['MYMIN_IMAGE']
    ellip = self['ELLIPTICITY'] 
    pa = self['THETA_IMAGE'] - 90. # Astronomical position angle.
    c = self['BOXY']
    sky = self['BACKGROUND']
    # END INPUTS
    
    # initialization of radial distances.
    
    # minimum radius : 1 pix
    # maximum radius : depends on the area given in stamp.
    # linear increase.
    
    n = stamp.shape
    q = 1. - ellip
    center = (ycenter,xcenter)
    radial_mask = dist_superellipse(n,center,q=q,pos_ang=pa,c=c)
    
    # Progression of radius
    LogRadial = self.execpars['doLogRadial'][0] == 1
    FineRadial = self.execpars['doFineRadial'][0] == 1
    step = self.execpars['RadialStep'][0]
    
    if LogRadial and not FineRadial:
        logminradius = 0.
        logmaxradius = num.log10(max(radial_mask.flat))
        logstep = 0.03
        fini_radial = num.arange(logminradius,logmaxradius+logstep,\
        logstep,type='Float32')
        fini_radial = 10.**fini_radial
        ini_radial = num.zeros(shape=(len(fini_radial)+1,),type='Float32')
        ini_radial[1:] = fini_radial.copy() # ADD r = 0.
        # r = [0.,1.,1.07,1.15,1.23,1.3,...]
    elif not LogRadial and not FineRadial:
        # ALTERNATIVE SHOT: Mixture of linear and geometric progression of 
        # radius.
        logminradius = 0.
        logmaxradius = min(num.log10(max(radial_mask.flat)),1.)
        logstep = 0.03
        fini_radial = num.arange(logminradius,logmaxradius+logstep,\
        logstep,type='Float32')
        fini_radial = 10.**fini_radial
        if max(radial_mask.flat) > 10. : 
            lini_radial = num.arange(10.,max(radial_mask.flat),type='Float32')
            if len(lini_radial) > 1 : lini_radial = lini_radial[1:]
        ini_radial = num.zeros(shape=(len(fini_radial)+len(lini_radial)+1,),\
        type='Float32')
        ini_radial[1:len(fini_radial)+1] = fini_radial.copy()
        ini_radial[len(fini_radial)+1:] = lini_radial.copy()
    elif FineRadial:
        maxradius = max(radial_mask.flat)
        ini_radial = num.arange(float(maxradius))
    
    radiusflag = num.ones(shape=ini_radial.shape,type='Int8')
    
    cumulflx = num.array(shape=ini_radial.shape,type='Float32')
    radii = num.array(shape=ini_radial.shape,type='Float32')
    npix = num.array(shape=ini_radial.shape,type='Int32')
    npixout = num.array(shape=ini_radial.shape,type='Int32')
    cumulflx[:] = -1.0 ; radii[:] = -1.0
    npix[:] = -1 ; npixout[:] = -1
    
    index = 0
    
    for radius in ini_radial:
        englobed = num.where(radial_mask <= radius)
        nenglobed = len(englobed[0])
        dummie = mask + (radial_mask > radius).astype('Int8') * 2
        active = num.where(dummie == 0)
        boundary = (0 in active[0]) or (stamp.shape[0]-1 in active[0]) \
        or (0 in active[1]) or (stamp.shape[1]-1 in active[1])
        if boundary :
            radiusflag[index:] = 0
            break
        else:
            nactive = len(active[0])
            null = num.where(dummie == 1) ; nnull = len(null[0])
            del dummie
            
            effradius = effective_radius(nenglobed,q,c=0)
            # effradius = radius
            
            cumul = num.sum(stamp[active]) - sky * nactive
            try : corrfactor = float(nenglobed)/float(nactive)
            except ZeroDivisionError : corrfactor = 1.
            if ExtrapolBadPix : cumul = cumul * corrfactor
            #print cumul/nenglobed,1./corrfactor
            if (effradius in radii) or (nnull > MaxBadRatio * nenglobed) or \
            (nenglobed in npix) or (cumul in cumulflx): 
                radiusflag[index] = 0
                # "nenglobed in npix" and "cumul in cumulflx" avoids 
                # duplication of values which make splrep crash
            else : 
                radii[index] = effradius
                cumulflx[index] = cumul
                npix[index] = nenglobed
                npixout[index] = nnull
            index += 1
    
    radiusflags = 0L
    
    goodradii = num.where(radiusflag != 0)[0]
    radii = radii[goodradii]
    cumulflx = cumulflx[goodradii]
    npix = npix[goodradii]
    npixout = npixout[goodradii]
    
    if not IgnoreBadValues:
        
        # CHECK FOR NEGATIVE GROWTH
        
        maxcumul = cumulflx.max()
        maxcumulindex = cumulflx.argmax()
        if maxcumulindex < len(cumulflx)-1:
            postmax = (num.arange(maxcumulindex+1,len(cumulflx)),)
            cumpostmax = cumulflx[postmax]
            below90 = num.where(cumpostmax < maxcumul * 0.9)
            if len(below90[0]) != 0:
               lastgood = postmax[0][below90][0]
               radii = radii[0:lastgood]
               cumulflx = cumulflx[0:lastgood]
               npix = npix[0:lastgood]
               npixout = npixout[0:lastgood]
               print 'Growth curve going negative!'
               radiusflags = addflag(radiusflags,allflags['NEGGROWTH'])
            else : pass
        else : pass
        
        if len(radii) < 8 : radiusflags = addflag(radiusflags,allflags['BADRADIAL'])
        if len(radii) < 4 : 
            radiusflags = addflag(radiusflags,allflags['NORADIAL'])
            radii = None ; cumulflx = None ; intens2 = None ; 
            npix = None ; npixout = None ; 
        else:
            # spline = interpol.splrep(npix,cumulflx,s=0,k=3)
            intens2 = deriv(npix,cumulflx)
    else:
        radiusflags = addflag(radiusflags,allflags['NONCHECKEDRADIAL'])
        if len(npix) > 3:
            intens2 = deriv(npix,cumulflx)
        else : 
            radiusflags = addflag(radiusflags,allflags['NORADIAL'])
            radii=None ; cumulflx=None
            intens2=None ; npix=None ; npixout=None
    
    larousse = {'radii':radii,'cumulflx':cumulflx,'intens':intens2,'npix':npix,'npixout':npixout,'radiusflags':radiusflags}
    
    self['flags'] = addflag(self['flags'],radiusflags)
    self['M_RADIAL'] = larousse
    
    if dograph: self.radial_graph()
    
    return None
    
def back_radial(intens,radii,skyradius):
    """Estimates Background from the intensity curve."""
    # IMPORT STUFF
    import numpy as num
    # END IMPORT
    
    skyreg = num.where(radii > skyradius)
    if len(skyreg) > 0:
        skyval = num.median(intens[skyreg])
    else : skyval = None
    
    return skyval
    
    
