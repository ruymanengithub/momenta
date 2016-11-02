#! /usr/bin/env python

# IMPORT STUFF
from ellipse import effective_radius, area_superellip, dist_superellipse
import numpy as num
from flags import allflags, addflag
from pdb import set_trace as stop
from time import time
# END IMPORT STUFF

def radial_v1(self,MaxBadRatio=0.25,ExtrapolBadPix=False,\
    IgnoreBadValues=False,dograph=False):
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
    # If less than 25% [optional] of the pixels in an apperture or more are
    # flagged out, the apperture is not valid.
    # The initial list of appertures is fixed. The final list of appertures 
    # depends on where is the center of the appertures and the flagged pixels.
    
    # DIFFERENCES WITH RESPECT TO VERSION 0:
    # npix is now the number of pixels inside the elliptical aperture, including  
    # the flagged pixels.
    
    # IMPORT STUFF
    from scipy import interpolate as interpol
    from momlib import mergemasks
    from deriv import deriv
    # END IMPORT STUFF
    
    # INPUTS
    stamp = self['STAMP'].copy()
    # image masking
    mask = self['MASKOTHER'].copy()
    objectmask = self['MASK'].copy()
    xcenter = self['X_IMAGE']  ; ycenter = self['Y_IMAGE']
    xcenter = xcenter - self['MXMIN_IMAGE']
    ycenter = ycenter - self['MYMIN_IMAGE']
    ellip = self['ELLIPTICITY'] 
    pa = self['THETA_IMAGE'] - 90. # Astronomical position angle.
    c = self['BOXY']
    sky = self['BACKGROUND']
    doRadSky = self.execpars['doRadSky'][0] == 1
    # END INPUTS
    
    # initialization of radial distances.
    
    # minimum radius : 1 pix
    # maximum radius : depends on the area given in stamp.
    # logarithmic increase of factor 0.1 [optional]
    
    n = stamp.shape
    q = 1. - ellip
    center = (ycenter,xcenter)
    radial_mask = dist_superellipse(n,center,q=q,pos_ang=pa,c=c)
    
    # Progression of radius
    LogRadial = self.execpars['doLogRadial'][0] == 1
    try: step = self.execpars['RadialStep'][0]
    except KeyError:
        if LogRadial : step = 0.03 # delta_r = 3%
        else : step = 1. # delta_r = 1 pix
    
    if LogRadial:
        logminradius = 0.
        logmaxradius = num.log10(max(radial_mask.flat))
        logstep = step
        fini_radial = num.arange(logminradius,logmaxradius+logstep,\
        logstep,dtype='float32')
        fini_radial = 10.**fini_radial
        ini_radial = num.zeros(shape=(len(fini_radial)+1,),dtype='float32')
        ini_radial[1:] = fini_radial.copy()
    elif not LogRadial:
        linstep = step
        maxradius = max(radial_mask.flat)
        ini_radial = num.arange(0.,float(maxradius)+linstep,linstep,\
        dtype='float32')
        # r = [0.,1.,2.,3,...]
    
    radiusflag = num.ones(shape=ini_radial.shape,dtype='int8')
    cumulflx = num.array(shape=ini_radial.shape,dtype='float32')
    radii = num.array(shape=ini_radial.shape,dtype='float32')
    npix = num.array(shape=ini_radial.shape,dtype='int32')
    npixout = num.array(shape=ini_radial.shape,dtype='Int32')
    cumulflx[:] = -1.0 ; radii[:] = -1.0
    npix[:] = -1 ; npixout[:] = -1
    
    index = 0
    
    for radius in ini_radial:
        englobed = num.where(radial_mask <= radius)
        nenglobed = len(englobed[0])
        dummie = mask + (radial_mask > radius).astype('int8') * 2
        active = num.where(dummie == 0)
        boundary = (0 in englobed[0]) or (stamp.shape[0]-1 in englobed[0]) \
        or (0 in englobed[1]) or (stamp.shape[1]-1 in englobed[1])
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
    
    
    if not IgnoreBadValues and not doRadSky:
        
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
            
            intens = deriv(npix,cumulflx)
            
            # Background estimation from intensity profile, substracion 
            # and update of BACKGROUND.
            if doRadSky:    
                areaobject = float(len(num.where(objectmask==0)[0]))
                skyradius_old = (areaobject / (num.pi * (1.-ellip)))**0.5
                
                ozero = num.where(objectmask==0)
                skyradius = float(max(max(ozero[0])-min(ozero[0]),\
                max(ozero[1])-min(ozero[1])))/2.
                skyradius *= 1.5 # skyradius a 50% bigger than object radius, +/-
                
                back2 = back_radial(intens,radii,skyradius)
                
                if not num.isnan(back2):
                    backtotal = sky+back2
                    self['BACKGROUND'] = backtotal
                    cumulflx = cumulflx - back2 * (npix-npixout)
                    intens = intens - back2
                else: radiusflags = addflag(radiusflags,allflags['SHORTSKYRAD'])
     
        else : 
            radiusflags = addflag(radiusflags,allflags['NORADIAL'])
            radii=None ; cumulflx=None
            intens=None ; npix=None ; npixout=None
    
    if len(radii) < 8 : radiusflags = addflag(radiusflags,allflags['BADRADIAL'])
    
    larousse = {'radii':radii,'cumulflx':cumulflx,'intens':intens,
    'npix':npix,'npixout':npixout,'radiusflags':radiusflags}
    
    self['flags'] = addflag(self['flags'],radiusflags)
    self['M_RADIAL'] = larousse
    
    if dograph: self.radial_graph()
    
    return None
 

def radial_v3(self,MaxBadRatio=0.25,ExtrapolBadPix=False,\
    IgnoreBadValues=True,dograph=False):
    """Updates in self the radial profile of source.
    This version computes intensity profile as difference in accumulated
    fluxes. It retrieves the same outputs version 1 does, plus errors 
    in intensities and acccumulated flux."""
    # IMPORT STUFF
    # END IMPORT STUFF
    
    # INPUTS
    stamp = self['STAMP'].copy()
    # image masking
    mask = self['MASKOTHER'].copy()
    objectmask = self['MASK'].copy()
    # geometrical parameters
    xcenter = self['X_IMAGE']  ; ycenter = self['Y_IMAGE']
    xcenter = xcenter - self['MXMIN_IMAGE']
    ycenter = ycenter - self['MYMIN_IMAGE']
    ellip = self['ELLIPTICITY']
    pa = self['THETA_IMAGE'] - 90. # Astronomical position angle
    c = self['BOXY']
    sky = self['BACKGROUND']
    skysigma = self.execpars['sigma_sky'][0]
    doRadSky = self.execpars['doRadSky'][0] == 1
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
    try: step = self.execpars['RadialStep'][0]
    except KeyError:
        if LogRadial : step = 0.03 # delta_r = 3%
        else : step = 1. # delta_r = 1 pix

    if LogRadial:
        logminradius = 0.
        logmaxradius = num.log10(max(radial_mask.flat))
        fini_radial = num.arange(logminradius,logmaxradius+step,\
        step,dtype='float32')
        fini_radial = 10.**fini_radial
        ini_radial = num.zeros(shape=(len(fini_radial)+1,),dtype='float32')
        ini_radial[1:] = fini_radial.copy() # ADD r = 0.
        # r = [0.,1.,1.07,1.15,1.23,1.3,...]
    else:
        maxradius = max(radial_mask.flat)
        ini_radial = num.arange(0.,float(maxradius)+step,step,\
        dtype='float32')
        # r = [0.,1.,2.,3,...] for step = 1
    
    stamp -= sky # SKY SUBSTRACTION
    
    cumulflx = [] ; ecumulflx = []
    intens = [] ; eintens = []
    npix = [] ; npixout = []
    radii = []
    
    for i in range(len(ini_radial)):
        
        # Initial radii : previous, actual and posterior
        radius = ini_radial[i]
        
        prevrad = ini_radial[max(0,i-1)]
        try: postrad = ini_radial[i+1]
        except IndexError : postrad = radius + (radius-prevrad)
        innerann = radius - max(1.,float(num.around(radius-prevrad)/2.))
        outerann = radius + max(1.,float(num.around(postrad-radius)/2.))
        if innerann<0.: innerann = 0.#; outerann = 0.
        
        englobed = radial_mask < radius
        englobedx = num.where(englobed)
        nenglobed = len(englobedx[0])
        cumactive = num.where(englobed & (mask == 0))
        ncumactive = len(cumactive[0])
        # next englobed area and previous : intensity area
        # if step > 1 pix, then it takes an anulus which goes half-way in both
        # directions.
        annular = (radial_mask >= innerann) &\
        (radial_mask<outerann)
        if radius == 0.:
            annular = radial_mask == radial_mask.min()
        
        annularx = num.where(annular)
        
        nannular = len(annularx[0])
        annactive = num.where(annular & (mask==0))
        nannactive = len(annactive[0])
        
        boundary = (0 in annularx[0]) or (stamp.shape[0]-1 in annularx[0]) \
        or (0 in annularx[1]) or (stamp.shape[1]-1 in annularx[1])
        	
        try: 
            BadCRatio = float(nenglobed-ncumactive)/float(nenglobed)
            isGoodCRatio = BadCRatio <= MaxBadRatio
        except ZeroDivisionError : isGoodCRatio = False
	
        try:
            BadARatio = float(nannular-nannactive)/float(nannular)
            isGoodARatio = BadARatio <= MaxBadRatio
        except ZeroDivisionError : isGoodARatio = False
	
        if not boundary:
	    radii.append(radius)
        else: break
        
        if ncumactive > 0 and isGoodCRatio:
            
            npix.append(float(nenglobed))
            npixout.append(float(nenglobed-ncumactive))
            
            cum = num.sum(stamp[cumactive])
            #print time() - time0
            cumcorr = float(nenglobed)/float(ncumactive)
            if ExtrapolBadPix : cum = cum * cumcorr
            ecum = skysigma / num.sqrt(ncumactive)
            
            cumulflx.append(cum)
            ecumulflx.append(ecum)	    
        elif nenglobed == 0:
            cumulflx.append(0.)
            ecumulflx.append(0.)
            npix.append(0.)
            npixout.append(0.)
        else:
            cumulflx.append(num.nan)
            ecumulflx.append(num.nan)
            npix.append(0.)
            npixout.append(0.)
        	
        if nannactive > 0 and isGoodARatio:
	
            sinten = float(num.median(stamp[annactive]))          
            #if radius == 0: stop()
            
            if nannular > 1:
                esinten = float(num.std(stamp[annactive]) / num.sqrt(nannactive))
            else:
                esinten = skysigma
            
            intens.append(sinten)
            eintens.append(esinten)
            
        else:
            
            intens.append(num.nan)
            eintens.append(num.nan)
    
    cumulflx = num.array(cumulflx,dtype='float64')
    ecumulflx = num.array(ecumulflx,dtype='float64')
    intens = num.array(intens,dtype='float64') 
    eintens = num.array(eintens,dtype='float64')
    npix = num.array(npix) ; npixout = num.array(npixout)
    radii = num.array(radii)
    
    notnan = num.where((num.isnan(cumulflx) == False) & \
        (num.isnan(intens) == False))
    
    cumulflx = cumulflx[notnan] ; ecumulflx = ecumulflx[notnan]
    intens = intens[notnan] ; eintens = eintens[notnan]
    npix = npix[notnan] ; npixout = npixout[notnan]
    radii = radii[notnan]
    
    radiusflags = 0L
    
    if not IgnoreBadValues and not doRadSky:
        
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
               ecumulflx = ecumulflx[0:lastgood]
               intens = intens[0:lastgood] ; eintens = eintens[0:lastgood]
               npix = npix[0:lastgood]
               npixout = npixout[0:lastgood]
               print 'Growth curve going negative!'
               radiusflags = addflag(radiusflags,allflags['NEGGROWTH'])
            else : pass
        else : pass
        
        if len(radii) < 4 : 	    
            radiusflags = addflag(radiusflags,allflags['NORADIAL'])
            radii = None ; cumulflx = None ; ecumulflx = None
            intens = None ; eintens = None
            npix = None ; npixout = None ; 
        else:
            pass
    else:
        
        radiusflags = addflag(radiusflags,allflags['NONCHECKEDRADIAL'])
        
        if len(npix) > 3:
            # Background estimation from intensity profile, substracion and 
            # update of BACKGROUND.
            if doRadSky:    
                areaobject = float(len(num.where(objectmask==0)[0]))
                skyradius_old = (areaobject / (num.pi * (1.-ellip)))**0.5
                ozero = num.where(objectmask==0)
                skyradius = float(max(max(ozero[0])-min(ozero[0]),\
                max(ozero[1])-min(ozero[1])))/2.
                skyradius *= 1.5 # skyradius a 50% bigger than object radius, +/-
                
                back2 = back_radial(intens,radii,skyradius)
                
                if not num.isnan(back2):
                    backtotal = sky+back2
                    self['BACKGROUND'] = backtotal
                    cumulflx = cumulflx - back2 * (npix-npixout)
                    intens = intens - back2
                else: 
                    radiusflags = addflag(radiusflags,allflags['SHORTSKYRAD'])
        else : 
            radiusflags = addflag(radiusflags,allflags['NORADIAL'])
            radii=None ; cumulflx=None ; ecumulflx = None
            intens=None ; eintens = None
            npix=None ; npixout=None
    
    try: 
        if len(radii) < 8 : 	    
	    radiusflags = addflag(radiusflags,allflags['BADRADIAL'])
    except : 
        pass
        
    
    larousse = {'radii':radii,'cumulflx':cumulflx,'ecumulflx':ecumulflx,\
    'intens':intens,'eintens':eintens,'npix':npix,'npixout':npixout,\
    'radiusflags':radiusflags}
    self['flags'] = addflag(self['flags'],radiusflags)
    self['M_RADIAL'] = larousse
    
    
    if dograph: self.radial_graph()
    
    return None
    
def back_radial(intens,radii,skyradius):
    """Estimates Background from the intensity curve.
    BEWARE: if profile is shorter than skyradius, it retrieves num.nan."""
    # IMPORT STUFF
    # END IMPORT
    
    skyreg = num.where(radii > skyradius)
    if len(skyreg[0]) >= 3:
        skyval = float(num.median(intens[skyreg]))
    else : 
        #skyval = float(num.median(intens[-3:]))
        skyval = num.nan
    
    return skyval
    
def SaveRadial(self,filename):
    """Saves a Radial Profile to a text file."""
    # IMPORT STUFF
    import time as t
    import string
    # END IMPORT
    
    t_str = t.ctime(t.time())
    
    radial = self['M_RADIAL']
    # radial = {'radii':Vector,'cumulflx':Vector,'intens':Vector,'npix':Vector,\
    #'npixout':Vector,'radiusflags':Scalar}
    radiusflags = radial['radiusflags']
    header = \
    """# Radial profile saved on %s
# radiusflags = %i
# radii cumulflx ecumul intens eintens npix npixout\n""" % (t_str,radiusflags)
    
    f = open(filename,'w')
    
    f.write(header)
    
    doradial = True
    try : nradial = len(radial['radii'])
    except TypeError: doradial = False
    
    if doradial:
        for i in range(nradial):
            # 'radii','cumulflx','intens','npix','npixout'
            radii = string.rjust('%.3e' % radial['radii'][i],8)
            cumulflx = string.rjust('%.3e' % radial['cumulflx'][i],15)
            try : ecumulflx = string.rjust('%.3e' % radial['ecumulflx'][i],15)
            except KeyError:
                ecumulflx = string.rjust('%.3e' % -1.0,15)
            intens = string.rjust('%.3e' % radial['intens'][i],15)
            try: eintens = string.rjust('%.3e' % radial['eintens'][i],15)
            except KeyError:
                eintens = string.rjust('%.3e' % -1.0,15)
            npix = string.rjust('%i' % radial['npix'][i],10)
            npixout = string.rjust('%i' % radial['npixout'][i],10)
            line = '%s %s %s %s %s %s %s\n' % (radii,cumulflx,ecumulflx,\
            intens,eintens,npix,npixout)
            f.write(line)
    
    f.close()
    
    self['radial_file'] = filename
    
    return None
    
def RadialModel(radii,intens,dims,center,q,pa,c=0.):
    """Retrieves an image with an axisymmetric 2D model for a radial
    distribution of intensity.
    
    Inputs: radii, intens,dims,center,q,pa,c"""
    # IMPORT STUFF
    from scipy import interpolate as interp
    # END IMPORT
    
    n = tuple(dims[-1::-1]) # y,x!
    ecenter = tuple(center[-1::-1]) # y,x!
    radial_mask = dist_superellipse(n,ecenter,q=q,pos_ang=pa,c=c)
    arg = num.where(radial_mask<=radii.max())
    ipoints = radial_mask[arg]
    Model = radial_mask * 0.
    
    Spline = interp.splrep(radii,intens,k=1,s=0)
    ESpline = interp.splev(ipoints,Spline)
    Model[arg] = ESpline.copy()
    
    
    return Model
