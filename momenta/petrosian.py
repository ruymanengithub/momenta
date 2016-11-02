#! /usr/bin/env python

"""Tasks related to Petrosian radius"""

# IMPORT STUFF
import numpy as num
import sys
from flags import allflags, addflag, isflagon
from pdb import set_trace as stop 
from time import time
# END IMPORT

def petrosian(self):
    """Gets petrosian radius and related quantities."""
    # Key points:
    # start from growth curve and intensity curves from radial.py
    # build equation to solve for petrosian radius and interpolate it.
    # Compute petrosian by finding roots of preceding function with
    # a root finder for B-splines (such as scipy.interpolate.sproof).
    # IMPORT STUFF
    import scipy.interpolate as interp
    from pylab import plot,show
    # END IMPORT

    rprof = self['M_RADIAL'] # radial profile. list
    cumulflx = rprof['cumulflx'] # Accumulated flux.
    intflx = rprof['intens'] # Intensity.
    radial = rprof['radii'] # Radius.
    npix = rprof['npix'] # Number of pixels.
    eta_petro = self.execpars['eta_petro'][0]
    
    try: usePetroFirst = self.execpars['usePetroFirst'][0] == 1
    except KeyError: 
        stop('NEW: WARNING AT petrosian.petrosian!')
        usePetroFirst = False
    
    try: minMuPetro = self.execpars['minMuPetro'][0]
    except KeyError: 
        stop('NEW: WARNING AT petrosian.petrosian!')
        minMuPetro = -1.0
    
    if minMuPetro > 0.:
        try:
            magzero = self.execpars['magzero'][0]
            scale = self.execpars['scale'][0]
            minIPetro = (scale**2.)*10.**((magzero-minMuPetro)/2.5)
        except:
            print 'Some Requested Parameter not given... \nStop at petrosian.petrosian'
            stop()
    else : minIPetro = -99.0
    
    if isflagon(self['flags'],allflags['NORADIAL']):
        self['R_PETRO'] = -99.0
        self['I_PETRO'] = -99.0
        self['F_PETRO'] = -99.0
        self['flags'] = addflag(self['flags'],allflags['NOPETRO'])
        print 'No radial profile, no petrosian radius...'
        return None
    
    if minIPetro > 0.: 
        belowIPetro = num.where(intflx < minIPetro)
        nozero = num.where(radial>0.)
        if len(nozero[0])>0 and len(belowIPetro[0])>0:
            ixminBelow = min(belowIPetro[0])
            nozero = (num.array([ix for ix in nozero[0] if ix <= ixminBelow]),)
    else: 
        nozero = num.where((radial>0.) & (intflx > 0.))
    
    cumulflx = cumulflx[nozero] ; intflx = intflx[nozero] 
    radial = radial[nozero] ; npix = npix[nozero]
    
    average_cumul = cumulflx / npix
    tosolve_f = intflx / average_cumul - eta_petro
    
    if len(tosolve_f)>3:
        tosolve_spl = interp.splrep(radial,tosolve_f,s=0,k=1)
        f_spl = interp.splrep(radial,cumulflx,s=0,k=1)
        i_spl = interp.splrep(radial,intflx,s=0,k=1)
        try: roots = interp.sproot(tosolve_spl)
        except: roots = num.array([])
    else:
        roots = num.array([])
    
    petro_flags = 0L
    r_petro = -99.0 ; i_petro = -99.0 ; f_petro = -99.0
    
    if len(roots) == 0:  
        petro_flags = addflag(petro_flags,allflags['NOPETRO']) # no petro radius
        print 'No roots for Petrosian Radius'
        try: print 'band : %s' % self.execpars['band']
        except : pass
##        show._needmain=False
##        #plot(radial,tosolve_f,'r-')
##        posp = num.where((intflx>0.) & (average_cumul>0.))
##        plot(radial[posp],intflx[posp],'r-')
##        plot(radial[posp],average_cumul[posp],'b-')
##        plot(radial[posp],\
##        cumulflx[posp]*max(average_cumul[posp])/max(cumulflx[posp]),'g-')
##        show()
##        stop()
    else :
        if usePetroFirst: r_petro = float(min(roots))
        else: r_petro = float(max(roots))
        i_petro = float(interp.splev(r_petro,i_spl))
        f_petro = float(interp.splev(r_petro,f_spl))
        if i_petro < 0. : petro_flags = addflag(petro_flags,\
        allflags['NEGPETRO']) # negative intensity at petro_radius.
    if len(roots) > 1:
        petro_flags = addflag(petro_flags,allflags['MANYPETRO'])
        # more than one petro radius.
    
    self['R_PETRO'] = r_petro
    self['I_PETRO'] = i_petro
    self['F_PETRO'] = f_petro
    self['flags'] = addflag(self['flags'],petro_flags)
    
    print 'PetroRadius = %7.2f\n' % self['R_PETRO']
    
    return None



def petrosian2(self):
    """Gets petrosian radius and related quantities.
    Now with a slightly different algorithm."""
    # Key points:
    # start from growth curve and intensity curves from radial.py
    # Build function to solve (I/<F> - eta = 0)
    # Finds zeros of function
    # the rest is the same as "petrosian"...
    # IMPORT STUFF
    from pylab import plot,show
    import scipy.interpolate as interp
    from momlib import findzeros
    # END IMPORT

    rprof = self['M_RADIAL'] # radial profile. list
    cumulflx = rprof['cumulflx'] # Accumulated flux.
    intflx = rprof['intens'] # Intensity.
    radial = rprof['radii'] # Radius.
    npix = rprof['npix'] # Number of pixels.
    eta_petro = self.execpars['eta_petro'][0]
    
    try: usePetroFirst = self.execpars['usePetroFirst'][0] == 1
    except KeyError: 
        stop('NEW: WARNING AT petrosian.petrosian!')
        usePetroFirst = False
    
    try: minMuPetro = self.execpars['minMuPetro'][0]
    except KeyError: 
        stop('NEW: WARNING AT petrosian.petrosian!')
        minMuPetro = -1.0
    
    if minMuPetro > 0.:
        try:
            magzero = self.execpars['magzero'][0]
            scale = self.execpars['scale'][0]
            minIPetro = (scale**2.)*10.**((magzero-minMuPetro)/2.5)
        except:
            print 'Some Requested Parameter not given... \nStop at petrosian.petrosian'
            stop()
    else : minIPetro = -99.0
    
    if isflagon(self['flags'],allflags['NORADIAL']):
        self['R_PETRO'] = -99.0
        self['I_PETRO'] = -99.0
        self['F_PETRO'] = -99.0
        self['flags'] = addflag(self['flags'],allflags['NOPETRO'])
        print 'No radial profile, no petrosian radius...'
        return None
    
    if minIPetro > 0.: 
        belowIPetro = num.where(intflx < minIPetro)
        nozero = num.where(radial>0.)
        if len(nozero[0])>0 and len(belowIPetro[0])>0:
            ixminBelow = min(belowIPetro[0])
            nozero = (num.array([ix for ix in nozero[0] if ix <= ixminBelow]),)
    else: 
        nozero = num.where((radial>0.) & (intflx > 0.))
    
    cumulflx = cumulflx[nozero] ; intflx = intflx[nozero] 
    radial = radial[nozero] ; npix = npix[nozero]
    
    average_cumul = cumulflx / npix
    tosolve_f = intflx / average_cumul - eta_petro
    
    if len(tosolve_f)>3:
        
        f_i1d = interp.interp1d(radial,cumulflx)
        i_i1d = interp.interp1d(radial,intflx)
        
	roots = findzeros(radial,tosolve_f)
    
    else:
    
        roots = num.array([])
    
    petro_flags = 0L
    r_petro = -99.0 ; i_petro = -99.0 ; f_petro = -99.0
    
    if len(roots) == 0:  
        petro_flags = addflag(petro_flags,allflags['NOPETRO']) # no petro radius
        print 'No roots for Petrosian Radius'
        try: print 'band : %s' % self.execpars['band']
        except : pass
##        show._needmain=False
##        #plot(radial,tosolve_f,'r-')
##        posp = num.where((intflx>0.) & (average_cumul>0.))
##        plot(radial[posp],intflx[posp],'r-')
##        plot(radial[posp],average_cumul[posp],'b-')
##        plot(radial[posp],\
##        cumulflx[posp]*max(average_cumul[posp])/max(cumulflx[posp]),'g-')
##        show()
##        stop()
    else :
        if usePetroFirst: r_petro = float(min(roots))
        else: r_petro = float(max(roots))
        i_petro = float(i_i1d(r_petro))
        f_petro = float(f_i1d(r_petro))
        if i_petro < 0. : petro_flags = addflag(petro_flags,\
        allflags['NEGPETRO']) # negative intensity at petro_radius.
    if len(roots) > 1:
        petro_flags = addflag(petro_flags,allflags['MANYPETRO'])
        # more than one petro radius.
    
    self['R_PETRO'] = r_petro
    self['I_PETRO'] = i_petro
    self['F_PETRO'] = f_petro
    self['flags'] = addflag(self['flags'],petro_flags)
    
    print 'PetroRadius = %7.2f\n' % self['R_PETRO']
    
    return None

    
def petromsk(self):
    """Mask of pixels contained in petrosian mask. 1 is 
    flagged out"""
    # IMPORT STUFF
    from ellipse import dist_superellipse
    # END IMPORT
    
    if isflagon(self['flags'],allflags['NOPETRO']):
        self['M_PETROMSK'] = None
        return None
    
    # INPUTS
    stamp = self['STAMP'].copy()
    xcenter = self['X_IMAGE'] - self['MXMIN_IMAGE']
    ycenter = self['Y_IMAGE'] - self['MYMIN_IMAGE']
    ellip = self['ELLIPTICITY']
    pa = self['THETA_IMAGE'] - 90. # Astronomical position angle
    c = 0
    if self.execpars['useExPetro'][0] == 1:
        r_petro = self['EX_R_PETRO']
    else: 
        r_petro = self['R_PETRO']
    petrofactor = self.execpars['petrofactor'][0]
    # END INPUTS
    n = (len(stamp[:,0]),len(stamp[0,:]))
    center = (ycenter,xcenter)
        
    if self.execpars['docircularpetro'][0] == 1 : q = 1
    else: q = 1 - ellip
        
    petro_mask = dist_superellipse(n,center,q=q,pos_ang=pa,c=c)
    radius = r_petro * petrofactor
   
    petro_mask[num.where(petro_mask <= radius)] = -1
    petro_mask[num.where(petro_mask != -1)] = 1
    petro_mask[num.where(petro_mask == -1)] = 0
    
    active = num.where(petro_mask == 0)
    if 0 in active[0] or 0 in active[1] or len(active[0])-1 in active[0] \
    or len(active[1]-1) in active[1]:
        print 'Petrosian mask truncated by window'
        self['flags'] = addflag(self['flags'],allflags['PETROMSKTRUNC'])
    
    self['M_PETROMSK'] = petro_mask.astype('Int8')
    
    return None
