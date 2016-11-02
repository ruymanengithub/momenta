#! /usr/bin/env python

"""Several tasks used by MassClumps.py"""

# IMPORT STUFF
import numpy as num
from pdb import set_trace as stop
from copy import copy
# END IMPORT

def MassProf(self):
    """Returns a Mass Profile, from Intensity profiles."""
    # IMPORT STUFF
    # END IMPORT

def MassProfMap(self):
    """Returns a Mass Profile, given a Mass Map."""
    # IMPORT STUFF
    # END IMPORT
    
    dummy = copy(self)
    
    dummy.execpars['sigma_sky'] = [1.]
    dummy.self['BACKGROUND'] = 0.
    dummy.self['MASKOTHER'] = self[''].copy
    
    dummy.radial_v3(MaxBadRatio=0.50,ExtrapolBadPix=False,\
    IgnoreBadValues=True,dograph=False)
    
    self['MassProf'] = copy(dummy['M_RADIAL'])
    
    return None
    
def MassAz(self):
    """Returns an azimuthally averaged model of a galaxy's Mass Distribution.
    It uses a Mass Profile and geometric parameters (center,ellipticity,thetha),
    and a mask."""
    # IMPORT STUFF
    from Moments.ellipse import dist_superellipse
    from Moments.radial import RadialModel
    # END IMPORT
    
    # INPUTS
    massmap = self[''].copy()
    massmask = self[''].copy()
    massprofile = self[''].copy()
    xcenter = self['X_IMAGE']  ; ycenter = self['Y_IMAGE']
    xcenter = xcenter - self['MXMIN_IMAGE']
    ycenter = ycenter - self['MYMIN_IMAGE']
    ellip = self['ELLIPTICITY'] 
    pa = self['THETA_IMAGE'] - 90. # Astronomical position angle.
    c = self['BOXY']
    petrofactor = self.execpars['petrofactor'][0]
    maxradius = petrofactor * self['R_PETRO']
    # END INPUTS
    
    radii = massprofile['radii']
    imass = massprofile['intensity']
    
    n = massmodel.shape
    dims = n[-1::-1]
    q = 1. - ellip
    icenter = (xcenter,ycenter)
    
    MassAz = RadialModel(radii,imass,dims,icenter,q,pa,c=0.)
    
##    # ALTERNATIVE:
##    MassAz = massmap * 0.
##    center = (ycenter,xcenter)
##    radialmap = dist_superellipse(n,center,q=q,pos_ang=pa,c=c)
##    
##    for i in range(0,len(radii)):
##        if i == 0:
##           selected = num.where(radialmap<radii[i])
##        if i == (len(radii)-1):
##           selected = num.where(radialmap>=radii[i])
##        else:
##           radiusbef = (radii[i-1]+radii[i])/2.
##           radiusaft = (radii[i]+radii[i+1])/2.
##           selected = num.where((radialmap>=radiusbef) & \
##           (radialmap<radiusaft))
##        
##        MassAz[selected] = imass[i]
    
    MassAz[num.where(massmask)] = 0.
    external = num.where(radialmap>=maxradius)
    MassAz[external] = 0.
    negative = num.where(MassAz<0.)
    MassAz[negative] = 0.
    
    self['MassAz'] = MassAz
    
    return None
    
def ClumpsMap(self):
    """Returns an image with Map clumps, as residuals of a Mass map and
    and an azimuthally averaged model of that distribution."""
    # IMPORT STUFF
    # END IMPORT
    # INPUTS
    massmap = self[''].copy()
    massmask = self[''].copy()
    massaz = self[''].copy()
    # END INPUTS
    
    clumpsmap = massmap - massaz
    negative = num.where(clumpsmap<0.)
    clumpsmap[negative] = 0.
    masked = num.where(massmask)
    clumpsmap[masked] = 0.
    
    self['ClumpsMap'] = clumpsmap.copy()
    return None
    
def MassClumps_Funct(self):
    """Analyzes the mass clumps."""
    # IMPORT STUFF
    
    # END IMPORT
    
    # Segmentation of clumps.
    
    # Mass statistics.
    
    # Sizes statistics.
    
    # Radial distribution.
