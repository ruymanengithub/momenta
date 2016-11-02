#! /usr/bin/env python

from pdb import set_trace as stop
import numpy as num

def Excentricity(self,dograph=False):
    """Computes Excentricity parameter."""
    # IMPORT STUFF
    from Excentricity import Ex_GetMinCirc_Aper
    from algorithms import moments
    # END IMPORT
    
    # INPUTS
    sky = self['BACKGROUND']
    image = self['STAMP'].copy() - sky
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    maskother = self['MASKOTHER'].copy()
    sigma_sky = self.execpars['sigma_sky'][0]
    # END INPUTS
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')
    
    mincirc,ccenter,cradius = Ex_GetMinCirc_Aper(image,mask)
    
    active = num.where((maskother == 0) & (mincirc == 0))
    
    xd,yd = moments(image,x0=0.,y0=0.,order=1,active=active)
    
    distance = ((ccenter[1]-xd)**2. + (ccenter[0]-yd)**2.)**0.5
    
    E = 100.*distance / cradius
    
    print 'Excentricity = %.2f o/o' % E
    
    self['M_E'] = E
    
    if dograph:
        ccenterp = (ccenter[1],ccenter[0])
        icenter = (xd,yd)
        self.Excentricity_graph(ccenterp,cradius,icenter)
    
    return None
    
def Ex_GetMinCirc_Aper(image,mask=-1):
    """Finds the minimum area circle which encloses a nonzero distribution
    of pixels"""
    # IMPORT STUFF
    from Moments.ellipse import dist_superellipse
    from Excentricity import Ex_DistMatrix
    # END IMPORT
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')
    active = num.where(mask == 0)
    
    minx = min(active[1]) ; maxx = max(active[1])
    miny = min(active[0]) ; maxy = max(active[0])
    
    extremax = []
    extremay = []
    for i in range(len(active[0])):
        if active[0][i] in (miny,maxy) or active[1][i] in (minx,maxx):
            extremax.append(active[1][i])
            extremay.append(active[0][i])
    extrema = (num.array(extremay),num.array(extremax))
    
    distances  = Ex_DistMatrix(extrema[1],extrema[0])
    
    maxd = distances.max()
    maxdp = num.where(distances == maxd)
    
    indx0 = maxdp[0][0] ; indx1 = maxdp[1][0]
    x0,y0 = extrema[1][indx0],extrema[0][indx0]
    x1,y1 = extrema[1][indx1],extrema[0][indx1]
    
    radius = maxd /2.
    centerx = (x0+x1)/2. ; centery = (y0+y1)/2.
    center = (centery,centerx)
    circmask = dist_superellipse(image.shape,center,q=1,pos_ang=0,c=0)
    circmask[num.where(circmask<=radius)] = 0
    circmask[num.where(circmask!=0)] = 1
    
    return circmask,center,radius
    
def Ex_DistMatrix(x,y):
    """Gets Matrix of Distances amongst a set of points"""
    n = len(x)
    distances = num.zeros(shape=(n,n),type='Float32')
    
    for i in range(0,n):
        for j in range(i+1,n):
            distances[i,j] = num.sqrt((x[i]-x[j])**2.+(y[i]-y[j])**2.)
    
    return distances
    
