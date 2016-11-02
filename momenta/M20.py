#! /usr/bin/env python

"""Different implementations of the M20 parameter.

M20 = 
"""

# IMPORT STUFF
import numpy as num
from pdb import set_trace as stop
import sys
from pyraf import iraf
import pyfits
from time import time
from Moments.graphs import Convert, doStamp
# END IMPORT

def M20Lotz(self,dograph=False):
    """Returns M20 as defined by Lotz et al. 2004
    M20 = M2_20 / M2
    M2 = num.sum(i(x,y) * r**2)
    r = sqrt((x-xcenter)+(y-center))
    """
    # IMPORT STUFF
    import sys
    from flags import addflag,allflags,isflagon
    # END IMPORT
    # INPUTS
    sky = self['BACKGROUND']
    img = self['STAMP'].copy() - sky
    try: mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK'].copy()
    x0 = self['MXMIN_IMAGE']
    y0 = self['MYMIN_IMAGE']
    xcenter = self['X_IMAGE']  ; ycenter = self['Y_IMAGE'] # as in RADIAL
    # END INPUTS
    
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')

    active = num.where((mask == 0) & (img > 0.))
    
    ord_int = num.sort(img[active])
    ord_cumul = num.zeros(shape=ord_int.shape,type='Float32')
    ord_cumul[-1] = ord_int[-1]
    for i in range(len(ord_int)-2,-1,-1):
        ord_cumul[i] = ord_int[i] + ord_cumul[i+1]
    
    ord_cumul = ord_cumul / ord_cumul[0]
    ord_cumul_20 = ord_cumul - 0.2
    
    M20 = -99
    try: 
        def get_active(img,ord_cumul,ord_int,mask):
            # trimmer = num.where(ord_cumul ==
            # ord_cumul[num.where(ord_cumul>0)].min())[0][0]
            trimmer = num.where(ord_cumul > 0)[0].max()
            cut = ord_int[trimmer]
            active = num.where((mask == 0) & (img>cut))
            
            return active,cut
            
        active20,cut20 = get_active(img,ord_cumul_20,ord_int,mask)
        print 'Cut20 at %f, background included' % (cut20+sky,)
        try:
            
            def GetM2(image,active,origin,center):
                i = image[active]
                x = active[1] + origin[1]
                y = active[0] + origin[0]
                r = num.sqrt((x-center[1])**2.+(y-center[0])**2.)
                M2 = num.sum(i * r**2.)
                return M2
            
            origin = (y0,x0)
            center = (ycenter,xcenter)
            
            M2_20 = GetM2(img,active20,origin,center)
            M2_100 = GetM2(img,active,origin,center)
        
        except ZeroDivisionError: raise num.libnumarray.error
        
        if M2_20 <= 0. or M2_100 <= 0.:
            M20 = -99.
            raise num.libnumarray.error
        else: M20 = num.log10(M2_20/M2_100)
        
    except num.libnumarray.error:
        self['flags'] = addflag(self['flags'],allflags['NOM20'])
    
    print 'M20 = %.2f' % M20
    self['M20'] = M20
    self['M_I20'] = cut20
    
    if dograph:
        if not isflagon(self['flags'],allflags['NOM20']):
            self.M20_graph(img,active,active20)
        else:
            pass
    
def M20Azzo(self,dograph=False):
    """Returns M20 in a fashion which is a slight modification to the definition
    by Lotz et al. 2004."""
    # IMPORT STUFF
    import numpy as num
    from pdb import set_trace as stop
    import sys
    from flags import addflag,allflags,isflagon
    # END IMPORT
    # INPUTS
    sky = self['BACKGROUND']
    img = self['STAMP'].copy() - sky
    try: mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK'].copy()
    x0 = self['MXMIN_IMAGE']
    y0 = self['MYMIN_IMAGE']
    # END INPUTS

    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')

    active = num.where((mask == 0) & (img > 0.))
    #x2, y2, xy = moments(img,x0,y0,order=2,active=active)
    
    ord_int = num.sort(img[active])
    ord_cumul = num.zeros(shape=ord_int.shape,type='Float32')
    ord_cumul[-1] = ord_int[-1]
    for i in range(len(ord_int)-2,-1,-1):
        ord_cumul[i] = ord_int[i] + ord_cumul[i+1]
    
    ord_cumul = ord_cumul / ord_cumul[0]
    ord_cumul_20 = ord_cumul - 0.2
    ord_cumul_90 = ord_cumul - 0.9
    M20 = -99
    try: 
        def get_active(img,ord_cumul,ord_int,mask):
            #trimmer = num.where(ord_cumul ==
            #ord_cumul[num.where(ord_cumul>0)].min())[0][0]
            trimmer = num.where(ord_cumul > 0)[0].max()
            cut = ord_int[trimmer]
            active = num.where((mask == 0) & (img>cut))
            
            return active,cut
            
        active20,cut20 = get_active(img,ord_cumul_20,ord_int,mask)
        active90,cut90 = get_active(img,ord_cumul_90,ord_int,mask)    
        print 'Cut20 at %f, Cut90 at %f, background included' % \
        (cut20+sky,cut90+sky)
        try:
            x2_20, y2_20, xy_20 = moments(img,x0,y0,order=2,active=active20)
            x2_90, y2_90, xy_90 = moments(img,x0,y0,order=2,active=active90)
        except ZeroDivisionError: raise num.libnumarray.error
        #M2 = num.sqrt(num.abs(x2)+num.abs(y2))
        M2_20 = num.sqrt(num.abs(x2_20)+num.abs(y2_20))
        M2_90 = num.sqrt(num.abs(x2_90)+num.abs(y2_90))
        
        M20 = num.log10(M2_20/M2_90)
        
    except num.libnumarray.error:
        self['flags'] = addflag(self['flags'],allflags['NOM20'])
    
    print 'M20 = %4.2f' % M20
    self['M20'] = M20
    self['M_I20'] = cut20
    
    if dograph:
        if not isflagon(self['flags'],allflags['NOM20']):
            self.M20_graph(img,active90,active20)
        else:
            pass
    
