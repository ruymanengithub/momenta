#! /usr/bin/env python

# IMPORT STUFF
import numpy as num
from pdb import set_trace as stop
import sys
# END IMPORT

def gini_core(x):
    """Computes gini parameter."""
    n = len(x)
    x = num.sort(x)
    # G = (<X> n(n-1) )^-1 * SUM(i->n) (2i - n - 1) Xi
    mean_x = x.mean()
    summatory = sum(( 2.* num.arange(1.,n+1.) - n - 1.) * x )
    G = summatory / (mean_x * n * (n-1.))
    return G
    
def gini(self,dograph=False):
    """Returns the gini parameter of a distribution of pixels.
    Algorithm should be checked.
    """
    # IMPORT STUFF
    # END IMPORT
    # INPUTS
    sky = self['BACKGROUND']
    image = self['STAMP'].copy() - sky
    try: mask = self['MASK'].copy()
    except AttributeError: mask = self['MASK']
    # END INPUTS
    
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')    
    active = num.where((mask == 0) & (image > 0.))
    x = image[active]
    # if x.min() < 0. : x -= x.min() # !!!
    # x = x[num.where(x >=0)]
    x = num.abs(x)
    
    G = gini_core(x)
    
    if G > 1. or G < 0. : 
        print 'G = %8.4f, stopping...' % G
        sys.exit()
    
    print 'G = %.2f' % G
    self['M_GINI'] = G
    
    if dograph:
        self.gini_graph(x)
