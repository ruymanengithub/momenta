#! /usr/bin/env python

"""Adaptation to Python of 'ladfit.pro', by Ruyman Azzollini

This function fits the paired data (X(i), Y(i)) to the linear model,
y = A + Bx, using a 'robust' least absolute deviation method. The
result is a two-element vector containing the model parameters, A
and B.

"""

# IMPORT STUFF
from numpy import median,std
import numpy as num
import sys
from pdb import set_trace as stop
# END IMPORT

def MDfunc(b,x,y,a,absdev,eps):
    """Function called by LadFit"""
    a[0] = float(median(y - b*x))
    d = y - (b*x +a[0])
    
    
    absdev[0] = num.sum(num.abs(d))
    nz = num.where(y != 0.0)
    nzcount = len(nz[0])
    if nzcount != 0: d[nz] = d[nz] / abs(y[nz]) # Normalize
    nz = num.where(num.abs(d) > eps)
    nzcount = len(nz[0])
    
    if nzcount != 0:
        return num.sum(x[nz] * ((d[nz]>0.) - (d[nz]<0.)))
    else: return 0.0

def LadFit(x,y,absdev=[0.],sigma=[0.],maxiter=1000):
    """Function to be called
    results = [y0,slope]
    absdev is a list... take care."""
    
    x = num.array(x,type='Float64')
    y = num.array(y,type='Float64')
    
    nX = len(x)
    
    if nX != len(y) :
        sys.exit('X and Y must be vectors of equal length.')
    
    sx = num.sum(x)
    sy = num.sum(y)
    
    sxy = num.sum(x*y)
    sxx = num.sum(x*x)
    del_ = float(nX) * sxx - sx**2.
    
    if (del_ == 0.) :              # All X's are the same
        result = [median(y),0.0]   # Bisect the range w/ a flat line
        return result
    
    aa = (sxx * sy - sx * sxy) / del_ # Least Squares solution y = x * aa + bb
    bb = (float(nX) * sxy - sx * sy) / del_
    chisqr = num.sum((y-(aa+bb*x))**2.)
    sigb = num.sqrt(chisqr/del_) # Standard deviation
    
    b1 = bb
    eps = 1.E-7
    
    aa = [aa]
    f1 = MDfunc(b1,x,y,aa,absdev,eps=eps)
    
    
    done = False
    
    # Quick return. The initial least squares gradient is the LAD solution.
    
    if (f1 == 0.):
        bb = b1
        done = True
    
    quit = 0
    
    if not done:
        if f1 >= 0.: delb = 3. * sigb
        else : delb = -3. * sigb
     
        b2 = b1 + delb
        
        f2 = MDfunc(b2,x,y,aa,absdev,eps=eps)
        
	niter1 = 1
        while (f1*f2 > 0.) : # Bracket the zero of the function
            b1 = b2
            f1 = f2
            b2 = b1 + delb
            f2 = MDfunc(b2,x,y,aa,absdev,eps=eps)
	    niter1 += 1
	    if niter1 > maxiter:
	        f1 = 1. ; f2 = -1.
		quit = 1
	
	if quit != 1:	

            # In case we finish early.
            bb = b2
            f = f2

            # Narrow tolerance to refine 0 of fcn.
            sigb = 0.01 * sigb
	
            niter2 = 1
            while ((abs(b2-b1) > sigb) and (f != 0.)):
                bb = 0.5 * (b1+b2)
                if ((bb == b1) or (bb == b2)): break
    
                f = MDfunc(bb,x,y,aa,absdev,eps=eps)
    
                if (f*f1 >= 0):
                    f1 = f
                    b1 = bb
                else:
                    f2 = f
                    b2 = bb
	    niter2 += 1
	    if niter2 > maxiter:
	        f = 0.
		quit = 1

    done = True
    
    if quit : return [None,None]
    else:
        if done : absdev[0] = absdev[0] / float(nX)
        sigma[0] = float(std(y-(aa[0]+bb*x)))
        return [aa[0],bb]
