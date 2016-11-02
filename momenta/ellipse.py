#! /usr/bin/env python

import numpy as num

def dist_superellipse(n,center,q=1,pos_ang=num.array(0),c=num.array(0)):

    """Form an array in which the value of each element is equal to the \n
    semi-major axis of the superellipse of specified center, axial ratio,\n
    position  angle, and c parameter which passes through that element.\n  
    Useful for super-elliptical aperture photometry.
    
    Inspired on dist_ellipse.pro from AstroLib (IDL).
    
    Note: this program doesn't take into account the change in the order of axes
    from IDL to Python. That means, that in 'n' and in 'center', the order of the
    coordinates must be reversed with respect to the case for dist_ellipse.pro, in
    order to get expected results. Nonetheless, the polar angle means the 
    counter-clock wise angle with respect to the 'y' axis.
    
    """

    # IMPORT STUFF
    from pdb import set_trace as stop
    # END
    
    # CHECK INPUTS
    
    criterio1 = type(n) in (type(1),type(1.0)) or \
    type(n) == type(()) and len(n) == 2
    
    if not criterio1: print 'n must be an scalar or a 2 elements tuple'
    
    if type(n) in (type(1),type(1.0)) : n = float(n)    
    if type(n) == 'tuple' : n = tuple([float(x) for x in n])

    criterio2 = 'int' in str(type(q)) or 'float' in str(type(q))
    if not criterio2: 
        print 'q must be an integer or float'
        criterio3 = False
    else: 
        criterio3 = q>=0 and q<=1
    
    if not criterio3: print 'q = %f is out of range (0,1)' % q

    if criterio2 and criterio3: q = float(q)

    criterio4 = 'int' in str(type(pos_ang)) or 'float' in str(type(pos_ang))
    if not criterio4:
        print 'pos_ang must be an integer or float'
        criterio5 = 'False'
    else:
        criterio5 = -180 <= pos_ang <= 180 
    
    if not criterio5: print 'pos_ang out of range (-180,180)'

    if criterio4 and criterio5: pos_ang = float(pos_ang)
    
    criterio6 = 'int' in str(type(c)) or 'float' in str(type(c))
    if not criterio6: print 'c must be an integer or float'	     
    else : c = float(c)

    criterio7 = type(center) == type(()) and len(center) == 2
    if not criterio7: print 'center must be a 2 element tuple'
    else : center = tuple([float(x) for x in center])
        
    # END CHECK OF INPUTS
        
    criterios = reduce(num.logical_and,num.array([criterio1,criterio2,criterio3,\
    criterio4,criterio5,criterio6,criterio7]))
    
    if criterios : 
    
        radeg = 180. / num.pi
    
        ang = pos_ang / radeg
        cosang = num.cos(ang)
        sinang = num.sin(ang)
    
        if len(n) == 2 : 
            nx = n[1] 
            ny = n[0]
        else : 
            nx = ny = n
        
        x = num.arange(nx,dtype='Float32') - center[1]
        y = num.arange(ny,dtype='Float32') - center[0]
        im = num.zeros(shape=(ny,nx),dtype='Float32')
        xcosang = x * cosang
        xsinang = x * sinang
        
        for i in range(ny) :
            xtemp = xcosang + y[i] * sinang
            ytemp = -xsinang + y[i] * cosang
            im[i,:] = ( (num.abs(xtemp/q))**(c+2.) + \
			(num.abs(ytemp))**(c+2.) )**(1./(c+2.))

    else: im = num.array([0.0])
    
    return im

    
def area_superellip(r,q,c=0):
    """Returns area of superellipse, given the semi-major axis length"""
    # IMPORT STUFF
    import numpy as num
    from scipy.special import gamma
    #from pdb import set_trace as stop
    # END IMPORT
    
    a_dummie = 4.**(1. - (c+2)**(-1.))
    b_dummie = r**2. * q * ((num.pi)**(0.5))
    c_dummie = gamma(1. + (c+2.)**(-1.)) / gamma(0.5 + (c+2.)**(-1.))
    area = a_dummie * b_dummie * c_dummie
    return area
    
    
def effective_radius(area,q=0,c=0):
    """Returns semi-major axis length of superellipse, given the area"""
    # IMPORT STUFF
    import numpy as num
    from scipy.special import gamma
    #from pdb import set_trace as stop
    # END IMPORT
    a_dummie = 4.**(1. - (c+2.)**(-1.))
    c_dummie = gamma(1. + (c+2.)**(-1.)) / gamma(0.5 + (c+2.)**(-1.))
    
    b_dummie = area / (a_dummie * c_dummie)
    r = (b_dummie/(q * (num.pi ** 0.5)))**(0.5)
    
    return r
