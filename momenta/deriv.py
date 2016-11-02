#! /usr/bin/env python


def deriv(x,y=None):
    """Perform numerical differentiation using 3-point, Lagrangian 
    interpolation. It is a "translation" of deriv.pro from IDL 6.2
    
    CALLING SEQUENCE:
       Dy = Deriv(Y)	           ; Dy(i)/di, point spacing = 1.
       Dy = Deriv(X, Y)        ; Dy/Dx, unequal point spacing.

    INPUTS:
       Y:  Variable to be differentiated. (numpy array)
       X:  Variable to differentiate with respect to.  If omitted, unit 
        spacing for Y (i.e., X(i) = i) is assumed.
           (numpy array)
    
    OPTIONAL INPUT PARAMETERS:
       As above.

    OUTPUTS:
       Returns the derivative.

    PROCEDURE:
       See Hildebrand, Introduction to Numerical Analysis, Mc Graw
       Hill, 1956.  Page 82.

    MODIFICATION HISTORY:
       Written, DMS, Aug, 1984 for IDL
       Corrected formula for points with unequal spacing.  DMS, Nov, 1999.
       
       Translated to Python by Sir Ruyman Azzollini, Oct 2005.
    
    """
    # IMPORT STUFF
    import numpy as num
    import sys
    from pdb import set_trace as stop
    # END IMPORT 
    
    def shift(array,inc):
        import numpy as num
        shifted = num.concatenate((array[-inc:],array[0:-inc]),axis=0)
        return shifted

    dummie = num.array([1])

    n = len(x)
    if n < 3 : 
        print 'Parameters must have at least 3 points'
        raise RuntimeError
        sys.exit()
    
    if type(x) != type(dummie) :
        print "input vectors must be of numpy.ndarray type"
        raise RuntimeError
        sys.exit()

    if y!=None:
        if type(y) != type(dummie):
            print "input vectors must be of numpy.ndarray type"
            raise RuntimeError
            sys.exit()
        if len(y) != n : 
            print 'Vectors must have same size'
            raise RuntimeError
            sys.exit()
        
        # df/dx = y0*(2x-x1-x2)/(x01*x02)+y1*(2x-x0-x2)/(x10*x12)+
        # y2*(2x-x0-x1)/(x20*x21) 
        # Where: x01 = x0-x1, x02 = x0-x2, x12 = x1-x2, etc.
        
        typein = x.type()
     
     	if typein != 'Float32': x = x.astype('Float32')
     
     	x12 = x - shift(x,-1) # x1 - x2
     	x01 = shift(x,1) - x # x0 - x1
     	x02 = shift(x,1) - shift(x,-1) # x0 - x2
     
        d = shift(y,1) * (x12 / (x01 * x02)) + \
        y * (1./x12 - 1./x01) - shift(y,-1) * (x01 / (x02 * x12)) 
     
        # Formulae for the first and last points:
        # First point
        d[0] = y[0] * (x01[1]+x02[1])/(x01[1]*x02[1]) - \
        y[1] * x02[1]/(x01[1]*x12[1]) + \
        y[2] * x01[1]/(x02[1]*x12[1])
     
        n2 = n-2
        # Last point
        	
        d[n-1] = -y[n-3] * x12[n2]/(x01[n2]*x02[n2]) + \
        y[n-2] * x02[n2]/(x01[n2]*x12[n2]) - \
        y[n-1] * (x02[n2]+x12[n2]) / (x02[n2]*x12[n2])
 
    else : # Equally spaced point case
        d = (shift(x,-1) - shift(x,1))/2.	    	    
        d[0] = (-3.0*x[0] + 4.0*x[1] - x[2])/2.
        d[n-1] = (3.*x[n-1] - 4.*x[n-2] + x[n-3])/2.

    return d
    
