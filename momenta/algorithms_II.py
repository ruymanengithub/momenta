#! /usr/bin/env python

import numpy as num
from pdb import set_trace as stop

def rot180numpy(image,center):
    """Rotates an image 180 degrees respect to a center, using numpy."""    
    # IMPORT STUFF
    from algorithms_II import shrotate
    from time import time
    # END IMPORT
    #print time()
    
    #shape = image.shape
    xin = center[1] + 1. ; yin = center[0] + 1. # RIGHT?
    
    #xout = (shape[1]-1.)/2. + 1
    #yout =  (shape[0]-1.)/2. + 1

    #nullangle = 0.
    #image1 = shrotate(image.copy(), nullangle, xout, yout, xin, yin, order = 3,\
    #mode = 'constant', cval = 0.0,prefilter = False)
    angle = 180.
    image1 = image.copy()
    #image2 = shrotate(image1, angle, xout, yout, order = 3,\
    #mode = 'constant', cval = 0.0,prefilter = False)
    image2 = shrotate(image1, angle, xin, yin, order = 3,\
    mode = 'constant', cval = 0.0,prefilter = False)
    
    return image1,image2


def bend180numpy(image,x,y,pa):
    """Bends an image over an axis, using numpy."""
    # IMPORT STUFF
    from algorithms_II import shrotate
    # END IMPORT
    
    angle = pa
    xin = x + 1. ; yin = y + 1. # RIGHT?
    imagerot = shrotate(image, angle, xin, yin, prefilter = False)
    
    indexes = num.where(image == image)
    
    y_indexes = indexes[0].copy()
    
    y_indexes = 2* num.around(num.array([y]))[0] - y_indexes
    
    inside = num.where((y_indexes > 0) & (y_indexes < len(indexes[0])))
    
    bindexes = (y_indexes[inside],indexes[1][inside].copy())
    oindexes = (indexes[0][inside].copy(),indexes[1][inside].copy())
    
    imagebend = num.zeros(image.shape,dtype='Float32')
    
    imagebend[bindexes] = imagerot[oindexes]
    
    return imagerot, imagebend
    
    
def shrotate(input, angle, xin, yin, xout = None, yout = None,
             axes = (-1, -2), reshape = False, output_type = None,
             output = None, order = 3, mode = 'constant', cval = 0.0,
             prefilter = True):


    import math

    """Rotate and shift an array around a center.

    Adapted from rotate in numpy, to include a rotation center
    defined by (xin,yin) and and optional output center (xout,yout)
    F Menanteau, JHU. May 2005.

    The array is rotated in the plane definde by the two axes given by
    the axes parameter using spline interpolation of the requested
    order. The angle is given in degrees. Points outside the
    boundaries of the input are filled according to the given
    mode. The output type can optionally be given. If not given it is
    equal to the input type. If reshape is true, the output shape is
    adapted so that the input array is contained completely in the
    output. Optionally an output array can be provided that must match
    the requested output shape and type. The parameter prefilter
    determines if the input is pre-filtered before interpolation, if
    False it is assumed that the input is already filtered.
    
    NOTE (by Azzollini): The angle is clock-wise.
    
    """
    
    input = numpy.asarray(input)

    angle = numpy.pi / 180 * angle

    if axes[0] < axes[1]:
        a1 = axes[0]
        a2 = axes[1]
    else:
        a1 = axes[1]
        a2 = axes[0]
        
    oshape = list(input.shape)
    ix = input.shape[a1]
    iy = input.shape[a2]
    if reshape:
        # Fix, we now take abs value to avoid crash, when giving
        # negative integer values
        ox = abs(ix * math.cos(angle) + iy * math.sin(angle) + 0.5)
        oy = abs(iy * math.cos(angle) + ix * math.sin(angle) + 0.5)
        ox = int(ox)
        oy = int(oy)
        oshape[a1] = ox
        oshape[a2] = oy
    else:
        ox = oshape[a1]
        oy = oshape[a2]

    m11 = math.cos(angle)
    m12 = math.sin(angle)
    m21 = -math.sin(angle)
    m22 = math.cos(angle)
    matrix = numpy.identity(input.rank, type = numpy.float64)
    matrix[a1, a1] = m11
    matrix[a1, a2] = m12
    matrix[a2, a1] = m21
    matrix[a2, a2] = m22

    # Fix the indices to start from 1 rather than zero
    xin = xin - 1
    yin = yin - 1

    # Use inputs if output centers not defined
    if xout and yout:
        xout = xout - 1
        yout = yout - 1
    else:
        xout = xin
        yout = yin
     
    offset = numpy.zeros((input.rank,), dtype = numpy.float64)
    # offset is just xo,yo
    offset[a1] = yin
    offset[a2] = xin
    # Multiply it by the rotation matrix
    offset = numpy.matrixmultiply(matrix, offset)

    # New center 
    tmp = numpy.zeros((input.rank,), dtype = numpy.float64)
    tmp[a1] = yout
    tmp[a2] = xout
    offset = tmp - offset
    
    return numpy.nd_image.affine_transform(input, matrix, offset, oshape,\
    output_type,output, order, mode, cval, prefilter)
