def FindPeaks(self,norm=-1,dograph=False):
    """Returns number of 'significative peaks' in an image."""
    # IMPORT STUFF
    import numpy as num
    from pdb import set_trace as stop
    from numpy.nd_image import shift
    from numpy.nd_image.filters import uniform_filter
    import pyfits
    import os
    #from Moments.algorithms import get_stat
    # END IMPORT
    
    # INPUTS
    sky = self['BACKGROUND']
    image = self['STAMP'].copy() - sky
    try : mask = self['MASK'].copy()
    except AttributeError : mask = self['MASK']
    sigma_sky = self.execpars['sigma_sky'][0]
    # END INPUTS
    
##    gauss33 = num.array([[0.54,0.73,0.54],[0.73,1.,0.73],[0.54,0.73,0.54]])
##    ngauss33 = gauss33.sum()
    
##    gauss55 = num.array([[0.09,0.21,0.29,0.21,0.09],\
##    [0.21,0.54,0.73,0.54,0.21],\
##    [0.29,0.73,1.,0.73,0.29],\
##    [0.21,0.54,0.73,0.54,0.21],\
##    [0.09,0.21,0.29,0.21,0.09]])
##    ngauss55 = gauss55.sum()
    
##    gauss77 = num.array([[0.004,0.02,0.05,0.06,0.05,0.02,0.004],\
##    [0.02,0.09,0.21,0.29,0.21,0.09,0.02],\
##    [0.05,0.21,0.54,0.73,0.54,0.21,0.05],\
##    [0.06,0.29,0.73,1.0,0.73,0.29,0.06],\
##    [0.05,0.21,0.54,0.73,0.54,0.21,0.05],\
##    [0.02,0.09,0.21,0.29,0.21,0.09,0.02],\
##    [0.004,0.02,0.05,0.06,0.05,0.02,0.004]])
##    ngauss77 = gauss77.sum()
    
    if mask is -1 : mask = num.zeros(shape=image.getshape(),type='Int8')    
    image[num.where(mask != 0)] = 0.
    # active = num.where((mask == 0) & (image > 0.))
    #sigma = get_stat(image,'stddev',minimum=-99)
    #if norm != -1 : image /= norm
    
    filtered3 = num.zeros(shape=image.shape,type='Float32')
    filtered7 = num.zeros(shape=image.shape,type='Float32')
    filtered15 = num.zeros(shape=image.shape,type='Float32')
    
    uniform_filter(image,(3,3),output=filtered3,\
    mode='constant',cval=0)
    uniform_filter(image,(7,7),output=filtered7,\
    mode='constant',cval=0)
    uniform_filter(image,(15,15),output=filtered15,\
    mode='constant',cval=0)
    
    #detect = 100. * (filtered3 - filtered7) / num.abs(filtered3)
    #+num.abs(filtered7))
    
    # detect = (9 * filtered3 - 49 * filtered7) / sigma
    #bigapper = 49*filtered7
    #smallapper = 9 * filtered3
    #detect = smallapper / bigapper
    #detect = filtered3.copy()
    #detect = filtered3 - filtered7
    detect = 49. * filtered7 - (49.*(225.*filtered15-49.*filtered7)/176.)
    #detsigma = get_stat(detect,'stddev',minimum=-99)
    #print sigma, detsigma
    
    maxima = num.ones(shape=image.shape,type='Bool')
    
    #gaussianity = num.zeros(shape=image.shape,type='Float32')
    
    for i in range(-2,3,1):
        for j in range(-2,3,1):
            if i==0 and j==0:
                #gaussianity += image.copy()
                pass
            else:
                tmpshift1 = detect.copy() * 0.
                #tmpshift2 = image.copy() * 0.
                shift(detect,(i,j),output = tmpshift1)
                #shift(image,(i,j),output = tmpshift2)
                maxima = maxima & (detect > tmpshift1)
                #if num.abs(i) <= 1 and num.abs(j) <=1:
                #    maxima = maxima & (tmpshift2 > sigma)
                #maxima = maxima & (num.abs(detect - tmpshift1)< \
                #0.50 * num.abs(detect))
                # maxima = maxima & (image > tmpshift2)
                # gaussianity += (tmpshift2.copy() / gauss33[i+1,j+1])
                #gaussianity += tmpsshift2.copy() / gauss55[i+2,j+2]
    
    #effsigma = 8.
    #gaussianity = (gaussianity / 9.) / filtered3
    #flux3 =  9. * filtered3 - (9.*(49.*filtered7-9.*filtered3)/(49.-9.))
    #relevance = (flux3 / (num.sqrt(9. * effsigma**2.))) >= 3
    #fquot = (9.*filtered3 / (49.*filtered7))
    #relevance = fquot >= 0.2
    #relevance = relevance & (fquot <= 0.5)
    relevance = filtered3 / sigma_sky > 1.
    relevance2 = detect > 3. * 13.13 * sigma_sky #detsigma
    maxima = maxima & relevance & relevance2
    #maxima = maxima & (num.abs(gaussianity-1.5)/1.5 < 0.1)
    # maxima = maxima & (detect > sigma)
    # maxima = maxima & (detect < -10.0)
    # maxima = maxima & (detect > -100.0)
    self['MAXIMA'] = maxima
    self['DETECTIMG'] = detect
    self['M_NPEAKS'] = len(num.where(maxima)[0])
    
    if dograph:
        self.FindPeaks_graph()
    
    return None
    
