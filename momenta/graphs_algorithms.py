#! /usr/bin/env python

from CommonTools.graphs import Convert, Plot, doStamp, Paint
from AxisAsymmetry import bend180iraf
from Asymmetry import rot180iraf
import numpy as num
import os
from pdb import set_trace as stop
from momlib import quitpath
import string
commtextpos = (20,20)

def _getGraphId(self):
    if self.execpars['father_img_name'] != 'None':
        imgid = quitpath(self.execpars['father_img_name'])
    else:
        imgid = quitpath(self['IMAGE'])
    imgid = imgid[0:string.rfind(imgid,'.')]
    name = '%s' % self['name']
    commid = '%s_%s' % (name,imgid)
    return commid

def gini_graph(self,x):
    """Does graph for Gini"""
    
    name = self['name']
    id = self._getGraphId()
    figname = 'G_%s.eps' % id
    G = self['M_GINI']
    sxlabel = 'Pixel Order' ; sylabel = 'Cumul Flux' 
    title = 'Gini, %s, G=%4.2f' % (name,G)
    x = num.sort(x)
    cumul = num.cumsum(x)
    lorenz_y = 100.* cumul/cumul[-1]
    lorenz_x = 100. * (num.arange(len(x))+1)/float(len(x))
    fakearray = num.array([0.,100.])
    xy = ((lorenz_x,lorenz_y),(fakearray,fakearray))
    Plot(xy,figname,sxlabel,sylabel,title)
    self['figures']['gini'] = figname

def M20_graph(self,img,active_90,active_20):
    """Does Graph for M20"""
    fake = num.zeros(shape=img.shape,type='Float32')
    fake[active_90] = 1. ; fake[active_20] = 2.
    id = self._getGraphId()
    figname = 'G_%s.eps' % id
    root = 'M20_%s' % id
    pngname = root + '.png' ; jpgname = root + '.jpg'
    epsname = root + '.eps'
    doStamp(fake,pngname,format='PNG',minhisto=0,maxhisto=100)
    Convert(pngname,jpgname)
    Painted = Paint(jpgname)
    Painted.load()
    text = 'M20=%5.2f' % self['M20']
    #Painted.Graffiti(text,commtextpos)
    Painted.save(jpgname)
    Painted.release()
    
    Convert(jpgname,epsname)
    os.system('rm %s %s' % (pngname,jpgname))
    self['figures']['M20'] = epsname
    self['figcomms']['M20'] = text

def axis_asymmetry_graph(self,label,masked,pa,Ax_center):
    """Does graph for axis asymmetry"""
    from algorithms import bend180iraf
    x = Ax_center[1] ; y = Ax_center[0]
    image1, image2 = bend180iraf(masked.copy(),x,y,pa)
    AxImg = num.abs(image1 - image2)
    
    id = self._getGraphId()
    root = 'Ax_%s' % (label,id)
    pngname = root + '.png' ; epsname = root + '.eps'
    jpgname = root + '.jpg'
    doStamp(AxImg,pngname,format='PNG')
    Convert(pngname,jpgname)
    
    Painted = Paint(jpgname)
    Painted.load()
    text = 'AxsAs%s=%5.2f' % (label,self['M_AXS%s'%label])
    #Painted.Graffiti(text,commtextpos)
    Painted.save(jpgname)    
    Painted.release()
    
    Convert(jpgname,epsname)
    os.system('rm %s' % (pngname,jpgname))
    self['figures']['Ax%s'%label] = epsname
    self['figcomms']['Ax%s'%label] = text
    
def asymmetry_graph(self,masked,A_center):
    """Does graph for asymmetry"""
    from algorithms import rot180iraf
    
    image1, image2 = rot180iraf(masked.copy(),A_center)
    AImg = num.abs(image1.copy() - image2.copy())
    
    id = self._getGraphId()
    root = 'AS_%s' % (id,)
    pngname = root + '.png' ; epsname = root + '.eps'
    jpgname = root + '.jpg'
    doStamp(AImg,pngname,format='PNG')
    Convert(pngname,jpgname,cleanafter=False)
    
    Painted = Paint(jpgname)
    
    Painted.load()
    
    text = 'A=%5.2f' % (self['M_AS'])
    # Painted.Graffiti(text,commtextpos)
    Painted.save(jpgname)    
    Painted.release()
    
    Convert(jpgname,epsname)
    os.system('rm %s %s' % (pngname,jpgname))
    self['figures']['A'] = epsname
    self['figcomms']['A'] = text
    
def concent_graph(self):
    """Does graph for concentration"""
    r_big = self['M_RSMALL']
    r_small = self['M_RBIG']
    C = self['M_C']
    
    xcenter = self['X_IMAGE']  ; ycenter = self['Y_IMAGE']
    xcenter = xcenter - self['MXMIN_IMAGE']
    ycenter = ycenter - self['MYMIN_IMAGE']
    center = (xcenter,ycenter)
    
    ellip = self['ELLIPTICITY'] 
    q = 1. - ellip
    pa = self['THETA_IMAGE']  # Astronomical position angle.

    stamp = self['STAMP'].copy()
    mask = self['MASKOTHER'].copy()
    sky = self['BACKGROUND']
    Img = stamp - sky
    Img[num.where(mask != 0)] = 0.
    
    id = self._getGraphId()
    root = 'C_%s' % (id,)
    pngname = root + '.png' ; epsname = root + '.eps'
    jpgname = root + '.jpg'
    doStamp(Img,pngname,format='PNG')
    Convert(pngname,jpgname)
    
    Painted = Paint(jpgname)
    Painted.load()
    Painted.DrawEllipse(center,r_big,q,pa,color='red',linewidth=2)
    Painted.DrawEllipse(center,r_small,q,pa,color='green',linewidth=2)
    
    text = 'C=%5.2f' % (self['M_C'])
    # Painted.Graffiti(text,commtextpos)
    Painted.save(jpgname)    
    Painted.release()
    
    Convert(jpgname,epsname)
    os.system('rm %s %s' % (pngname,jpgname))
    
    self['figures']['C'] = epsname
    self['figcomms']['C'] = text
    
def clumpy_graph(self,img):
    """Does graph for clumpyness"""
    
    id = self._getGraphId()
    root = 'S_%s' % (id,)
    pngname = root + '.png' ; epsname = root + '.eps'
    jpgname = root + '.jpg'
    doStamp(img.copy(),pngname,format='PNG')
    Convert(pngname,jpgname)
    
    Painted = Paint(jpgname)
    try : Painted.load()
    except IOError : stop()
    text = 'S=%5.2f' % (self['M_S'])
    # Painted.Graffiti(text,commtextpos)
    Painted.save(jpgname)    
    Painted.release()
    
    Convert(jpgname,epsname)
    os.system('rm %s %s' % (pngname,jpgname))
    self['figures']['S'] = epsname
    self['figcomms']['S'] = text
    
def petropeakcenter_graph(self):
    """Does graph for Petrosian Radius, Peak, and Center"""
    
    R_PETRO = self['R_PETRO']
    
    xcenter = self['X_IMAGE']  ; ycenter = self['Y_IMAGE']
    xcenter = xcenter - self['MXMIN_IMAGE']
    ycenter = ycenter - self['MYMIN_IMAGE']
    center = (xcenter,ycenter)
    
    ellip = self['ELLIPTICITY'] 
    q = 1. - ellip
    pa = self['THETA_IMAGE']  # Astronomical position angle.
    
    if 'M_XPEAK' in self:
        xpeak = self['M_XPEAK'] - self['MXMIN_IMAGE']
        ypeak = self['M_YPEAK'] - self['MYMIN_IMAGE']
        peak = (xpeak,ypeak)
    
    stamp = self['STAMP'].copy()
    mask = self['MASK'].copy()
    sky = self['BACKGROUND']
    Img = stamp - sky
    Img[num.where(mask != 0)] = 0.
    
    id = self._getGraphId()
    root = 'petropeakcenter_%s' % (id,)
    pngname = root + '.png' ; epsname = root + '.eps'
    jpgname = root + '.jpg'
    doStamp(Img,pngname,format='PNG')
    Convert(pngname,jpgname)
    
    Painted = Paint(jpgname)
    Painted.load()
    if R_PETRO >0. : 
        Painted.DrawEllipse(center,R_PETRO,q,pa,color='red',linewidth=2)
    Painted.DrawCross(center,length=20,color='green')
    try : Painted.DrawCross(peak,length=20,color='red')
    except : pass
    
    text = 'RP=%6.2f' % R_PETRO
    
    # Painted.Graffiti(text,commtextpos)
    
    Painted.save(jpgname)    
    Painted.release()
    
    Convert(jpgname,epsname)
    os.system('rm %s %s' % (pngname,jpgname))
    
    self['figures']['petropeakcenter'] = epsname
    self['figcomms']['petropeakcenter'] = text
    
    
def radial_graph(self):
    """Does graph for radial"""
    
    if self['M_RADIAL']['intens'] != None:
        name = self['name']
        id = self._getGraphId()
        figname = 'RADIAL_%s.eps' % id
        sxlabel = 'Pixel Radius' ; sylabel = 'Intens' 
        title = 'Radial profile, %s' % (name,)
        y = self['M_RADIAL']['intens']
        x = self['M_RADIAL']['radii']
        xy = ((x,y),)
        Plot(xy,figname,sxlabel,sylabel,title)
        self['figures']['radial'] = figname
    else : pass

def FindPeaks_graph(self):
    """Does graph for FindPeaks"""
    import string
    
    maxima = self['FP_LOC'].copy()
    maxima = num.where(maxima)
    maxima = (maxima[1],maxima[0])
    detectimg = self['FP_DETECT'].copy()
    
    id = self._getGraphId()
    root = 'FindPeaks_%s' % (id,)
    pngname = root + '.png' ; epsname = root + '.eps'
    jpgname = root + '.jpg'

    doStamp(detectimg,pngname,format='PNG')
    Convert(pngname,jpgname)
    
    Painted = Paint(jpgname)
    Painted.load()
    Painted.DrawCross(maxima,length=7,color='green')
    
    strpeaks = string.strip('%i'% (self['M_NPEAKS']))
    text = 'NP=%s' % strpeaks 
    
    # Painted.Graffiti(text,commtextpos)
    
    Painted.save(jpgname)
    Painted.release()
    
    Convert(jpgname,epsname)
    os.system('rm %s %s' % (pngname,jpgname))
    self['figures']['FindPeaks'] = epsname
    self['figcomms']['FindPeaks'] = text

def FindClumps_graph(self):
    """Does graph for FindClumps"""
    # IMPORT STUFF
    import string
    # END IMPORT
    
    maxima = self['CL_LOC'].copy()
    maxima = num.where(maxima)
    maxima = (maxima[1],maxima[0])
    detectimg = self['STAMP'].copy()
    
    id = self._getGraphId()
    root = 'FindClumps_%s' % (id,)
    pngname = root + '.png' ; epsname = root + '.eps'
    jpgname = root + '.jpg'

    doStamp(detectimg,pngname,format='PNG')
    Convert(pngname,jpgname)
    
    Painted = Paint(jpgname)
    Painted.load()
    Painted.DrawCross(maxima,length=7,color='green')
    
    strpeaks = string.strip('%i'% (self['M_NUM_CL']))
    text = 'NC=%s' % strpeaks 
    
    # Painted.Graffiti(text,commtextpos)
    
    Painted.save(jpgname)
    Painted.release()
    
    Convert(jpgname,epsname)
    os.system('rm %s %s' % (pngname,jpgname))
    self['figures']['FindClumps'] = epsname
    self['figcomms']['FindClumps'] = text
    
def stamp_graph(self,key,label):
    """Does stamp graph"""
    stamp = self[key].copy()
    id = self._getGraphId()
    root = '%s_%s' % (label,id)
    pngname = root + '.png'
    epsname = root + '.eps'
    if label == 'mask': 
        doStamp(stamp,pngname,format='PNG',minhisto=0,maxhisto=100)
    else : doStamp(stamp,pngname,format='PNG')
    
    Convert(pngname,epsname)
    os.system('rm %s' % (pngname,))
    self['figures'][label] = epsname

    
def Excentricity_graph(self,ccenter,cradius,icenter):
    """Does Excentricity graph. Draw the circle that encloses the object, its 
    center, and the center of brightness (first order moment)."""
    
    ellip = 0.
    q = 1. - ellip
    pa = 0.
    
    stamp = self['STAMP'].copy()
    mask = self['MASK'].copy()
    sky = self['BACKGROUND']
    Img = stamp - sky
    Img[num.where(mask != 0)] = 0.
    
    id = self._getGraphId()
    root = 'Excentricity_%s' % (id,)
    pngname = root + '.png' ; epsname = root + '.eps'
    jpgname = root + '.jpg'
    doStamp(Img,pngname,format='PNG')
    Convert(pngname,jpgname)
    
    Painted = Paint(jpgname)
    Painted.load()
    
    Painted.DrawEllipse(ccenter,cradius,q,pa,color='red',linewidth=2)
    Painted.DrawCross(ccenter,length=20,color='red')
    Painted.DrawCross(icenter,length=20,color='green')
    
    text = 'E=%.2f o/o' % self['M_E']
    
    # Painted.Graffiti(text,commtextpos)
    
    Painted.save(jpgname)    
    Painted.release()
    
    Convert(jpgname,epsname)
    os.system('rm %s %s' % (pngname,jpgname))
    
    self['figures']['Excentricity'] = epsname
    self['figcomms']['Excentricity'] = text

    
def FFactor_graph(self):
    """Does FFactor graph"""
    pass
def AC_graph(self,image,ACmask):
    """Does AC graph"""
    
    stop()
    image[num.where(ACmask==0)] = 0.
    
    id = self._getGraphId()
    root = 'AC_%s' % (id,)
    pngname = root + '.png' ; epsname = root + '.eps'
    jpgname = root + '.jpg'
    doStamp(image,pngname,format='PNG')
    Convert(pngname,jpgname)
    
    Painted = Paint(jpgname)
    Painted.load()
    
    text = 'AC8=%.2f' % self['M_AC8']
    
    # Painted.Graffiti(text,commtextpos)
    
    Painted.save(jpgname)
    Painted.release()
    
    Convert(jpgname,epsname)
    os.system('rm %s %s' % (pngname,jpgname))
    
    self['figures']['AC'] = epsname
    self['figcomms']['AC'] = text
