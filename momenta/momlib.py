#! /usr/bin/env python 

import numpy as num
from pdb import set_trace as stop

def mergemasks(masks):
    """Merges masks in input tuple"""
    # IMPORT STUFF
    # END IMPORT
    
    shapes = map(num.shape,masks)
    criterio1 = num.all(map(lambda x : x == masks[0].shape,shapes))
    if criterio1 == 1 : criterio1 = True
    else: criterio1 = False
    
    types = map(type,masks)
    dummie2 = num.array([0])
    criterio2 = num.all(map(lambda x : x == type(dummie2),types))
    if criterio2 == 1 : criterio2 = True
    else: criterio2 = False
    
    if criterio1 and criterio2 :
        allmasks = masks[0].copy() * 0
        for mask in masks : allmasks += mask
        nonactive = num.where(allmasks != 0)
        allmasks[nonactive] = 1
        allmasks = allmasks.astype('Int32')
        return allmasks
    else : return None
    

def make_stamp(self,imgname=None,img=None,name='Unknown',
    mode='pyfits',extension=0):
    """Returns a stamp of an image, given the section."""
    # IMPORT STUFF
    import pyfits
    from pyraf import iraf
    from pdb import set_trace as stop
    import sys
    import os
    from time import time
    import string
    # END IMPORT STUFF
    
    xmin = self['MXMIN_IMAGE'] ; xmax = self['MXMAX_IMAGE']
    ymin = self['MYMIN_IMAGE'] ; ymax = self['MYMAX_IMAGE']
    isthere = os.path.exists
    
    tmpcut = 'tmpstamp%s.fits' % string.replace('%f'%time(),'.','_')
    if isthere(tmpcut) : os.system('rm %s'%tmpcut)
    
    if img == None: 
        if mode == 'pyfits':
            stamp = pyfits.getdata(imgname, \
            ext=extension).astype('Float32')[ymin-1:ymax-1+1,\
            xmin-1:xmax-1+1]
        elif mode == 'fitscut':
            x0 = xmin ; y0 = ymin
            ncols = xmax - xmin + 1
            nrows = ymax - ymin +1
            os.system('fitscut --x0=%i --y0=%i -r %i -c %i %s %s' \
            % (x0-1,y0-1,nrows,ncols,imgname,tmpcut))
            try: stamp = pyfits.getdata(tmpcut).astype('Float32')
            except IndexError : stop()
        elif mode == 'pyraf':
            workname = imgname+('[%i:%i,%i:%i]'%(xmin,xmax,ymin,ymax))
            iraf.imcopy(input=workname,output=tmpcut,verbose='no',mode='h')
            try: stamp = pyfits.getdata(tmpcut).astype('Float32')
            except IndexError : stop()
    elif img != None : 
        stamp = img[ymin-1:ymax-1+1,xmin-1:xmax-1+1]
    else :
        print 'There is no input image. What are you doing?'
        print 'Execution halted at make_stamp'
        sys.exit()
    
    if isthere(tmpcut) : os.system('rm %s'%tmpcut)
    
    # stamp = num.array(stamp)
    self[name] = stamp.copy()
    
    # return stamp
    

def make_mask(self,maskstamp,name,mask_in=None,mode='nosky'):
    """Makes a mask image for 'maskstamp'. 
    
    if mode is 'nosky', returned array has 0 on those pixels which 'are' idnum and 1 on the rest.
    if mode is 'withsky', returned array has 0 on those pixels which idnum or 0, and zero elsewhere.
    """
    
    #IMPORT STUFF
    #END IMPORT STUFF
    
    idnum = self['NUMBER']
    
    if mask_in is not None:
        mask = mask_in
    else : mask = num.ones(shape=maskstamp.shape,dtype='UInt8')
    
    if mode == 'nosky':
        indexes = num.where(maskstamp == idnum)
    if mode == 'withsky':
        issue1 = maskstamp == idnum
        issue2 = maskstamp < 1
        indexes = num.where(issue1 | issue2)
        
    mask[indexes] = 0
    mask[num.where(mask != 0)] = 1
    
    self[name] = mask
    
    # return mask
    
    

def wrap2mef(self) :
    """Creates a MEF file with stamp and mask for an object"""
 
    # header: MXMIN_IMAGE, MYMIN_IMAGE, MXMAX_IMAGE, MYMAX_IMAGE, 
    # IMAGENAME, CATNAME, SEGMAPNAME, OUTFILE, DATE, version

    # IMPORT STUFF
    import pyfits
    import time as t
    from os import access, F_OK
    import string 
    from momlib import quitpath
    from CommonTools.ds9 import ds9class
    # END IMPORT STUFF
    
    father_img_name = self.execpars['father_img_name']
    father_img_name_root = string.replace(father_img_name,'.fits','')
    appendix = str(self['NUMBER'])
    name = self['name']
    try : tag = self.execpars['tag'][0]
    except KeyError: tag = ''
    outname = '%s_MEF%s.fits' % (name,tag)
    
    stamp = self['STAMP'].copy()                # numpy
    segstamp = self['SEGSTAMP'].copy()     # numpy
    sexmask = self['SEXMASK'].copy()        # numpy
    maskother = self['MASKOTHER'].copy() # numpy
    
    xmin_image, ymin_image, xmax_image, ymax_image = self['MXMIN_IMAGE'], \
    self['MYMIN_IMAGE'], self['MXMAX_IMAGE'], self['MYMAX_IMAGE']
 
    father_cat_name = self.execpars['father_cat_name']
    father_seg_name = self.execpars['father_seg_name']
    outfile = self.execpars['outfile_name']
    version = self.execpars['version']
    t_str = t.ctime(t.time())
 
    fitsobj = pyfits.HDUList()
    hdu0 = pyfits.PrimaryHDU()
    
    hdu1 = pyfits.ImageHDU(name='STAMP')
    hdu1.data = stamp.copy()
    hdu1.update_header() 
 
    hdu2 = pyfits.ImageHDU(name='SEGSTAMP')
    hdu2.data = segstamp.copy()
    hdu2.update_header()
    
    hdu3 = pyfits.ImageHDU(name='SEGMASK')
    hdu3.data = sexmask.copy()
    hdu3.update_header()
    
    hdu4 = pyfits.ImageHDU(name='MASKOTHER')
    hdu4.data = maskother.copy()
    hdu4.update_header()
    
    fitsobj.append(hdu0) ; fitsobj.append(hdu1)
    fitsobj.append(hdu2) ; fitsobj.append(hdu3)
    fitsobj.append(hdu4)

    fitsobj[0].update_header()
    fitsobj.verify(option = 'warn')
    
    prihdr = fitsobj[0].header
    
    prihdr.update('NAME',quitpath(outname),comment="Name of this file")
    prihdr.update('OUTFILE',quitpath(outfile),comment="Output catalogue")
    prihdr.update('VERSION',version,comment="Version")
    prihdr.update('TIME',t_str,comment="Local Time")
    prihdr.update('FIMAGE',quitpath(father_img_name),comment="Father Image Name")
    prihdr.update('FSEG',quitpath(father_seg_name),comment="Father Segmentation-Map Name")
    prihdr.update('FCAT',quitpath(father_cat_name),comment="Father Catalogue Name")
    prihdr.update('MXMIN',xmin_image,comment="bottom left 'x' coor. in father image")
    prihdr.update('MXMAX',xmax_image,comment="top right 'x' coor. in father image")
    prihdr.update('MYMIN',ymin_image,comment="bottom left 'y' coor. in father image")
    prihdr.update('MYMAX',ymax_image,comment="top right 'y' coor. in father image")
    
    fitsobj.verify(option='fix')
    
    # AVOID OVERWRITING of MEF
    
    isthere = access(outname,F_OK)
    rootname = string.replace(outname,'.fits','')
    
    if isthere == True :
        
        satisfaction = 0			  
        inc = 0 				  
        while satisfaction != 1 :		  
            incstr = str(inc)			    
            testname = rootname + '.'+incstr+'.fits'
            isalsothere = access(testname,F_OK)     
            if isalsothere == True :		    
                inc += 1  			    
            else :				    
                satisfaction = 1  		    
                outname = testname
                prihdr = fitsobj[0].header
                prihdr.update('NAME',outname,comment="Name of this file")
                fitsobj.verify(option='fix')
                fitsobj.writeto(outname)  	    
    else : 
        fitsobj.writeto(outname)
   
    fitsobj.close()
    
    if self.execpars['showds9'][0] == 1:
        myds9 = ds9class()
        isopen = myds9.isOpen()
        if isopen : myds9.xpaset('exit')
        myds9.launch()
        myds9.xpaset('file %s[1]' % outname) ; myds9.xpaset('scale mode 99.5' )
        myds9.xpaset('frame 2')
        myds9.xpaset('file %s[2]' % outname) ; myds9.xpaset('scale mode 99.5' )
        myds9.xpaset('frame 3')
        myds9.xpaset('file %s[3]' % outname) ; myds9.xpaset('scale mode 99.5' )
        myds9.xpaset('frame 1')
        raw_input('Kick any key to continue...')
    
    return outname
  


def read_inputs(infile):
    """Reads inputs for moments.py"""
    # IMPORT STUFF
    from read_sex import read_filesex, read_param
    # END IMPORT STUFF   

    def updatedict(child,parent):
        for key in parent.keys():
           if key not in child : child[key] = parent[key]
        return child
    
    def updatelist(child,parent):
        for item in parent:
           if item not in child : child.append(item)
        return child
    
    # MANDATORY
    exec_mandatory = ['M_RADIAL','M_PETRO','M_PETROMSK','SNR','M_MOM'] # always computed, by default
    sexpars_mandatory = ['XMIN_IMAGE','YMIN_IMAGE','XMAX_IMAGE',\
    'YMAX_IMAGE','NUMBER','THETA_IMAGE','ELLIPTICITY','X_IMAGE','Y_IMAGE']
    to_output_mandatory = ['NUMBER'] # always output, by default
    # END MANDATORY
    
    prim = read_filesex(infile)
    
    sexpars = read_param(prim['sexpars_f'][0])
    toexec = read_param(prim['toexec_f'][0])
    to_output = read_param(prim['to_output_f'][0])
    
    execpars = read_filesex(prim['execpars_f'][0])
    defaults = read_filesex(prim['defaults_f'][0])    
    
    execpars = updatedict(execpars,defaults)
    sexpars = updatelist(sexpars,sexpars_mandatory)
    toexec = updatelist(toexec,exec_mandatory)
    to_output = updatelist(to_output,to_output_mandatory)
    
    exec_dict = {}
    
    exec_dict['father_img_name'] = prim['father_img_name'][0]
    exec_dict['father_cat_name'] = prim['father_cat_name'][0]
    exec_dict['father_sex_name'] = prim['father_sex_name'][0]
    exec_dict['father_seg_name'] = prim['father_seg_name'][0]
    if 'father_mask_name' in prim:  exec_dict ['father_mask_name'] = prim['father_mask_name'][0]
    exec_dict['outfile_name'] = prim['outfile_name'][0]
    exec_dict['sexpars'] = sexpars
    exec_dict['toexec'] = toexec
    exec_dict['to_output'] = to_output
    for key in execpars : exec_dict[key] = execpars[key]
    
    return exec_dict # dictionary with all the inputs for moments.py



def quitpath(filename):
    """Gets complete pointer to file rid of the path. Only file name remains."""
    return filename[filename.rfind('/')+1:]
    
#ruler75###################################################################    
    

def getwindow(self,dimensions=(-1,-1)):
    """Computes a window for the stamps"""
    # iniwindow = (xmin,ymin,xmax,ymax)
    # increase = 50 # %
    # dimensions = (ydim,xdim)
    # IMPORT STUFF
    # END IMPORT
    
    increase = self.execpars['window_inc'][0]
    xmin = self['XMIN_IMAGE']; xmax = self['XMAX_IMAGE']
    ymin = self['YMIN_IMAGE'] ; ymax = self['YMAX_IMAGE']
    iniwindow = (xmin,ymin,xmax,ymax)
    
    if increase == 0 : 
        window = iniwindow
    else:
        
        deltax = ( iniwindow[2] - iniwindow[0] +1) * increase / 100.
        deltay = ( iniwindow[3] - iniwindow[1] +1) * increase / 100.
        delta = num.array([-deltax,-deltay,deltax,deltay])
        window = num.zeros(shape=(4,),dtype='Int32')
        for i in range(4) : window[i] = num.around(iniwindow[i] + delta[i])
        
        window[0] = int(num.clip(window[0],1,window[0]))
        window[1] = int(num.clip(window[1],1,window[1]))
        window[2] = int(num.clip(window[2],iniwindow[2],
        min(window[2],dimensions[1])))
        window[3] = int(num.clip(window[3],iniwindow[3],
        min(window[3],dimensions[0])))
    
    self['MXMIN_IMAGE'] = window[0] ; self['MXMAX_IMAGE'] = window[2]
    self['MYMIN_IMAGE'] = window[1] ; self['MYMAX_IMAGE'] = window[3]

def secs_to_dhms(seconds):
    m,s = divmod(seconds,60.)
    h,m = divmod(m,60.)
    d,h = divmod(h,24.)
    return (d,h,m,s)


def findzeros(x,f):
    """Retrieves zeros of the function f(x)"""
    # IMPORT STUFF
    import scipy.interpolate as interp
    # END IMPORT
    
    test = f[0:-1] * f[1:]
    prezeros = num.where(test < 0.)
    zeros = []
    
    for ix in range(len(prezeros[0])):
        start = prezeros[0][ix]
	xi = f[start:start+2]
	yi = x[start:start+2]
        order = xi.argsort()
	xi = xi[order] ; yi = yi[order]
	inter = interp.interp1d(xi,yi)		
	zero = float(inter(0.))
	zeros.append(zero)
    return num.array(zeros)
