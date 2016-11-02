#! /usr/bin/env python

import numpy as num

# INPUTS
cl_parfor = {'father_cat_name':'Char', 'outfile_name':'Char',
 'execpars_f':'Char', 'defaults_f':'Char',
 'father_img_name':'Char', 'father_seg_name':'Char',
 'father_mask_name':'Char','sexpars_f':'Char',
 'father_sex_name':'Char','toexec_f':'Char', 'to_output_f':'Char',
 'Empty':'Char','eta_petro':'Float32',
 'petrofactor':'Float32','box2petro':'Float32',
 'window_inc':'Int32','usepetro':'Int8',
 'docircularpetro':'Int8','minMuPetro':'Float32','usePetroFirst':'Int8',
 'delay':'Int32','sigma_sky':'Float32','showds9':'Int8',
 'doLogRadial':'Int8','useMANBACK':'Int8','MANBACK':'Float32',
 'saveradial':'Int8','useDo':'Int8','useExPetro':'Int8','cutmode':'Char',
 'dograph':'Int8','Ndump':'Int32','doFineRadial':'Int8',
 'RadialStep':'Float32','doRadSky':'Int8','tag':'Char',
 'clsubradial':'Int8','cldetectother':'Int8','interactive':'Int8',
 'CLTHRESH':'Float32','CLMINAREA':'Float32','magzero':'Float32',
 'scale':'Float32','clmaxradius':'Int8','maskoffinside':'Int32',
 'doMEF':'Int32','clglobal':'Int32','inradial':'Char',
 'doRadialModel':'Int8'}

cl_ordered = num.array(['NUMBER','UNQ_NUMBER','name','X_IMAGE',
 'Y_IMAGE','XPEAK_IMAGE','YPEAK_IMAGE','THETA_IMAGE','ELLIPTICITY',
 'ELONGATION','M_A','A_IMAGE','M_B','B_IMAGE','M_THETA','M_ELONG',
 'M_ELLIP','M_XPEAK','M_YPEAK','EX_R_PETRO','R_PETRO',
 'I_PETRO','F_PETRO','APFLXRADIUS','M_APFLX','M_X','M_Y',
 'M_X2','M_Y2','M_XY','X2_IMAGE','Y2_IMAGE','XY_IMAGE','M_M2','BACKGROUND',
 'M_NPIX','M_FLUX','M_AVINT', 'CL_cat_in','CL_seg_in','Ncl','Fcl',
 'Acl','maxFcl','minFcl','medianFcl', 'meanFcl', 'sigmaFcl',
 'p25Fcl', 'p75Fcl', 'maxAcl', 'minAcl', 'medianAcl', 'meanAcl', 
 'sigmaAcl', 'p25Acl', 'p75Acl','maxRcl', 'minRcl', 'FweightRcl', 'FmaxRcl', 
 'medianRcl','phiAcl','phiFcl','FL_IN','N_IN','F_OUT','N_OUT','F_ALL',
 'N_ALL','MG_IN','AR_IN','MG_OUT','AR_OUT','MG_ALL',
 'AR_ALL','CLTHRESH','CLMINAREA','CL_cat','CL_seg',
 'radial_file','pdf','IMAGE','SEGIMAGE',
 'MASK','flags','DO'])

cl_comm = {}
 
from verbose import sex_regfor
 
cl_regfor = {'M_NPIX':'Int32','M_FLUX':'Float32',
 'M_AVINT':'Float32','M_ELLIP':'Float32','EX_R_PETRO':'Float32',
 'R_PETRO':'Float32','I_PETRO':'Float32','F_PETRO':'Float32',
 'M_A':'Float32','M_B':'Float32','M_THETA':'Float32',
 'M_ELONG':'Float32','M_ELLIP':'Float32',
 'APFLXRADIUS':'Float32','M_APFLX':'Float32', 'flags':'Int64',
 'MXMIN_IMAGE':'Int32','MYMIN_IMAGE':'Int32',
 'MXMAX_IMAGE':'Int32','MYMAX_IMAGE':'Int32',
 'M_X':'Float32','M_Y':'Float32','M_X2':'Float32','M_X2':'Float32',
 'M_Y2':'Float32','M_XY':'Float32','M_M2':'Float32','radial_file':'Char',
 'pdf':'Char','SEGIMAGE':'Char','IMAGE':'Char','MASK':'Char',
 'SKY_SIGMA':'Float32','SKY_MEDIAN':'Float32','CL_cat_in':'Char',
 'CL_seg_in':'Char','Ncl':'Int32','Fcl':'Float32','Acl':'Float32',
 'maxFcl':'Float32','minFcl':'Float32','medianFcl':'Float32',
 'meanFcl':'Float32', 'sigmaFcl':'Float32','p25Fcl':'Float32',
 'p75Fcl':'Float32', 'maxAcl':'Float32','minAcl':'Float32', 
 'medianAcl':'Float32', 'meanAcl':'Float32','sigmaAcl':'Float32', 
 'p25Acl':'Float32', 'p75Acl':'Float32','phiAcl':'Float32',
 'phiFcl':'Float32','CL_cat':'Char','CL_seg':'Char','CL_cat_in':'Char',
 'CL_seg_in':'Char','CLTHRESH':'Float32','CLMINAREA':'Float32',
 'maxRcl':'Float32','minRcl':'Float32', 'FweightRcl':'Float32', 
 'FmaxRcl':'Float32', 'medianRcl':'Float32',
 'MG_IN':'Float32','FL_IN':'Float32','N_IN':'Int32','AR_IN':'Float32',
 'MG_OUT':'Float32','FL_OUT':'Float32','N_OUT':'Int32','AR_OUT':'Float32',
 'MG_ALL':'Float32','FL_ALL':'Float32','N_ALL':'Int32','AR_ALL':'Float32'}

for flag in sex_regfor: cl_regfor[flag] = sex_regfor[flag]

# DATA

Lclsexcols = ['NUMBER','X_IMAGE','Y_IMAGE','A_IMAGE', 
 'B_IMAGE', 'THETA_IMAGE','ELLIPTICITY', 'FLUX_ISO','ISOAREA_IMAGE']

Dclsexformats = {'NUMBER':'Int32','X_IMAGE':'Float32',
 'Y_IMAGE':'Float32','A_IMAGE':'Float32', 'B_IMAGE':'Float32', 
 'THETA_IMAGE':'Float32','ELLIPTICITY':'Float32', 'FLUX_ISO':'Float32',
 'ISOAREA_IMAGE':'Float32'}

Lphotpars = ['CATALOG_TYPE','CATALOG_NAME','PARAMETERS_NAME',
 'DETECT_TYPE','FLAG_TYPE','FLAG_IMAGE','FITS_UNSIGNED',
 'DETECT_MINAREA','DETECT_THRESH','ANALYSIS_THRESH','THRESH_TYPE',
 'FILTER','FILTER_NAME','DEBLEND_NTHRESH','DEBLEND_MINCONT','CLEAN',
 'CLEAN_PARAM','MASK_TYPE','PHOT_APERTURES','PHOT_AUTOPARAMS',
 'PHOT_FLUXFRAC','MAG_ZEROPOINT','GAIN','SATUR_LEVEL','MAG_GAMMA',
 'PIXEL_SCALE','SEEING_FWHM','STARNNW_NAME','BACK_SIZE',
 'BACK_FILTERSIZE','BACKPHOTO_TYPE','BACKPHOTO_THICK','BACK_TYPE',
 'BACK_VALUE','CHECKIMAGE_TYPE','CHECKIMAGE_NAME','MEMORY_OBJSTACK',
 'MEMORY_PIXSTACK','MEMORY_BUFSIZE','VERBOSE_TYPE','INTERP_MAXXLAG',
 'INTERP_MAXYLAG','INTERP_TYPE']

Ddefaultphotpars = {'CATALOG_TYPE':'ASCII',
 'DETECT_TYPE':'CCD','THRESH_TYPE':'ABSOLUTE',
 'FILTER':'N','DEBLEND_NTHRESH':32,
 'DEBLEND_MINCONT':0.005,'CLEAN':'Y',
 'CLEAN_PARAM':1.0,'MASK_TYPE':'NONE','PHOT_APERTURES':5,
 'PHOT_AUTOPARAMS':[2.5, 3.5],'PHOT_FLUXFRAC':0.5,
 'GAIN':1,'BACK_SIZE':[100,100],'BACK_FILTERSIZE':3,
 'BACKPHOTO_TYPE':'GLOBAL','BACK_TYPE':'MANUAL','BACK_VALUE':[0,0],
 'CHECKIMAGE_TYPE':'SEGMENTATION','MEMORY_OBJSTACK':3000,
 'MEMORY_PIXSTACK':300000,'MEMORY_BUFSIZE':1024,'VERBOSE_TYPE':'NORMAL',
 'INTERP_TYPE':'NONE'}
