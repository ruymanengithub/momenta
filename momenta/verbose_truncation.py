#! /usr/bin/env python

import numpy as num

mom_ordered = num.array(['NUMBER','UNQ_NUMBER','name','X_IMAGE',\
 'Y_IMAGE','XPEAK_IMAGE','YPEAK_IMAGE','MXMIN_IMAGE','MYMIN_IMAGE',
 'MXMAX_IMAGE','MYMAX_IMAGE','THETA_IMAGE','ELLIPTICITY',\
 'ELONGATION','M_A','A_IMAGE','M_B','B_IMAGE','M_THETA','M_ELONG',
 'M_ELLIP','M_X','M_Y','M_X2','M_Y2','M_M2','M_NPIX','M_FLUX','M_AVINT',
 'TType','Tx','Ty','Tback1','Tback2','Tskyradius','Tb1x','Tb1y',
 'Tb2x','Tb2y','Tb3x','Tb3y','Tb4x','Tb4y','Tm1','eTm1','Tc1',
 'eTc1','Trpearson1',
 'Th1','eTh1','Tm2','eTm2','Tc2','eTc2','Trpearson2','Th2','eTh2','Tr',
 'eTr','Tmufit','eTmufit','Tmu','Tzero','Tscale','TBelieve','Tradialf',
 'Tpdf','BACKGROUND','MAG_BEST','FLUX_BEST','IMAGE','SEGIMAGE','MASK',
 'flags','DO'])

#ruler75############################################

mom_regfor = {'M_ELLIP':'Float32','M_A':'Float32','M_B':'Float32',
 'M_THETA':'Float32','M_ELONG':'Float32','M_ELLIP':'Float32',
 'M_M2':'Float32','M_X':'Float32','M_Y':'Float32','M_X2':'Float32',
 'M_Y2':'Float32','M_XY':'Float32','TType':'Char','Tx':'Float32',
 'Ty':'Float32','Tback1':'Float32',
 'Tskyradius':'Float32','Tback2':'Float32','Tb1x':'Float32',
 'Tb1y':'Float32','Tb2x':'Float32','Tb2y':'Float32','Tb3x':'Float32',
 'Tb3y':'Float32','Tb4x':'Float32','Tb4y':'Float32','Tm1':'Float32',
 'eTm1':'Float32','Tc1':'Float32','eTc1':'Float32','Trpearson1':
 'Float32','Th1':'Float32','eTh1':'Float32','Tm2':'Float32','eTm2':
 'Float32','Tc2':'Float32','eTc2':'Float32','Trpearson2':'Float32',
 'Th2':'Float32','eTh2':'Float32','Tr':'Float32','eTr':'Float32',
 'Tmufit':'Float32','eTmufit':'Float32','Tmu':'Float32',
 'Tzero':'Float32','Tscale':'Float32','TBelieve':'Char',
 'Tradialf':'Char','Tpdf':'Char','SEGIMAGE':'Char','IMAGE':'Char',
 'MASK':'Char','M_NPIX':'Int32','M_FLUX':'Float32','M_AVINT':'Float32',\
 'flags':'Int32','MXMIN_IMAGE':'Int32','MYMIN_IMAGE':'Int32',
 'MXMAX_IMAGE':'Int32','MYMAX_IMAGE':'Int32',}

sex_regfor = {'UNQ_NUMBER':'Int32','FLUX_AUTO':'Float32',\
 'FLUXERR_AUTO':'Float32','MAG_AUTO':'Float32',\
 'MAGERR_AUTO':'Float32', 'KRON_RADIUS':'Float32',\
 'BACKGROUND':'Float32', 'FLUX_RADIUS' :'Float32',\
 'SOMETHING1':'Float32', 'SOMETHING2':'Float32', \
 'THRESHOLD':'Float32', 'MU_THRESHOLD':'Float32',\
 'FLUX_MAX':'Float32','MU_MAX':'Float32', 'ISOAREA_IMAGE':'Float32',\
 'ISOAREA_WORLD':'Float32', 'XMIN_IMAGE':'Int32', \
 'YMIN_IMAGE':'Int32', 'XMAX_IMAGE':'Int32', 'YMAX_IMAGE':'Int32', \
 'X_IMAGE':'Float32','Y_IMAGE':'Float32', 'X_WORLD':'Float32', \
 'Y_WORLD':'Float32', 'A_IMAGE':'Float32', 'B_IMAGE':'Float32',\
 'A_WORLD':'Float32','B_WORLD':'Float32','THETA_IMAGE':'Float32',\
 'THETA_WORLD':'Float32','ELONGATION':'Float32', \
 'ELLIPTICITY':'Float32','ERRTHETA_IMAGE':'Float32',\
 'ERRTHETA_WORLD':'Float32', 'FWHM_IMAGE':'Float32', \
 'FWHM_WORLD':'Float32','FLAGS':'Int32','CLASS_STAR':'Float32',\
 'NUMBER':'Int32','UNKNOWN':'Float32','VECTOR_ASSOC':'Float32',\
 'ALPHA_J2000':'Float32','DELTA_J2000':'Float32',\
 'XPEAK_IMAGE':'Float32','YPEAK_IMAGE':'Float32',\
 'YPEAK_WORLD':'Float32','XPEAK_WORLD':'Float32',\
 'ISOAREAF_IMAGE':'Float32','CLASS_STAR':'Float32',\
 'IMAFLAGS_ISO':'Float32', 'NIMAFLAGS_ISO':'Int32',\
 'MAG_ISO':'Float32','MAGERR_ISO':'Float32',\
 'FLUX_ISO':'Float32','FLUXERR_ISO':'Float32',\
 'MAG_ISOCOR':'Float32','MAGERR_ISOCOR':'Float32',\
 'FLUX_ISOCOR':'Float32','FLUXERR_ISOCOR':'Float32',\
 'MAG_BEST':'Float32','MAGERR_BEST':'Float32',\
 'FLUX_BEST':'Float32','FLUXERR_BEST':'Float32',\
 'MAG_APER':'Float32','MAGERR_APER':'Float32',\
 'FLUX_APER':'Float32','FLUXERR_APER':'Float32',\
 'X2_IMAGE':'Float32','Y2_IMAGE':'Float32',
 'XY_IMAGE':'Float32','ERRX2_IMAGE':'Float32',\
 'ERRY2_IMAGE':'Float32','ERRXY_IMAGE':'Float32',\
 'ERRA_IMAGE':'Float32','ERRB_IMAGE':'Float32','name':'Char',\
 'DO':'Int8','EX_R_PETRO':'Float32','SEGIMAGE':'Char','IMAGE':'Char',\
 'MASK':'Char','PDF':'Char'}

for flag in sex_regfor: mom_regfor[flag] = sex_regfor[flag] 
    
#ruler75##############################################

mom_comm = {
 'M_A':'(My) A (geom.)                                       [adim]',\
 'M_B':'(My) B (geom.)                                       [adim]',\
 'M_THETA':'(My) Polar angle to NAXIS1 (ccw)                     [deg]',\
 'M_ELONG':'(My) Elongation                                       [adim]',\
 'M_ELLIP':'(My) Ellipticity                                     [adim]',\
 'M_X':'(My) First order X moment                            [pixel]',
 'M_Y':'(My) First order Y moment                            [pixel]',
 'M_X2':'(My) Second order X moment                           [pixel**2]',
 'M_Y2':'(My) Second order Y moment                           [pixel**2]',
 'M_XY':'(My) Second order XY moment                          [pixel**2]',
 'flags':'(My) flags                                           [adim]',\
 'Tradialf':'File with radial profile.',
 'Tpdf':'PDF file with graphs and results.',
 'TType':'Type of Truncation.',
 'Tx':'Recentered X',
 'Ty':'Recentered Y',
 'Tback1':'1st Background substracted (from image)',
 'Tskyradius':'Radius to measure 2nd Background (arcsecs)',
 'Tback2':'2nd Background substracted (from profile)',
 'Tb1x':'left edge 1st region, x',
 'Tb1y':'left edge 1st region, y',
 'Tb2x':'right edge 1st region, x',
 'Tb2y':'right edge 1st region, y',
 'Tb3x':'left edge 2nd region, x',
 'Tb3y':'left edge 2nd region, y',
 'Tb4x':'right edge 2nd region x',
 'Tb4y':'right edge 2nd region, y',
 'Tm1':'slope fit 1st region',
 'eTm1':'error in Tm1',
 'Tc1':'cut with y axis fit 1st region',
 'eTc1':'error in Tc1',
 'Trpearson1':'Pearson r coeff fit 1',
 'Th1':'scale fit 1 (arcsecs)',
 'eTh1':'error in Th1',
 'Tm2':'slope fit 2nd region',
 'eTm2':'error in Tm2',
 'Tc2':'cut with y axis fit 2nd region',
 'eTc2':'error in Tc2',
 'Trpearson2':'Pearson r coeff fit 2',
 'Th2':'scale fit 2 (arcsecs)',
 'eTh2':'erro in Th2',
 'Tr':'Truncation Radius (arcsecs)',
 'eTr':'error in Tr',
 'Tmufit':'mu of intersections (mag/arcsec^2)',
 'eTmufit':'error in Tmufit',
 'Tmu':'Intensity at Tr (mag/arcsec^2)',
 'Tzero':'zeropoint',
 'Tscale':'angular scale (arcsec/pix)',
 'TBelieve':'Trusted Result?',
 'M_NPIX':'Number of pixels inside MASK',
 'M_FLUX':'Flux in counts inside MASK',
 'M_AVINT':'Flux per pixel inside MASK'}

sex_comm = {
 'FLUX_AUTO':'Flux within a Kron-like elliptical aperture     [mag]',\
 'NUMBER':'',\
 'UNQ_NUMBER':'',\
 'FLUXERR_AUTO':'RMS error for AUTO flux			   [mag]',\
 'MAG_AUTO':'Kron-like elliptical aperture magnitude         [mag]',\
 'MAGERR_AUTO':'RMS error for AUTO magnitude  		  [mag]',\
 'KRON_RADIUS':'Kron apertures in units of A or B',\
 'BACKGROUND':'Background at centroid position		 [count]',\
 'FLUX_RADIUS' :'Fraction-of-light radii			   [pixel]',\
 'SOMETHING1':'[pixel]',\
 'SOMETHING2':'[pixel]',\
 'THRESHOLD':'Detection threshold above background		[count]',\
 'MU_THRESHOLD':'Detection threshold above background 	   [mag * arcsec**(-2)]',\
 'FLUX_MAX':'Peak flux above background		       [count]',\
 'MU_MAX':'Peak surface brightness above background	     [mag * arcsec**(-2)]',\
 'ISOAREA_IMAGE':'Isophotal area above Analysis threshold	    [pixel**2]',\
 'ISOAREA_WORLD':'Isophotal area above Analysis threshold[deg**2]',\
 'XMIN_IMAGE':'Minimum x-coordinate among detected pixels	 [pixel]',\
 'YMIN_IMAGE':'Minimum y-coordinate among detected pixels	 [pixel]',\
 'XMAX_IMAGE':'Maximum x-coordinate among detected pixels	 [pixel]', \
 'YMAX_IMAGE':'Maximum y-coordinate among detected pixels	 [pixel]',\
 'X_IMAGE':'Object position along x			           [pixel]',\
 'Y_IMAGE':'Object position along y			           [pixel]',\
 'X_WORLD':'Barycenter position along world x axis	      [deg]',\
 'Y_WORLD':'Barycenter position along world y axis	      [deg]',\
 'A_IMAGE':'Profile RMS along major axis		      [pixel]',\
 'B_IMAGE':'Profile RMS along minor axis		      [pixel]',\
 'A_WORLD':'Profile RMS along major axis (world units)      Float32', \
 'B_WORLD':'Profile RMS along minor axis (world units)      Float32',\
 'THETA_IMAGE':'Position angle (CCW/x)			           [Deg]',\
 'THETA_WORLD':'Position angle (CCW/world-x)  		       [Deg]  ',\
 'ELONGATION':'A_IMAGE/B_IMAGE                                    [adim]',\
 'ELLIPTICITY':'1 - B_IMAGE/A_IMAGE                                  [adim]',\
 'ERRTHETA_IMAGE':'Error ellipse position angle (CCW/x)	     [deg]', \
 'ERRTHETA_WORLD':'Error ellipse pos. angle (CCW/world-x)	     [deg]',\
 'FWHM_IMAGE':'FWHM assuming a gaussian core  		 [pixel]',\
 'FWHM_WORLD':'FWHM assuming a gaussian core  		 [deg]  ', \
 'FLAGS':'Extraction flags',\
 'CLASS_STAR':'S/G classifier output',\
 'NUMBER':'',\
 'UNKNOWN':'',\
 'VECTOR_ASSOC':'',\
 'XPEAK_IMAGE':'',\
 'YPEAK_IMAGE':'',\
 'X2_IMAGE':'',\
 'Y2_IMAGE':'',
 'XY_IMAGE':'',\
 'MAG_BEST':'',\
 'FLUX_BEST':'',\
 'name':'',\
 'DO':'',\
 'PDF':''}

for key in sex_comm : mom_comm[key] = sex_comm[key]
    
#ruler75#############################################

mom_parfor = {'father_cat_name':'Char', 'outfile_name':'Char',
 'execpars_f':'Char', 'defaults_f':'Char',
 'father_img_name':'Char', 'father_seg_name':'Char',
 'father_mask_name':'Char','sexpars_f':'Char',
 'father_sex_name':'Char','toexec_f':'Char', 'to_output_f':'Char',
 'makefits':'Int8','showds9':'Int8','window_inc':'Int32','delay':'Int32',
 'sigma_sky':'Float32','doLogRadial':'Int8','useMANBACK':'Int8',
 'MANBACK':'Float32','saveradial':'Int8',
 'useDo':'Int8','cutmode':'Char','Ndump':'Int32','doFineRadial':'Int8',
 'RadialStep':'Float32','magzero':'Float32','scale':'Float32',
 'tag':'Char','dograph':'Int8','commfile':'Char',
 'manual':'Int8'}

sex_parfor = {'CHECKIMAGE_NAME':'Char', 'PARAMETERS_NAME':'Char',   
 'ANALYSIS_THRESH':'Float32','THRESH_TYPE':'Char',
 'PHOT_AUTOPARAMS':'Float32', 'INTERP_TYPE':'Char', 
 'WEIGHT_THRESH':'Int16','GAIN':'Float32',
 'MEMORY_BUFSIZE':'Int32', 'FILTER_NAME':'Char',
 'PHOT_AUTOAPERS':'Float32', 'CLEAN_PARAM':'Float32', 
 'BACK_FILTERSIZE':'Int16', 'FILTER':'Char', 
 'PHOT_FLUXFRAC':'Float32','DEBLEND_MINCONT':'Float32', 
 'CATALOG_NAME':'Char','FLAG_IMAGE':'Char','DETECT_THRESH':'Float32', 
 'CATALOG_TYPE':'Char', 'SATUR_LEVEL':'Float32', 
 'MASK_TYPE':'Char','STARNNW_NAME':'Char', 
 'MEMORY_PIXSTACK':'Int32', 'DETECT_TYPE':'Char',
 'INTERP_MAXYLAG':'Int16', 'WEIGHT_TYPE':'Char', 
 'MEMORY_OBJSTACK':'Int32', 'PIXEL_SCALE':'Float32',
 'MAG_GAMMA':'Float32','INTERP_MAXXLAG':'Int16', 
 'CHECKIMAGE_TYPE':'Char', 'CLEAN':'Char','WEIGHT_IMAGE':'Char',
 'MAG_ZEROPOINT':'Float32', 'PHOT_APERTURES':'Float32',
 'FITS_UNSIGNED':'Char', 'FLAG_TYPE':'Char', 
 'VERBOSE_TYPE':'Char', 'DEBLEND_NTHRESH':'Int16',
 'DETECT_MINAREA':'Int16',
 'BACKPHOTO_THICK':'Int16', 'BACKPHOTO_TYPE':'Char',
 'SEEING_FWHM':'Float32', 'BACK_SIZE':'Int16',
 'BACK_VALUE':'Float32','BACK_TYPE':'Char'}

for flag in sex_parfor : mom_parfor[flag] = sex_parfor[flag]

#ruler75##############################################
