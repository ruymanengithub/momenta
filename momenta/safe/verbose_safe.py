"""Set of functions for reading Sextractor files."""

def reg_order(outpars):
    """Puts a list of columns in a given order.
    
    outpars = reg_order(outpars)
    
    outpars : a list
    
    """
    # IMPORT STUFF    
    import numpy as num
    from pdb import set_trace as stop
    # END IMPORT
    
    outpars = numstr.array(outpars)
    ordered = numstr.array(['NUMBER','X_IMAGE','Y_IMAGE','XPEAK_IMAGE','YPEAK_IMAGE',\
    'THETA_IMAGE','ELLIPTICITY','ELONGATION','M_A','A_IMAGE','M_B','B_IMAGE','M_THETA',\
    'M_ELONG','M_ELLIP','M_XPEAK','M_YPEAK','M_GINI','R_PETRO','I_PETRO',\
    'F_PETRO','M_AS','M_AS_X','M_AS_Y','M_AS_SKY','M_C','M_S','M_S_SKY','M_AXS_MAJ',\
    'M_AXS_MIN','M_AXS_MAJ_SKY','M_AXS_MIN_SKY','M_X','M_Y','M_X2','M_Y2','M_XY',\
    'X2_IMAGE','Y2_IMAGE','XY_IMAGE','BACKGROUND','MAG_BEST','FLUX_BEST','SNR','flags'])

    order = num.array([num.where(item == ordered)[0][0] for item in outpars])
    order = order.argsort()
    newoutpars = outpars[order]
    
    return newoutpars

def reg_format(flags):
    """Asigns formats to each column in a Sextractor/moments.py catalog
    Needs improvements: Just updating.
    """

    #IMPORT STUFF
    import string
    #END IMPORT STUFF

    # print "UPDATE THE DICTIONARY OF FORMATS! YOU LAZY BUM..."

    allflags = {'M_GINI':'Float32','M_ELLIP':'Float32','R_PETRO':'Float32','I_PETRO':'Float32',\
    'F_PETRO':'Float32','M_A':'Float32','M_B':'Float32','M_THETA':'Float32','M_ELONG':'Float32',\
    'M_ELLIP':'Float32','M_XPEAK':'Float32','M_YPEAK':'Float32','M_AS':'Float32',\
    'M_AS_SKY':'Float32','M_AS_X':'Float32','M_AS_Y':'Float32','M_C':'Float32','M_S':'Float32',\
    'M_S_SKY':'Float32','M_AXS_MAJ':'Float32','M_AXS_MIN':'Float32','M_AXS_MAJ_SKY':'Float32',\
    'M_AXS_MIN_SKY':'Float32','flags':'Int64','MXMIN_IMAGE':'Int32','MYMIN_IMAGE':'Int32',\
    'MXMAX_IMAGE':'Int32','MYMAX_IMAGE':'Int32','M_X':'Float32','M_Y':'Float32','M_X2':'Float32',\
    'M_X2':'Float32','M_Y2':'Float32','M_XY':'Float32','SNR':'Float32'}

    sexflags = {'FLUX_AUTO':'Float32', 'FLUXERR_AUTO':'Float32','MAG_AUTO':'Float32',\
    'MAGERR_AUTO':'Float32', 'KRON_RADIUS':'Float32', 'BACKGROUND':'Float32', \
    'FLUX_RADIUS' :'Float32','SOMETHING1':'Float32', 'SOMETHING2':'Float32', \
    'THRESHOLD':'Float32', 'MU_THRESHOLD':'Float32', 'FLUX_MAX':'Float32', \
    'MU_MAX':'Float32', 'ISOAREA_IMAGE':'Float32', 'ISOAREA_WORLD':'Float32', \
    'XMIN_IMAGE':'Int32', 'YMIN_IMAGE':'Int32', 'XMAX_IMAGE':'Int32', \
    'YMAX_IMAGE':'Int32', 'X_IMAGE':'Float32','Y_IMAGE':'Float32', 'X_WORLD':'Float32', \
    'Y_WORLD':'Float32', 'A_IMAGE':'Float32', 'B_IMAGE':'Float32', 'A_WORLD':'Float32', \
    'B_WORLD':'Float32','THETA_IMAGE':'Float32', 'THETA_WORLD':'Float32', \
    'ELONGATION':'Float32', 'ELLIPTICITY':'Float32','ERRTHETA_IMAGE':'Float32', \
    'ERRTHETA_WORLD':'Float32', 'FWHM_IMAGE':'Float32', 'FWHM_WORLD':'Float32', \
    'FLAGS':'Int32','CLASS_STAR':'Float32','NUMBER':'Int32','UNKNOWN':'Float32',\
    'VECTOR_ASSOC':'Float32','ALPHA_J2000':'Float32','DELTA_J2000':'Float32',\
    'XPEAK_IMAGE':'Float32','YPEAK_IMAGE':'Float32','YPEAK_WORLD':'Float32',\
    'XPEAK_WORLD':'Float32','ISOAREAF_IMAGE':'Float32','CLASS_STAR':'Float32',\
    'IMAFLAGS_ISO':'Float32', 'NIMAFLAGS_ISO':'Int32', 'MAG_ISO':'Float32',\
    'MAGERR_ISO':'Float32','FLUX_ISO':'Float32','FLUXERR_ISO':'Float32','MAG_ISOCOR':\
    'Float32','MAGERR_ISOCOR':'Float32','FLUX_ISOCOR':'Float32','FLUXERR_ISOCOR':'Float32',\
    'MAG_BEST':'Float32','MAGERR_BEST':'Float32','FLUX_BEST':'Float32','FLUXERR_BEST':\
    'Float32','MAG_APER':'Float32','MAGERR_APER':'Float32','FLUX_APER':'Float32',\
    'FLUXERR_APER':'Float32','X2_IMAGE':'Float32','Y2_IMAGE':'Float32','XY_IMAGE':\
    'Float32','ERRX2_IMAGE':'Float32','ERRY2_IMAGE':'Float32','ERRXY_IMAGE':'Float32',\
    'ERRA_IMAGE':'Float32','ERRB_IMAGE':'Float32'}

    for flag in sexflags: allflags[flag] = sexflags[flag] 
    
    formats = {}							     
    for flag in flags: 						     
        check = string.rfind(flag,'.') 				     
        if check >= 0 : 						     
            formats[flag] = allflags[flag[0:check]] # ignores .n extensions!  
        else : formats[flag] = allflags[flag]				     

    return formats  						     


def reg_comment(flags):
    """Asigns comment to each parameter in the header of a in a moments.py catalog

    Needs improvements: Just updating.
    """
    #IMPORT STUFF
    import string
    #END IMPORT STUFF

    allflags = {'M_GINI':'(My) Gini Parameter                                  [adim]',\
    'F_PETRO':'(My) Flux inside Petrosian radius                    [count]',\
    'R_PETRO':'(My) Petrosian Radius                                [pixel]',\
    'I_PETRO':'(My) Intensity at Petrosian radius                   [count]',\
    'M_A':'(My) A (geom.)                                       [adim]',\
    'M_B':'(My) B (geom.)                                       [adim]',\
    'M_THETA':'(My) Polar angle to NAXIS1 (ccw)                     [deg]',\
    'M_ELONG':'(My) Elongation                                      [adim]',\
    'M_ELLIP':'(My) Ellipticity                                     [adim]',\
    'M_XPEAK':'(My) X COOR. PEAK                                    [pixel]',\
    'M_YPEAK':'(My) Y COOR. PEAK                                    [pixel]',\
    'M_AS':'(My) Asymmetry                                       [adim]',\
    'M_AS_SKY':'(My) Asymmetry of the sky            [adim]',\
    'M_AS_X':'(My) Asymmetry center, x                             [pixel]',\
    'M_AS_Y':'(My) Asymmetry center y                              [pixel]',\
    'M_C':'(My) Concentration                                   [adim]',\
    'M_S':'(My) Clumpyness                                      [adim]',\
    'M_S_SKY':'(My) Clumpyness for sky                [adim]',\
    'M_AXS_MAJ':'(My) Axis Asymmetry (major axis)                      [adim]',\
    'M_AXS_MIN':'(My) Axis Asymmetry (minor axis)                      [adim]',\
    'M_AXS_MAJ_SKY':'(My) Axis Asymmetry for sky (major axis)   [adim]',\
    'M_AXS_MIN_SKY':'(My) Axis Asymetry for sky (minor axis)     [adim]',\
    'flags':'(My) flags                                           [adim]',\
    'M_X':'(My) First order X moment                            [pixel]',\
    'M_Y':'(My) First order Y moment                            [pixel]',\
    'M_X2':'(My) Second order X moment                           [pixel**2]',\
    'M_Y2':'(My) Second order Y moment                           [pixel**2]',\
    'M_XY':'(My) Second order XY moment                          [pixel**2]',\
    'SNR':'(My) Average SNR per pixel                                 [adim]'}
    
    sexflags = {'FLUX_AUTO':'Flux within a Kron-like elliptical aperture     [mag]',\
   'NUMBER':'',\
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
   'FLUX_BEST':''}

    for key in sexflags : allflags[key] = sexflags[key]
   
    comments = {}
    for flag in flags: 
        check = string.rfind(flag,'.')
        if check >= 0 : 
            comments[flag] = allflags[flag[0:check]] # ignores .n extensions!
        else : comments[flag] = allflags[flag]

    return comments
def par_format(flags) :
    """Asigns formats to each entry in a Sextractor/moments.py configuration file
    (like .sex files)
    
    Needs improvements: Just updating.
    """

    #IMPORT STUFF
    import string
    from pdb import set_trace as stop
    #END IMPORT STUFF

    # print "UPDATE THE DICTIONARY OF FORMATS! YOU LAZY BUM..."

    allflags = {'father_cat_name':'Char', 'outfile_name':'Char', 'execpars_f':'Char', 'defaults_f':'Char', \
    'father_img_name':'Char', 'father_seg_name':'Char', 'father_mask_name':'Char','sexpars_f':'Char', \
    'father_sex_name':'Char','toexec_f':'Char', 'to_output_f':'Char','Empty':'Char','makefits':'Int8',\
    'eta_petro':'Float32','petrofactor':'Float32','box2petro':'Float32','r_big_proc':'Float32','r_small_proc':\
    'Float32','window_inc':'Int32','usepetro':'Int8','docircularpetro':'Int8','Asym_iterative':'Int8','delay':'Int32',\
    'sigma_sky':'Float32','doAsym_sky':'Int8','doAxisAsym_sky':'Int8','doClumpy_sky':'Int8','showds9':'Int8',\
    'doAxisAsym_iter':'Int8'}
    
    sexflags = {'CHECKIMAGE_NAME':'Char', 'PARAMETERS_NAME':'Char', 'ANALYSIS_THRESH':'Float32',\
    'THRESH_TYPE':'Char', 'PHOT_AUTOPARAMS':'Float32', 'INTERP_TYPE':'Char', 'WEIGHT_THRESH':'Int16',\
    'GAIN':'Float32', 'MEMORY_BUFSIZE':'Int32', 'FILTER_NAME':'Char', 'PHOT_AUTOAPERS':'Float32', \
    'CLEAN_PARAM':'Float32', 'BACK_FILTERSIZE':'Int16', 'FILTER':'Char', 'PHOT_FLUXFRAC':'Float32',\
    'DEBLEND_MINCONT':'Float32', 'CATALOG_NAME':'Char','FLAG_IMAGE':'Char', 'DETECT_THRESH':\
    'Float32', 'CATALOG_TYPE':'Char', 'SATUR_LEVEL':'Float32', 'MASK_TYPE':'Char','STARNNW_NAME':\
    'Char', 'MEMORY_PIXSTACK':'Int32', 'DETECT_TYPE':'Char','INTERP_MAXYLAG':'Int16', 'WEIGHT_TYPE':\
    'Char', 'MEMORY_OBJSTACK':'Int32', 'PIXEL_SCALE':'Float32', 'MAG_GAMMA':'Float32',\
    'INTERP_MAXXLAG':'Int16', 'CHECKIMAGE_TYPE':'Char', 'CLEAN':'Char', 'WEIGHT_IMAGE':'Char',\
    'MAG_ZEROPOINT':'Float32', 'PHOT_APERTURES':'Float32', 'FITS_UNSIGNED':'Char', 'FLAG_TYPE':\
    'Char', 'VERBOSE_TYPE':'Char', 'DEBLEND_NTHRESH':'Int16','DETECT_MINAREA':'Int16', \
    'BACKPHOTO_THICK':'Int16', 'BACKPHOTO_TYPE':'Char', 'SEEING_FWHM':'Float32', 'BACK_SIZE':\
    'Int16'}
    
    for flag in sexflags : allflags[flag] = sexflags[flag]
    
    formats = {}
    for flag in flags: 
        check = string.rfind(flag,'.')
        if check >= 0 : 
            formats[flag] = allflags[flag[0:check]] # ignores .n extensions!
        else : 
           try : formats[flag] = allflags[flag]
           except KeyError: stop()
    return formats 



