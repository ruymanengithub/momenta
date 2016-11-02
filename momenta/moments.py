#! /usr/bin/env python

import numpy as num
from pdb import set_trace as stop

def moments(infile):
    """Program for image analysis, based on Sextractor photometry.
        It is aimed to compute several parameters of the spatial distribution
        of light in galaxies, not available in sextractor."""
    
    version = 0.0
    #ruler75###########################################
    # ALL THE IMPORT STUFF 
    import read_sex
    import string
    from flags import addflag, allflags,isflagon
    from momsource import momsource, momcat
    from time import time
    import os, sys
    import pyfits
    from CommonTools.loadascii import loadascii
    import momlib
    quitpath = momlib.quitpath
    secs_to_dhms = momlib.secs_to_dhms
    from Moments.latex import LaTeX
    from copy import copy
    from algorithms import get_stat
    # END IMPORT STUFF
    
    isthere = os.path.exists
    
    t1 = time()
    
    execpars = momlib.read_inputs(infile) # execution parameters

    #ruler75############################################
    
    # READOUT OF CATALOGUE
    
    father_cat = momcat(execpars) # initialize input catalog
    father_cat.read_catsex()             # read input catalog
    
    for key in father_cat.flags:
        if father_cat.flags[key] == False :
            print """SOME REQUESTED PARAMETER NOT IN CATALOGUE: %s""" %\
            key
            stop()
    del key
    
    # READOUT OF SEX FILE
    # dictionary with keywords:values
    # father_sex = read_sex.read_filesex(execpars['father_sex_name'])
    # Not used by now.
    
    # Filtering of catalog using an external file which selects which objects
    # in the SEXTRACTOR catalog are to be analyzed.
    
    if father_cat.execpars['NUMBERs_f'][0] != 'None' :
        NUMBERs_f = father_cat.execpars['NUMBERs_f'][0]
        cols = ['NUMBER']
        formats = 'i'
        load = loadascii(NUMBERs_f,cols,formats,separator='blank')
        NUMBERs = load['NUMBER']
        indxs = []
        for number in NUMBERs:
            indxs.append(num.where(cat['NUMBER'] == number)[0])
        indxs = (num.array(indxs),)
        for key in father_cat:
            father_cat[key] = father_cat[key][indxs]
        
    elif 'DO' in father_cat and father_cat.execpars['useDo'][0] == 1:
        indxs = num.where(father_cat['DO'] == 1)
        NUMBERs = father_cat['NUMBER'][indxs]
        for key in father_cat:
            father_cat[key] = father_cat[key][indxs]
    
    if father_cat.execpars['delay'][0] != 0:
        delay = father_cat.execpars['delay'][0]
        for key in father_cat:
            father_cat[key] = father_cat[key][delay:]
    
    nobj_out = len(father_cat[father_cat.keys()[0]]) # number of objects to 
    # analyze
    
    father_cat.initialize() # initialize output catalog
    
    # LOADING IMAGES:
    
    extension = 0 # Image extension. Not ready to handle higer extensions.
    
    # ARE IMAGES TOO BIG? THEN USE FITSCUT TO MAKE STAMPS.
    stampmode = execpars['cutmode'][0]
    
    if 'SEGIMAGE' not in father_cat and 'IMAGE' not in father_cat:
    
        if stampmode == 'pyfits':
            father_img = pyfits.getdata(execpars['father_img_name'],\
            ext=extension).astype('Float32')
            father_seg = pyfits.getdata(execpars['father_seg_name'],\
            ext=extension).astype('Int32')
            try: father_mask = pyfits.getdata(execpars['father_mask_name'],\
            ext = extension).astype('Int32')
            except KeyError : pass
            father_img_dim = father_img.shape
        elif stampmode == 'fitscut' or stampmode == 'pyraf':
            father_img_name = execpars['father_img_name']
            father_seg_name  = execpars['father_seg_name']
            try : 
                father_mask_name = execpars['father_mask_name']
            except KeyError:
                pass
            father_img_dim = pyfits.getdata(father_img_name,\
            ext=extension).shape
    
    # LOOP OVER SOURCES AND SOURCE OBJECT CREATION
    
    father_img_name_prev = ''
    
    try: Ndump = father_cat.execpars['Ndump'][0]
    except KeyError: Ndump = 0
    
    # MAIN LOOP
    
    for counter in range(nobj_out):
        
        if 'SEGIMAGE' in father_cat and 'IMAGE' in father_cat:
            
            if father_cat['IMAGE'][counter] != father_img_name_prev:
                father_img_name = father_cat['IMAGE'][counter]
                father_seg_name = father_cat['SEGIMAGE'][counter]
                try: father_mask_name = father_cat['MASK'][counter]
                except KeyError: father_mask_name = 'None'
                
                if stampmode == 'pyfits':
                    father_img = pyfits.getdata(father_img_name,\
                    ext=extension).astype('Float32')
                    father_img_dim = father_img.shape
                    father_seg = pyfits.getdata(father_seg_name,\
                    ext=extension).astype('Int32')
                    if isthere(father_mask_name):
                        father_mask = pyfits.getdata(father_mask_name,\
                        ext=extension).astype('Int32')
                    else : pass
                    
                else:
                    father_img_dim = pyfits.getdata(father_img_name,\
                    ext= extension).shape
 
                father_img_name_prev = father_img_name
                
            else : pass
        
        # TIME CONTROL
        
        t3 = time()
        print '\nprocessing object %i of a total of %i\n' % \
        (counter+1,nobj_out)
        talready = t3-t1
        talready_fr = secs_to_dhms(talready)
        if talready_fr[0]>1:
            print '... %i days, %i hours, %i minutes, %f seconds since start' % \
            talready_fr[0:]
        else:
            print '...%i hours, %i minutes, %f seconds since start' % \
            talready_fr[1:]
        
        tahead = ((t3-t1)/(counter+1))* (nobj_out-counter+1)
        tahead_fr = secs_to_dhms(tahead)
        if tahead_fr[0] > 0 :
            print '%i days, %i hours, %.2f minutes to finish...' % \
            tahead_fr[0:3]
        else:
            print '%i hours, %.2f minutes to finish...' % tahead_fr[1:3]
        
        source = momsource(execpars) # inherits from dict class.
        for col in father_cat : source[col] = father_cat[col][counter]		      
        del col
	
        if source.execpars['useMANBACK'][0] == 1:
            source['BACKGROUND'] = source.execpars['MANBACK'][0]
        
        # WINDOW DEFINITION
        
        source.getwindow(father_img_dim)
        
        # MAKING STAMPS
        extension = 0 # default by now
        
	
        if stampmode == 'pyfits':
            source.make_stamp(imgname=None,img=father_img,\
            name='STAMP',extension=extension,mode=stampmode)
            source.make_stamp(imgname=None,img=father_seg,\
            name='SEGSTAMP',extension=extension,mode=stampmode)
            try: source.make_stamp(imgname=None,img=father_mask,\
            name='EXTMASK',extension=extension,mode=stampmode)
            except NameError : source['EXTMASK'] = None
            
        elif stampmode == 'fitscut' or stampmode=='pyraf':
            
            source.make_stamp(imgname=father_img_name,img=None,\
            name='STAMP',extension=extension,mode=stampmode)
            source.make_stamp(imgname=father_seg_name,img=None,\
            name='SEGSTAMP',extension=extension,mode=stampmode)
            if father_mask_name != 'None':
                source.make_stamp(imgname=father_mask_name,img=None,\
                name='EXTMASK',extension=extension,mode=stampmode)
            else : source['EXTMASK'] = None
        
        # "MASQUERADE..."
        
        # 1 is masked, 0 is non masked
        
        source.make_mask(source['SEGSTAMP'],name='SEXMASK',\
        mask_in=None,mode='nosky') # array.
        source.make_mask(source['SEGSTAMP'],name='SEXMASKOTHER',\
        mask_in=None,mode='withsky') # array.


        if source['EXTMASK'] != None : 
            source['MASK'] = momlib.mergemasks((source['SEXMASK'],\
                source['EXTMASK']))
            source['MASKOTHER'] = momlib.mergemasks((source['SEXMASKOTHER'],\
                source['EXTMASK']))
            source['SKYMASK'] = momlib.mergemasks((1-source['SEXMASK'],\
                source['SEXMASKOTHER'],source['EXTMASK']))
            
        else : 
            source['MASK'] = source['SEXMASK']
            source['MASKOTHER'] = source['SEXMASKOTHER']
            source['SKYMASK'] = momlib.mergemasks((1-source['SEXMASK'],\
               source['SEXMASKOTHER']))
            
        source.execpars["version"] = version
        if bool(source.execpars['makefits'][0]):
            source.wrap2mef()
        
        # If there's no pixel different from 0 in object, flag as "BLANK".
        nselected = len(num.where(source['MASK']==0)[0])
        
        nonblankprocent = 100.0 * len(num.where(((1-source['MASK']) * \
        source['STAMP']) != 0.)[0])/ nselected
        print '\nnonblankprocent = %.1f\n' % nonblankprocent
        
        if nonblankprocent < 90.:	    
            source['flags'] = addflag(source['flags'],allflags['BLANK'])
            isblank = True
            print '\n BLANK OBJECT!\n'	    
        else: isblank = False
        
        # Do verbose Graphical Output?
        
        dograph = execpars['dograph'][0] == 1
        if dograph: 
            source['figures'] = {}
            source['figcomms'] = {}
        
        # MANDATORY
        
        # mandatory: radial profile, petrosian radius, petrosian mask, 
        # 1st & 2nd order moments
        # radial returns radial profile (within MASK)
        
        source['BOXY'] = 0 # Elliptical apertures always.
        
        # Measure sky sigma: THIS SAVES SOME TIME, AS IT IS ALREADY
        # REQUIRED BY SEVERAL TASKS.
        
        sky_sigma = source.execpars['sigma_sky'][0]
        
        source.getsky()
	
        if source['SKY_MEDIAN'] != None and \
        source.execpars['useMANBACK'][0] == -1:
            print '\nMeasured background : %.2e\n' % source['SKY_MEDIAN']
            source['BACKGROUND'] = source['SKY_MEDIAN']
        
        if int(sky_sigma) == -1:
            sky_sigma = source['SKY_SIGMA']
            print 'sky_sigma = %.2e' % sky_sigma
            
            if sky_sigma != None : source.execpars['sigma_sky'] = [sky_sigma]
            del sky_sigma
        
        # RADIAL PROFILE
        tRAD_1 = time()
        if not isblank:
            
            if 'M_RADIAL' in execpars['toexec']: 
                if execpars['doFineRadial'][0]==0:
                    source.radial_v1(dograph=dograph)
                else:
                    source.radial_v3(dograph=dograph)
            
            tRAD_2 = time()
            print '%f seconds in making RADIAL profile' % (tRAD_2-tRAD_1,)
        
            # SAVE RADIAL PROFILE
            tRADSAVE_1 = time()
            if 'M_RADIAL' in execpars['toexec'] and \
            bool(source.execpars['saveradial'][0]) :
                try: imgid = father_img_name
                except NameError: imgid = \
                quitpath(source.execpars['father_img_name'])
                imgid = imgid[0:string.rfind(imgid,'.')]
                id = '%s' % source['name']
                radialfile = '%s_%s_RADIAL.txt' % (id,imgid)
                source.SaveRadial(radialfile)
                tRADSAVE_2 = time()
                print '%f seconds in saving RADIAL profile' % \
                (tRADSAVE_2-tRADSAVE_1,)
            
            # PETROSIAN returns petrosian radius, intensity and 
            # flux (within MASK)
            if 'M_PETRO' in execpars['toexec']:  
                #source.petrosian()
		source.petrosian2()
            
            # petromsk returns PETROSIAN MASK!
            if 'M_PETROMSK' in execpars['toexec']:
                source.petromsk()
            
            # Average Signal to Noise ratio
            if 'SNR' in execpars['toexec']:
                source.snr()
            
            # ellipticity parameters are computed before the mask may be
            # updated to Petrosian mask.
            # ellipse updates A, B, THETA, ELONGATION, ELLIP
            # (within MASK^SEGMAP)
            
            tELL_1 = time()
            if 'M_ELLIP' in execpars['toexec']:
                source.ellipse()
            tELL_2 = time()
            print '%f seconds in running ellipse' % (tELL_2-tELL_1,)
            
            tMOM_1 = time()
            if 'M_MOM' in execpars['toexec']:
                source.getmoments()
                tMOM_2 = time()
                print '%f seconds in running moments' % (tMOM_2-tMOM_1,)
        else:
            # RADIAL
            radiusflags = 0L
            radiusflags = addflag(radiusflags,allflags['NORADIAL'])
            source['M_RADIAL'] = {'radii':None,'cumulflx':None,\
            'intens':None,'npix':None,'npixout':None,\
            'radiusflags':radiusflags}
            # SAVE RADIAL
            source['radial_file'] = '0'
            # PETROSIAN
            source['R_PETRO'] = -99.0 ; source['I_PETRO'] = -99.0
            source['F_PETRO'] = -99.0
            # PETROMSK
            source['M_PETROMSK'] = None
            # SNR
            source['SNR'] = -99.
            # ELLIPSE
            source['M_A'] = -99.0 ; source['M_B'] = -99.0 
            source['M_THETA'] = -99.0 ; source['M_ELONG'] = -99.0 
            source['M_ELLIP'] = -99.0
            # Moments
            source['M_X'] = -99.0 ; source['M_Y'] = -99.0
            source['M_X2'] = -99.0 ; source['M_Y2'] = -99.0
            source['M_XY'] = -99.0
        
        # END MANDATORY
        
        # Use of petrosian mask
        # 
        if source.execpars['usepetro'][0] == 1 and \
        source['M_PETROMSK'] != None:
            if source['EXTMASK'] != None : 
                bunch_of_masks = (source['M_PETROMSK'],\
                source['SEXMASKOTHER'],source['EXTMASK'])
            else : 
                bunch_of_masks = (source['M_PETROMSK'],source['SEXMASKOTHER'])
            
            source['MASK'] = momlib.mergemasks(bunch_of_masks)
	    source['flags'] = addflag(source['flags'],allflags['USEPETROMSK'])
            
            # WARNING: 'THE PETROSIAN MASK MAY BE UNNOTICEDLY
            # TRUNCATED BY THE WINDOW!!'
        
        # OPTIONAL PARAMETERS
        
        if not isblank:
            
            # 'Basics' gets area, average intensity and total flux of object 
            # within MASK.
            
            if 'BASICS' in execpars['toexec']:
                source.Basics()
            
            # SECOND ORDER MOMENT (ALTERNATIVE TAKE)
            if 'M2' in execpars['toexec']:
                source.M2()
            
            # RADII WHICH CONTAIN SEVERAL FLUX RATIOS
            
            if 'RADII' in execpars['toexec']:
                #source.Radii()
		source.Radii2()
            
            # FLUX INSIDE A GIVEN RADIUS (APFLXRADIUS)
            
            if 'APFLX' in execpars['toexec'] and 'APFLXRADIUS' in source:
                source.ApFlx()
            
            # COORDINATES OF PEAK EMISSION
            
            # peak updates peak center (within MASK^SEGMAP)
            if 'M_PEAK' in execpars['toexec']: 
                source.peak() # Minimum boxwidth = 3 pix
              
            # GINI
            tG_1 = time()
            if 'M_GINI' in execpars['toexec'] : 
                source.gini(dograph)
                tG_2 = time()
                print '%f seconds in running Gini' % (tG_2-tG_1,)
            
            # ASYMMETRY
            
            tA_1 = time()
            if 'ASYM' in execpars['toexec'] : 
                source.asymmetry(dograph)
                tA_2 = time()
                print '%f seconds in running asymmetry' % (tA_2-tA_1,)
            
            # ANGULAR CONTRAST
            
            tAC_1 = time()
            if 'AC' in execpars['toexec']:
                source.AC(dograph)
                tAC_2 = time()
                print '%f seconds in running AC' % (tAC_2-tAC_1,)
            
            # CONCENTRATION
            
            tC_1 = time()
            if 'CONCENT' in execpars['toexec']:
                source.concent(dograph)
                tC_2 = time()
                print '%f seconds in running concent' % (tC_2-tC_1,)
            
            # FIND PEAKS (TO BE DROP SOON...)
            
            tNP_1 = time()
            if 'NPEAKS' in execpars['toexec']:
                source.FindPeaksII(dograph) # ON TESTS
                tNP_2 = time()
                print '%f seconds in running FindPeaks' % (tNP_2 - tNP_1,)
            
            # CLUMPS STATISTICS
            
            tNC_1 = time()
            if 'NCLUMPS' in execpars['toexec']:
                source.FindClumps(dograph) # ON TESTS
                tNC_2 = time()
                print '%f seconds in running FindClumps' % \
                (tNC_2-tNC_1,)
            
            # CLUMPINESS
            
            tCL_1 = time()
            if 'CLUMPY' in execpars['toexec'] and source['R_PETRO'] != -99:
                source.clumpy(dograph) # Minimum boxwidth = 3 pix
            else:
                source['M_S'] = -99.0 ; source['M_S_SKY'] = -99.0
                tCL_2 = time()
                print '%f seconds in running clumpy' % (tCL_2 - tCL_1,)
            
            # AXIS ASYMETRY (MAJOR AXIS)
            
            tMAXAXIS_1 = time()
            if 'MAJOR_SIM' in execpars['toexec']:
                source.axis_asymmetry(axis='major',dograph=dograph)
                tMAXAXIS_2 = time()
                print '%f seconds in MAJOR_SIM' % (tMAXAXIS_2-tMAXAXIS_1,)
             
            # AXIS ASYMMETRY (MINOR AXIS)
            
            tMINAXIS_1 = time()
            if 'MINOR_SIM' in execpars['toexec']:
                source.axis_asymmetry(axis='minor',dograph=dograph)
                tMINAXIS_2 = time()
                print '%f seconds in running MINOR_SIM' % \
                (tMINAXIS_2-tMINAXIS_1,)
             
            # M20
             
            tM20_1 = time()
            if 'M20' in execpars['toexec'] : 
                if execpars['M20mode'][0] == 'Lotz':
                    source.M20Lotz(dograph)
                elif execpars['M20mode'][0] == 'Azzo':
                    source.M20Azzo(dograph)
                tM20_2 = time() 
                print '%f seconds in running M20 (%s)' % (tM20_2-tM20_1,\
                execpars['M20mode'][0])
             
            # EXCENTRICITY
             
            tE_1 = time()
            if 'EXCENTRICITY' in execpars['toexec']:
                source.Excentricity(dograph)
                tE_2 = time()
                print '%f seconds in running Excentricity' % (tE_2-tE_1,)
            
            # FILLING FACTOR
            
            tFF_1 = time()
            if 'FFACTOR' in execpars['toexec']:
                source.FFactor(dograph)
                tFF_2 = time()
                print '%f seconds in running FFactor' % (tFF_2-tFF_1)
            
            # VISITORS
            ttr_1 = time()
            if 'TRUNC' in execpars['toexec']:
                scale = 0.03
                zeroT = 24.84315
                source.trunc(scale,zeroT,dograph)
                ttr_2 = time()
                print '%f seconds in running trunc' % (ttr_2-ttr_1)
            
        else : 
            # BASICS
            source['M_NPIX'] = -99 ; source['M_FLUX'] = -99.
            source['M_AVINT'] = -99.
            # M2
            source['M_M2'] = -99.
            # RADII
            source['M_R20'] = -99. ; source['M_R50'] = -99. ; 
            source['M_R80'] = -99. ;
            # APFLX
            source['M_APFLX'] = -99.
            # PEAK
            source['M_XPEAK'] = -99. ; source['M_YPEAK'] = -99.
            # GINI
            source['M_GINI'] = -99.
            # ASYMMETRY
            source['M_AS_X'] = -99. ; source['M_AS_Y'] = -99.
            source['M_AS'] = -99. ; source['M_AS_SKY'] = -99.
            # AC
            source['MAC8'] = -99. ; source['MAC8M'] = -99.
            source['MAC4'] = -99.
            # CONCENT
            source['M_C'] = -99.0
            # NPEAKS
            source['M_NPEAKS'] = -99
            # NCLUMPS
            source['M_NUM_CL'] = -99 
            source['M_MAX_CL'] = -99. ; source['M_MIN_CL'] = -99.
            source['M_ACC_CL'] = -99. ; source['M_FAR_CL'] = -99.
            # CLUMPY
            source['M_S'] = -99. ; source['M_S_SKY'] = -99.    
            # MAJOR_SIM
            source['M_AXS_MAJ'] = -99. ; source['M_AXS_SKY_MAJ'] = -99.
            # MINOR_SIM
            source['M_AXS_MIN'] = -99. ; source['M_AXS_SKY_MIN'] = -99.
            # getmoments
            source['M_X'] = -99. ; source['M_Y'] = -99.
            source['M_X2'] = -99. ; source['M_Y2'] = -99. ; source['M_XY'] = -99.
            # M20
            source['M20'] = -99.
            # EXCENTRICITY
            source['M_E'] = -99.
            # FFACTOR
            source['M_FF'] = -99.
        
        # SOME OTHER OPTIONAL GRAPHIC OUTPUTS
        
        if dograph:
            source.stamp_graph('STAMP','stamp')
            source.stamp_graph('MASK','mask')
            if 'SKYSTAMP' in source : 
                source.stamp_graph('SKYSTAMP','sky')
            source.petropeakcenter_graph()
            
            pdfid = source._getGraphId()
            latexfile = '%s.tex' % pdfid
            pdffile = '%s.pdf' % pdfid
            psfile = '%s.ps' % pdfid
            latex = LaTeX()
            figures = copy(source['figures'])
            figcomms = copy(source['figcomms'])
            header = '%s' % source['name']
            latex.DoBody(header,figures,figcomms)
            latex.Write(latexfile)
            latex.Compile(latexfile,cleanafter=True,\
            figures=copy(source['figures']))
            latex.Ps2Pdf(psfile,pdffile,cleanafter=True)
            
            source['PDF'] = pdffile
            
        # PACKAGING OF REMAINING DATA TO THE OUTPUT OBJECT
        for key in execpars['to_output']: 
            try:    
                father_cat[key][counter] = source[key]
            except : stop()
        
        del source # free memory
        
        #sys.exit()
        t4 = time()
        lapsus = t4 - t3
        print '%i seconds in analyzing 1 object\n\n\n' % lapsus
        
           
        # DUMPING OUTPUT to a FILE
        # it writes outputs to a file in sextractor fashion.
        # outfile_name = dumpcat(outfile_name,output,to_output)
        if (Ndump >=1) and (counter % Ndump == 0):
            if isthere(execpars['outfile_name']):
                os.system('rm %s' % execpars['outfile_name'])
            father_cat.dumpcat()
    
    if Ndump > 1:
        if isthere(execpars['outfile_name']):
            os.system('rm %s' % execpars['outfile_name'])
    
    father_cat.dumpcat()
    
    t2 = time()
    
    lapsus = t2 - t1
    
    lapsus_fr = secs_to_dhms(lapsus)
    
    if lapsus_fr[0]>0.:
        print """\n\n\n'Only' %i days, %i hours, %i minutes, %f seconds 
        in analyzing %i objects\n\n\n""" % \
        (lapsus_fr[0],lapsus_fr[1],lapsus_fr[2],lapsus_fr[3],nobj_out)
    else: 
        print """'Only' %i hours, %i minutes, %f seconds in analyzing 
        %i objects\n\n\n""" % (lapsus_fr[1],lapsus_fr[2],lapsus_fr[3],nobj_out)
    
    return None


###############################################################################################
## to execute the program from the console ####################################################

if __name__=="__main__":

    from optparse import OptionParser
    from pdb import set_trace as stop 
    
    parser = OptionParser()
    parser.add_option("-f","--file", dest="filename", default=None, help="Input file to read data from")

    (options, args) = parser.parse_args()

    if options.filename:
    	file = options.filename
    
    moments(file)
