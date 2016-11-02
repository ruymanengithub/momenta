#! /usr/bin/env python

def clumps(infile):
    """Spin-off of moments.py to measure CLUMPS in galaxies."""
    #
    version = 0.0
    #ruler75###########################################
    # ALL THE IMPORT STUFF 
    import numpy as num
    import read_sex
    from pdb import set_trace as stop
    import string
    from flags_truncation import addflag,allflags,isflagon
    from clumps_verbose import cl_ordered, cl_comm, cl_regfor
    from momsource import momcat
    from clumps_lib import Clumpssource
    from time import time
    import os, sys
    import pyfits
    from CommonTools.loadascii import loadascii
    import momlib
    quitpath = momlib.quitpath
    secs_to_dhms = momlib.secs_to_dhms
    from clumps_lib3 import read_inputs_CL
    from Moments.latex import LaTeX
    from copy import copy
    from algorithms import get_stat
    # END IMPORT STUFF
    
    isthere = os.path.exists
    
    t1 = time()
    
    execpars = read_inputs_CL(infile) # execution parameters
    
    from clumps_verbose import cl_ordered
    from clumps_verbose import cl_comm
    from clumps_verbose import cl_regfor
    
    #ruler75##########################################
    
    # READOUT OF CATALOGUE
    
    father_cat = momcat(execpars) # initialize input catalog
    father_cat.read_catsex(cl_regfor)  # read input catalog
    
    for key in father_cat.flags:
        if father_cat.flags[key] == False :
            print """SOME REQUESTED PARAMETER NOT IN CATALOGUE: %s""" %\
            key
            stop()
    del key
    
    
    # Filtering of catalog 
    
    if 'DO' in father_cat and father_cat.execpars['useDo'][0] == 1:
        indxs = num.where(father_cat['DO'] == 1)
        NUMBERs = father_cat['NUMBER'][indxs]
        for key in father_cat:
            father_cat[key] = father_cat[key][indxs]
    
    if father_cat.execpars['delay'][0] != 0:
        delay = father_cat.execpars['delay'][0]
        for key in father_cat:
            father_cat[key] = father_cat[key][delay:]
    
    # number of objects to analyze
    nobj_out = len(father_cat[father_cat.keys()[0]])
    
    father_cat.initialize(cl_regfor) # initialize output catalog
    
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
        
        if 'SEGIMAGE' in father_cat and 'IMAGE' in father_cat :
            
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
        
        source = Clumpssource(execpars) 
        # inherits from dict class and momsource class.
        for col in father_cat : source[col] = father_cat[col][counter]		      
        del col
        source['flags'] = 0L
        
        if bool(source.execpars['useMANBACK'][0]):
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
        # if bool(source.execpars['makefits'][0]) :
        #     source.wrap2mef()
        
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
        
        # MANDATORY
        source['BOXY'] = 0 # Elliptical apertures always.
        
        # Measure sky sigma: THIS SAVES SOME TIME, AS IT IS ALREADY
        # REQUIRED BY SEVERAL TASKS.
        
        sky_sigma= source.execpars['sigma_sky'][0]
        source.getsky()

        if source['SKY_MEDIAN'] != None and \
        source.execpars['useMANBACK'][0] == -1: 
            print '\nMeasured background : %f\n' % source['BACKGROUND']
            source['BACKGROUND'] = source['SKY_MEDIAN']
        
        if int(sky_sigma) is -1:
            sky_sigma = source['SKY_SIGMA']
            print 'sky_sigma = %.2e' % sky_sigma
            
            if sky_sigma != None : source.execpars['sigma_sky'] = [sky_sigma]
            del sky_sigma
        
        # RADIAL PROFILE
        tRAD_1 = time()
        
        if not isblank:
            try: useInRadial = string.lower(execpars['inradial'][0]) != 'none'
            except KeyError : useInRadial = False
            
            if 'M_RADIAL' in execpars['toexec'] and not useInRadial: 
                if execpars['doFineRadial'][0]==0:
                    source.radial_v1(dograph=False)
                else:
                    source.radial_v3(dograph=False)
            else:
                inradial = execpars['inradial'][0]
                rcols = ['radii','cumulflx','ecumulflx','intens','eintens',
                'npix','npixout']
                rformats = 'f,f,f,f,f,i,i'
                M_RADIAL = loadascii(inradial,rcols,rformats,separator='blank')
                source['M_RADIAL'] = copy(M_RADIAL)
            
            
            tRAD_2 = time()
            if useInRadial: 
                print '%f seconds in READING RADIAL profile' % (tRAD_2-tRAD_1,)
            else:
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
                source.petrosian()
            
            # petromsk returns PETROSIAN MASK!
            if 'M_PETROMSK' in execpars['toexec']:
                source.petromsk()
            
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
            # ELLIPSE
            source['M_A'] = -99.0 ; source['M_B'] = -99.0 
            source['M_THETA'] = -99.0 ; source['M_ELONG'] = -99.0 
            source['M_ELLIP'] = -99.0
            # MOMENTS
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
        
        if not isblank:
            
            # 'Basics' gets area, average intensity and total flux of object 
            # within MASK.
            
            if 'BASICS' in execpars['toexec']:
                source.Basics()
            
            # SECOND ORDER MOMENT (ALTERNATIVE TAKE)
            if 'M2' in execpars['toexec']:
                source.M2()
            
            # STARRING
            
            tcl_1 = time()
            if 'CLUMPS' in execpars['toexec']:
                source.Clumpscore()
                tcl_2 = time()
                print '%f seconds in running clumps...' % (tcl_2-tcl_1)
        else :
            # BASICS
            source['M_NPIX'] = -99 ; source['M_FLUX'] = -99.
            source['M_AVINT'] = -99.
            # M2
            source['M_M2'] = -99.
            # CLUMPS
            source['Acl'] = -99. ; source['Fcl'] = -99. ; source['Ncl'] = -99.
            source['maxAcl'] = -99. ; source['minAcl'] = -99.
            source['medianAcl'] = -99. ; source['meanAcl'] = -99.
            source['sigmaAcl'] = -99. ; source['p25Acl'] = -99.
            source['p75Acl'] = -99. ; source['phiAcl'] = -99.
            source['maxFcl'] = -99. ; source['minFcl'] = -99.
            source['medianFcl'] = -99. ; source['meanFcl'] = -99.
            source['sigmaFcl'] = -99. ; source['p25Fcl'] = -99.
            source['p75Fcl'] = -99. ; source['phiFcl'] = -99.
            source['CLTHRESH'] = -99.  ; source['CLMINAREA'] = -99.
            source['CL_cat'] = 'None' ; source['CL_seg'] = 'None'
            source['pdf']  = 'None'

        # PACKAGING OF REMAINING DATA TO THE OUTPUT OBJECT
        
        for key in execpars['to_output']: father_cat[key][counter] = source[key]
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
            father_cat.dumpcat(cl_ordered,cl_comm,cl_regfor)
         
    
    if isthere(execpars['outfile_name']):
        os.system('rm %s' % execpars['outfile_name'])
    father_cat.dumpcat(cl_ordered,cl_comm,cl_regfor)

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
    
    clumps(file)
