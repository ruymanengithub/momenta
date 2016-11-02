#! /usr/bin/env python

def truncation(infile):
    """Spin-off of moments.py to measure Truncation Radii in galaxies."""
    #
    version = 0.0
    #ruler75###########################################
    # ALL THE IMPORT STUFF 
    import numpy as num
    import read_sex
    from pdb import set_trace as stop
    import string
    from flags_truncation import addflag,allflags,isflagon
    from momsource import momsource, momcat
    from time import time
    import os, sys
    import pyfits
    from CommonTools.loadascii import loadascii
    import momlib
    quitpath = momlib.quitpath
    secs_to_dhms = momlib.secs_to_dhms
    from lib_truncation import read_inputs_T
    from Moments.latex import LaTeX
    from copy import copy
    from algorithms import get_stat
    # END IMPORT STUFF
    
    isthere = os.path.exists
    
    t1 = time()
    
    execpars = read_inputs_T(infile) # execution parameters
    
    from verbose_truncation import mom_ordered as trunc_order
    from verbose_truncation import mom_comm as trunc_comm
    from verbose_truncation import mom_regfor as trunc_regfor
    
    #ruler75############################################
    
    # READOUT OF CATALOGUE
    
    father_cat = momcat(execpars) # initialize input catalog
    father_cat.read_catsex(trunc_regfor)             # read input catalog
    
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
    
    nobj_out = len(father_cat[father_cat.keys()[0]]) 
    # number of objects to analyze
    
    father_cat.initialize(trunc_regfor) # initialize output catalog
    
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
        
        source = momsource(execpars) # inherits from dict class.
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
        
        
        if bool(source.execpars['makefits'][0]) :
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
        
        # MANDATORY
        source['BOXY'] = 0 # Elliptical apertures always.
        
        # Measure sky sigma: THIS SAVES SOME TIME, AS IT IS ALREADY
        # REQUIRED BY SEVERAL TASKS.
        
        sky_sigma= source.execpars['sigma_sky'][0]
        
        source.getsky()
        
        if int(sky_sigma) is -1:
            sky_sigma = source['SKY_SIGMA']
            print 'sky_sigma = %.2e' % sky_sigma
            
            if sky_sigma != None : source.execpars['sigma_sky'] = [sky_sigma]
            del sky_image,skymask,sky_sigma
            
        # END MANDATORY
        
        if not isblank:
            
            # 'Basics' gets area, average intensity and total flux of object 
            # within MASK.
            
            if 'BASICS' in execpars['toexec']:
                source.Basics()
            
            # SECOND ORDER MOMENT (ALTERNATIVE TAKE)
            if 'M2' in execpars['toexec']:
                source.M2()
            
            # PROTAGONIST
            ttr_1 = time()
            if 'TRUNC' in execpars['toexec']:
                try: manual = execpars['manual'][0] == 1
                except : manual = False
                if manual : source.trunc_man()
                else: source.trunc()
                ttr_2 = time()
                print '%f seconds in running trunc' % (ttr_2-ttr_1)
             
        else : 
            # BASICS
            source['M_NPIX'] = -99 ; source['M_FLUX'] = -99.
            source['M_AVINT'] = -99.
            # M2
            source['M_M2'] = -99.
            # TRUNC
            
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
            father_cat.dumpcat(trunc_order,trunc_comm,trunc_regfor)
         
    
    if isthere(execpars['outfile_name']):
        os.system('rm %s' % execpars['outfile_name'])
    father_cat.dumpcat(trunc_order,trunc_comm,trunc_regfor)

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
    
    truncation(file)
