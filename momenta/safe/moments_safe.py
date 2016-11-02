#! /usr/bin/env python

#def moments(): # ONLY FOR PROFILING!

def moments(infile) :
    """Program for image analysis, based on Sextractor photometry.
        It is aimed to compute several parameters of the spatial distribution
        of light in galaxies, not available in sextractor."""
    # infile = '/net/cipres/scratch1/THESIS/WORK/PROJECTIONS/first/SDSS/0.19/g/moments/infile_mom019g_0.txt'
    # PREVIOUS LINE ONLY FOR PROFILING!
    #
    # TO DO:
    #
    # IMPROVE RADIAL. THERE MIGHT BE EVEN SOME CONCEPTUAL "BUGS" IN IT! ????
    # THERE MIGHT BE SOME CONCEPTUAL BUG IN USING PIXELS WITH NEGATIVE VALUES IN THE 
    # ALGORITHMS.
    # CHECK "test_algor.py" FOR OTHER KNOWN PROBLEMS.
    
    version = 0.0

    #ruler75###################################################################     

    # ALL THE IMPORT STUFF 
    import numpy as num
    import read_sex
    import momlib
    from pdb import set_trace as stop
    import string
    from flags import addflag,allflags,isflagon
    from momsource import momsource, momcat
    from time import time
    import pyfits
    from CommonTools.loadascii import loadascii
    from momlib import quitpath
    # END IMPORT STUFF

    t1 = time()
    
    execpars = momlib.read_inputs(infile)

    #ruler75############################################

    # READOUT OF CATALOGUE
    
    father_cat = momcat(execpars)
    father_cat.read_catsex()

    for key in father_cat.flags:
        if father_cat.flags[key] == False :
            print 'SOME REQUESTED PARAMETER NOT IN CATALOGUE: %s' % key
            stop()
    del key
    
    # READOUT OF SEX FILE
    # dictionary with keywords:values

    father_sex = read_sex.read_filesex(execpars['father_sex_name'])     # dictionary. Not used by now.

    # FILTERING OF INPUT CATALOGUE: SELECTION OF OBJECTS TO ANALYZE
    #if ('clipassoc' in toexec) and ('NUMBER_ASSOC' in sexpars) and (flags_check['NUMBER_ASSOC'] == True):
    #    father_cat .obj_filter()
    
    # Filtering of catalog using an external file which selects which objects
    # in the SEXTRACTOR catalog are to be analyzed.
    
    if father_cat.execpars['NUMBERs_f'][0] != 'None' :
        NUMBERs_f = father_cat.execpars['NUMBERs_f'][0]
        cols = ['NUMBER']
        formats = 'i'
        load  = loadascii(NUMBERs_f,cols,formats,separator='blank')
        NUMBERs = load['NUMBER']
        indxs = []
        for number in NUMBERs:
            indxs.append(num.where(cat['NUMBER'] == number)[0])
        indxs = (num.array(indxs),)
        for key in father_cat:
            father_cat[key] = father_cat[key][indxs]
        
    elif 'DO' in father_cat and father_cat.execpars['useDo'][0] == 1:
        indxs = num.where(father_cat['DO'] == 1)
        NUMERs = father_cat['NUMBER'][indxs]
        for key in father_cat:
            father_cat[key] = father_cat[key][indxs]
    
    if father_cat.execpars['delay'][0] != 0:
        delay = father_cat.execpars['delay'][0]
        for key in father_cat:
            father_cat[key] = father_cat[key][delay:]
    
    nobj_out = len(father_cat[father_cat.keys()[0]])
    
    father_cat.initialize()
    
    # LOADING IMAGES:
    
    extension = 0 # default by now    
    
    if 'SEGIMAGE' not in father_cat and 'IMAGE' not in father_cat:
     
        father_img = pyfits.getdata(execpars['father_img_name'], ext=extension)
        father_seg = pyfits.getdata(execpars['father_seg_name'], ext=extension)
        try: father_mask = pyfits.getdata(execpars['father_mask_name'], ext=     extension)
        except KeyError : pass
        father_img_dim = father_img.shape
    
    # LOOP OVER SOURCES AND SOURCE OBJECT CREATION
    father_img_name_prev = ''
    
    for counter in range(nobj_out):
     
        # source object creation (f = father_cat)						      
        # source = {'NUMBER':f['NUMBER'][counter],\					      
        # 'X_WORLD':f['X_WORLD'][counter],...}						     
        
        if 'SEGIMAGE' in father_cat and 'IMAGE' in father_cat :
            if father_cat['IMAGE'][counter] != father_img_name_prev:
                father_img_name = father_cat['IMAGE'][counter]
                father_seg_name = father_cat['SEGIMAGE'][counter]
                try: father_mask_name = father_cat['MASK'][counter]
                except KeyError: father_mask_name = 'None'
                
                father_img = pyfits.getdata(father_cat['IMAGE'][counter],\
                ext= extension)
                father_img_dim = father_img.shape
                
                father_seg = pyfits.getdata(father_seg_name,\
                ext=extension)
                
                try: father_mask = pyfits.getdata(father_mask_name,\
                ext=extension)
                except IOError: pass
                
                father_img_name_prev = father_img_name
                
            else : pass
        
        t3 = time()
        print '\nprocessing object %i of a total of %i\n' % (counter+1,nobj_out) 
        print '%f seconds since start...\n' % (t3 - t1,)
        source = momsource(execpars) # inherits from dict class.
        for col in father_cat : source[col] = father_cat[col][counter]		      
        del col
        
        if bool(source.execpars['useMANBACK'][0]):
            source['BACKGROUND'] = source.execpars['MANBACK'][0]
        
        # WINDOW DEFINITION
        
        source.getwindow(father_img_dim)
        
        # MAKING STAMPS
        extension = 0 # default by now
        source.make_stamp(imgname=None,img=father_img,name='STAMP',extension=extension)
        source.make_stamp(imgname=None,img=father_seg,name='SEGSTAMP',extension=extension)
        try: source.make_stamp(imgname=None,img=father_mask,name='EXTMASK',extension=extension)
        except NameError : source['EXTMASK'] = None
        # MASKERADE...
        
        # 1 is masked, 0 is non masked
        source.make_mask(source['SEGSTAMP'],name='SEXMASK',mask_in=None,mode='nosky') # array.
        source.make_mask(source['SEGSTAMP'],name='SEXMASKOTHER',mask_in=None,mode='withsky') # array.
        if source['EXTMASK'] != None : 
            source['MASK'] = momlib.mergemasks((source['SEXMASK'],source['EXTMASK']))
            source['MASKOTHER'] = momlib.mergemasks((source['SEXMASKOTHER'],source['EXTMASK']))
            source['SKYMASK'] = momlib.mergemasks((1-source['SEXMASK'],source['SEXMASKOTHER'],\
            source['EXTMASK']))
        else : 
            source['MASK'] = source['SEXMASK']
            source['MASKOTHER'] = source['SEXMASKOTHER']
            source['SKYMASK'] = momlib.mergemasks((1-source['SEXMASK'],source['SEXMASKOTHER']))
            
        source.execpars["version"] = version
        if bool(source.execpars['makefits'][0]) :
            source.wrap2mef()
        
        # If there's no pixel different from 0 in object, flag as "BLANK".
        if not num.any(source['MASK'] * source['STAMP']): 
            source['flags'] = addflag(source['flags'],allflags['BLANK'])
            isblank = True
        else: isblank = False
        
        # mandatory: radial profile, petrosian radius, petrosian mask
        # radial returns radial profile (within MASK)
        source['BOXY'] = 0 # Elliptical apertures always.
        t3_1 = time()
        if not isblank:
            if 'M_RADIAL' in execpars['toexec']: source.radial_v1()
            #if 'M_RADIAL' in execpars['toexec']: source.radial_v2()
        else :
            radiusflags = 0L
            radiusflags = addflag(radiusflag,allflags['NORADIAL'])
            source['M_RADIAL'] = {'radii':None,'cumulflx':None,\
            'intens':None,'npix':None,'npixout':None,'radiusflags':\
            radiusflags}
        t3_2 = time()
        print '%f seconds in making RADIAL profile' % (t3_2-t3_1,)
        
        # Save radial profile
        t3_3 = time()
        if 'M_RADIAL' in execpars['toexec'] and \
        bool(source.execpars['saveradial'][0]) :
            try: imgid = father_img_name
            except NameError: imgid = quitpath(source.execpars['father_img_name'])
            imgid = imgid[0:string.rfind(imgid,'.')]
            id = '%s' % source['name']
            radialfile = '%s_%s_RADIAL.txt' % (id,imgid)
            source.SaveRadial(radialfile)
        t3_4 = time()
        print '%f seconds in saving RADIAL profile' % (t3_4-t3_3,)
        # petrosian returns petrosian radius, intensity and flux (within MASK)
        if 'M_PETRO' in execpars['toexec']:  source.petrosian()
        # petromsk returns petrosian mask!
        if 'M_PETROMSK' in execpars['toexec']: source.petromsk()
        # self['M_PETROMSK'] = petro_mask.astype('Int8')
        if not isblank:
            if 'SNR' in execpars['toexec']: source.snr()
        else : self['SNR'] = 0.
        # end mandatory
        
        # ellipticity parameters are computed before the mask may be updated to Petrosian mask...
        # ellipse updates A, B, THETA, ELONGATION, ELLIP (within MASK^SEGMAP)
        t3_4 = time()
        if not isblank:
            if 'M_ELLIP' in execpars['toexec'] : source.ellipse()
        else:
            self['M_A'] = -99. ; self['M_B'] = -99. ; self['M_THETA'] = -99
            self['M_ELONG'] = -99 ; self['M_ELLIP'] = -99
        t3_5 = time()
        print '%f seconds in running ellipse' % (t3_5-t3_4,)
        # Use of petrosian mask?
        # 
        if source.execpars['usepetro'][0] == 1 and source['M_PETROMSK'] != None:
            if source['EXTMASK'] != None : 
                bunch_of_masks = (source['M_PETROMSK'],source['SEXMASKOTHER'],source['EXTMASK'])
            else : 
                bunch_of_masks = (source['M_PETROMSK'],source['SEXMASKOTHER'])
            source['MASK'] = momlib.mergemasks(bunch_of_masks)
            source['flags'] = addflag(source['flags'],allflags['USEPETROMSK'])
            
            # print 'THE PETROSIAN MASK MAY BE UNNOTICEDLY TRUNCATED BY THE WINDOW!!'
            # print 'THIS SHOULD BE CORRECTED... BUT DO NOT MESS IT UP, PLEASE!'
        
        # optional: ellipse, peak, gini, asymmetry, concent, clumpy, axis_symmetry
        
        # peak updates peak center (within MASK^SEGMAP)
        if 'M_PEAK' in execpars['toexec']: source.peak() # Minimum boxwidth = 3 pix
        t3_6 = time()
        if 'M_GINI' in execpars['toexec'] : source.gini()
        t3_7 = time()
        print '%f seconds in running gini' % (t3_7-t3_6,)
        
        t3_8 = time()
        if 'ASYM' in execpars['toexec'] : source.asymmetry()
        t3_9 = time()
        print '%f seconds in running asymmetry' % (t3_9-t3_8,)
        
        t3_10 = time()
        if 'CONCENT' in execpars['toexec']: source.concent() #
        t3_11 = time()
        print '%f seconds in running concent' % (t3_11-t3_10,)
        
        t3_12 = time()
        if 'CLUMPY' in execpars['toexec'] and source['R_PETRO'] != -99 : source.clumpy() # Minimum boxwidth = 3 pix
        else : 
            source['M_S'] = -99. ; source['M_S_SKY'] = -99.
        t3_13 = time()
        print '%f seconds in running clumpy' % (t3_13-t3_12,)
        
        t3_14 = time()
        if 'MAJOR_SIM' in execpars['toexec'] : source.axis_asymmetry(axis='major')
        t3_15 = time()
        print '%f seconds in MAJOR_SIM' % (t3_15-t3_14,)
        
        # source['MAJOR_SIM'], flags = axis_symmetry(source,execpars,axis='major')
        t3_16 = time()
        if 'MINOR_SIM' in execpars['toexec']: source.axis_asymmetry(axis='minor')
        # source['MINOR_SIM'], flags = axis_symmetry(source,execpars,axis='minor')
        t3_17 = time()
        print '%f seconds in running MINOR_SIM' % (t3_17-t3_16,)
        t3_18 = time()
        if 'M_MOM' in execpars['toexec']: source.getmoments()
        t3_19 = time()
        print '%f seconds in running moments' % (t3_19-t3_18,)
        t3_20 = time()
        if 'M20' in execpars['toexec'] : source.M20()
        t3_21 = time() 
        print '%f seconds in running M20' % (t3_21-t3_20,)
        # PACKAGING OF REMAINING DATA TO THE OUTPUT OBJECT
        for key in execpars['to_output']: father_cat[key][counter] = source[key]
        del source # free memory
        
        t4 = time()
        lapsus = t4 - t3
        print '%i seconds in analyzing 1 object\n\n\n' % lapsus
        
    # DUMPING OUTPUT to a FILE
    # it writes outputs to a file in sextractor fashion.
    # outfile_name = dumpcat(outfile_name,output,to_output)
    
    father_cat.dumpcat()

    t2 = time()
    
    lapsus = t2 - t1
    
    print '\n\nThis is amazing!'
    print '"Only" %i seconds in analyzing %i objects\n\n\n' % (lapsus,nobj_out)
    
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
