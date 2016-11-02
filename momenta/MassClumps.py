#! /usr/bin/env python

"""Script to produce, for a set of galaxies, a Mass map, a radial profile 
of Mass, and a survey on Mass Clumps."""

# IMPORT STUFF
import numpy as num
from pdb import set_trace as stop
# END IMPORT

def MassClumps(infile):
    """Script to produce, for a set of galaxies, a Mass map, a radial profile 
    of Mass, and a survey on Mass Clumps."""
    version = 0.0
    #ruler75###########################################
    # ALL THE IMPORT STUFF 
    import read_sex
    import string
    from flags_MassClumps import addflag,allflags,isflagon
    from momsource import momsource, momcat
    from time import time
    import os, sys
    import pyfits
    from CommonTools.loadascii import loadascii
    import momlib
    quitpath = momlib.quitpath
    secs_to_dhms = momlib.secs_to_dhms
    from MassClumps_lib import read_inputs_MC
    from Moments.latex import LaTeX
    from copy import copy
    from algorithms import get_stat
    # END IMPORT STUFF
    
    isthere = os.path.exists
    
    t1 = time()
    
    execpars = read_inputs_MC(infile) # execution parameters
    
    from verbose_MassClumps import mom_ordered as MC_order
    from verbose_MassClumps import mom_comm as MC_comm
    from verbose_MassClumps import mom_regfor as MC_regfor
    
    #ruler75############################################
    
    # READOUT OF CATALOGUE
    
    father_cat = momcat(execpars)                    # initialize input catalog
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
    
    # LOADING IMAGES: NOW 4 IMAGES + A MASK.
    # TO CIRCUMVENT MEMORY EXHAUSTION, ALWAYS USE FITSCUT.
    
    extension = 0 # Image extension. Not ready to handle higer extensions.
    
    # LOOP OVER SOURCES AND SOURCE OBJECT CREATION
    
    try: Ndump = father_cat.execpars['Ndump'][0]
    except KeyError: Ndump = 0
    
    nbands = father_cat.execpars['nbands'][0]
    
    
    # MAIN LOOP
    for counter in range(nobj_out):
        
        father_img_names = []
        for i in range(nbands):
            father_img_names.append(father_cat['IMAGE_%i'%(i+1)][counter])
        father_seg_name = father_cat['SEGIMAGE'][counter]
        try: father_mask_name = father_cat['MASK'][counter]
        except KeyError: father_mask_name = 'None'
            
        # TIME CONTROL
        
        t3 = time()
        print '\nprocessing object %i of a total of %i\n' % \
        (counter+1,nobj_out) 
        talready = t3-t1
        talready_fr = secs_to_dhms(talready)
        if talready_fr[0]>1:
            print '... %i days, %i hours, %i minutes, %f seconds since start' % talready_fr[0:]
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
        
        source = MCsource(execpars) # inherits from dict class.
        for col in father_cat : source[col] = father_cat[col][counter]		      
        del col
        source['flags'] = 0L
        
        if bool(source.execpars['useMANBACK'][0]):
            source['BACKGROUND'] = source.execpars['MANBACK'][0]
        
        # WINDOW DEFINITION
        source.getwindow() # It may go off bounds... taken into account later.
        
        # MAKING STAMPS
        extension = 0 # default by now
        
        ! check here that window is in bounds
        source.make_stamp(imgname=father_seg_name,img=None,\
        name='SEGSTAMP',extension=extension,mode='fitscut')

        
        for i in range(nbands):
            father_img_name = father_img_names[i]
            source.make_stamp(imgname=father_img_name,
            img=None,name='STAMP_%i'%(i+1),extension=extension,
            mode='fitscut')
         
        if father_mask_name != 'None':
            source.make_stamp(imgname=father_mask_name,img=None,\
            name='EXTMASK',extension=extension,mode='fitscut')
        else : source['EXTMASK'] = None
        
        # "MASQUERADE..."
        
        # 1 is masked, 0 is non masked, always.
        
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
        
        
        # MANDATORY
        source['BOXY'] = 0 # Elliptical apertures always.
        
        # PROTAGONIST
        ttr_1 = time()
        if 'MC' in execpars['toexec']:
            source.MassClumps_F()
            ttr_2 = time()
            print '%f seconds in running MassClumps' % (ttr_2-ttr_1)
            
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
