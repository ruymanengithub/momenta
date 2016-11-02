#! /usr/bin/env python

import numpy as num
from pdb import set_trace as stop

def reg_format(flags):
    """Asigns formats to each column in a sdss_ready catalogue
    Needs improvements: Just updating.
    """

    #IMPORT STUFF
    import string
    #END IMPORT STUFF

    print "UPDATE THE DICTIONARY OF FORMATS! YOU LAZY BUM..."

    allflags = {'RUN':'Int32','CAMCOL':'Int32','RA':'Float32'}
  
    formats = {}							     
    for flag in flags: 						     
        formats[flag] = allflags[flag]				     

    return formats  			

def feedin(self):
    """Reads a catalog as a dictionary of type value:numpy.array."""
    # IMPORT STUFF
    import string
    from catlib import reg_format
    # END IMPORT
    
    catname = self.execpars['incat']
    
    # readout of catalogue

    file = open(catname,'r')
    all_lines = file.readlines() # all catalogue as a list of strings
    file.close()

   
    first2 = []
    
    for item in all_lines : first2.append(item[0:3])
    all_lines = [all_lines[n] for n in range(len(first2)) if first2[n] != '#!'] # drops lines which start by '#!'

    nlines = len(all_lines)		    # text lines in catalogue
    header = [header[n] for n in range(len(first2)) if first2[n] == '# '] # header info
    nheader = len(header)		    # text lines in header
    nobj = nlines - nheader		    # objects in catalogue
    ncols = len(string.split(all_lines[-1])) # object columns in catalogue
    
    
    cols = []
    
    for line in header: 
        cols.append(words[2])

    # readout of objects. Asignation of identities to columns, creation
    # of a record array with all columns and their names

    formats = reg_format(cols) # dictionary	   

    dummie = num.array(['0']*nobj,dtype='|S60')

    record = {}     # dictionary

    for col in cols : 
        record[col] = dummie
        if formats[col] != 'Char': record[col] = dummie.fasteval(type=formats[col])
    del dummie      
    del col

    for i in range(nobj):
    	items = string.split(all_lines[i+nheader])
    	for j in range(ncols): record[cols[j]][i] = num.array(items[j]).fasteval(type=formats[cols[j]])[0]

    for item in record : self[item] = record[item]
    
    return None
