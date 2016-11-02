"""Module for reading Sextractor files."""

# IMPORT STUFF
from CommonTools.IOtools import par_format, reg_comment
# END IMPORT

def read_catsex(self,biblioformat=-1,verbose=True):
    _verbose = verbose
    del verbose
    """Reads a list of columns (flags) from a Sextractor catalog
       flags: list of columns, fed in as a list
       catalog: name of catalog, fed in as a string
       records: dictionary with arrays of data
       flags_check: dictionary with the check of availability of each
       flag in the catalog (Boolean values).
       
    """
    # IMPORT STUFF
    import numpy as num
    import string
    from CommonTools.IOtools import par_format
    from pdb import set_trace as stop
    # END IMPORT STUFF
    
    if biblioformat is -1 : from verbose import mom_regfor as biblioformat
    
    # readout of catalog
    catalog = self.execpars['father_cat_name']
    flags = self.execpars['sexpars']
    
    file = open(catalog,'r')
    all_lines = file.readlines() # all catalog as a list of strings
    file.close()

    all_lines = [all_lines[n] for n in range(len(all_lines)) if all_lines[n][0:2] != '#!' and all_lines[n][0:2] != '\n']
    
    if all_lines[-1][0] == '#' :
        print '\nNo data on catalog!\n'
        raise RuntimeError
        return None

    nlines = len(all_lines)		    # text lines in catalog
    header = [all_lines[n] for n in range(len(all_lines)) if all_lines[n][0:2] == '# '] # header info
    nheader = len(header)		    # text lines in header
    nobj = nlines - nheader		    # objects in catalog
    ncols = len(string.split(all_lines[-1])) # object columns in catalog
    cols = ['empty'] * ncols
    
    head_comms = {}
    
    for line in header: 
        words = string.split(line)
        index = int(words[1])
        try: cols[index-1] = words[2]
        except IndexError: stop()
        head_comms[words[2]] = string.join(words[3:],' ')

    empties = [n for n in range(len(cols)) if cols[n] == 'empty']
    nonempties = [n for n in range(len(cols)) if cols[n] != 'empty']
    
    for i in range(ncols):
        if cols[i] == 'empty':
            lastnonempty = max([nonempties[n] for n in range(len(nonempties)) if nonempties[n]< i])
            colname = cols[lastnonempty]
            rootindex = colname.find('.')
            if rootindex == -1: rootindex = len(colname)
            colname = colname[0:rootindex]
            difference = i - lastnonempty
            if difference == 1:
                cols[lastnonempty] = colname+'.0'
                cols[i] = colname+'.1'
            else : 
               cols[i] = colname+'.'+str(difference)
    
    head_cols = []
    for col in cols:
        point = string.rfind(col,'.')
        if point == -1:
            head_cols.append(col)
        else:
            if col[point+1:] == '0' : head_cols.append(col[0:point])
            else : pass
    
    # readout of objects. Asignation of identities to columns, creation
    # of a record array with all columns and their names
    
    formats = par_format(cols,biblioformat,verbose=_verbose) # dictionary	   
    

    record = {}     # dictionary
    
    for col in cols :
        if formats[col] == 'Char' :
            record[col] = num.zeros((nobj,),dtype='|S%i' % 70)
            record[col][:] = '0'
        else :
            try:record[col] = num.zeros((nobj,),dtype=formats[col])
            except TypeError:stop('\n something wrong at Moments.read_sex\n')
            record[col][:] = 0
    del col
	
    for i in range(nobj):
        items = string.split(all_lines[i+nheader])
        
        for j in range(len(items)):
            if formats[cols[j]] == 'Char':
            	record[cols[j]][i] = items[j]
            else: 
                try: record[cols[j]][i] = \
                    num.array(items[j]).astype(formats[cols[j]])
                except : stop()
    
    # selection of columns according to flags

    sel_record = {}
    flags_check = {}
    multirecords = []
    for item in record.keys():
        point = string.rfind(item,'.')
        if point > -1 : multirecords.append(item[0:point])
    multirecords = num.array(multirecords)
    
    if flags == [] : flags = [col for col in cols if col != 'empty']
    
    for flag in flags:
        try:
            findit = num.where(multirecords == flag)
	    try: didfind = len(findit[0]) > 0
	    except IndexError: didfind = False
            if didfind:
                counter = 0
                all = 0
                while all != 1:
                    subflag = flag+'.%i' % counter
                    if subflag in record: 
                        sel_record[subflag] = record[subflag]
                        counter += 1
                    else : all = 1
            else:	       
                sel_record[flag] = record[flag]
        except KeyError:
            print "KeyError on %s" % flag
            flags_check[flag] = False
        else: flags_check[flag] = True  		
    
    for key in sel_record : self[key] = sel_record[key]
    for key in flags_check: self.flags[key] = flags_check[key]
    
    self.ancillary['COMMENTS'] = head_comms
    self.ancillary['ORDER'] = head_cols
    
def read_filesex(inputf,biblioformat=-1):
    """Reads a sex file and returns it as a dictionary with pairs 
    parameter:list"""
    # ALL IMPORT STUFF
    import string
    from pdb import set_trace as stop
    from CommonTools.IOtools import par_format
    # END IMPORT
    
    if biblioformat is -1 : from verbose import mom_parfor as biblioformat
    
    # Start of program
    
    file = open(inputf,'r')
    all_lines = file.readlines() # all file as a list
    file.close()

    dictionary = {}

    for line in all_lines:
        words = string.split(line)
        if len(words) >0:
            try: 
                comment_index = [word[0] for word in words].index('#')
                words = words[0:comment_index]
            except ValueError : pass
            
            if len(words) >0:
                param = words[0]
                rest = words[1:]
                rest = [string.replace(n,',',' ') for n in rest]
                rest = string.join(rest," ")
                rest = string.split(rest)	# is a list
                dictionary[param] = rest
        
    # type conversion
    
    formats = par_format(dictionary.keys(),biblioformat)
    
    for param in dictionary :
        format = formats[param]
        if format[0:3] == 'Int': dictionary[param] = [int(item) for item in dictionary[param]]
        try:
            if format[0:5] == 'Float': dictionary[param] = [float(item) for item in dictionary[param]]
        except ValueError: stop()
        if format[0:4] == 'Char': pass
        
    return dictionary

def obj_filter(catalog):
    """Filters a catalog, clipping al non associated objects."""

    # IMPORT STUFF
    import numpy as num
    # END IMPORT
    
    associated = num.where(catalog['NUMBER_ASSOC'][:] > 0)
    output = {}
    for col in catalog :
        output[col] = catalog[col][associated[0]]

    return output
    

def dumpcat(self,biblioorder=-1,bibliocomment=-1,biblioformat=-1,\
    bibliolength=-1,verbose=True):
    """Writes selected arrays in a dictionary to a sextractor-like catalog."""
    _verbose = verbose
    del verbose
    # IMPORT STUFF
    from os import access
    from os import F_OK
    import numpy as num
    import string
    from pdb import set_trace as stop
    from copy import copy
    from CommonTools.IOtools import reg_order,reg_comment,reg_length
    # END IMPORT STUFF
    
    # INPUTS    
    if biblioorder is -1 : from verbose import mom_ordered as biblioorder
    if bibliocomment is -1 : from verbose import mom_comm as \
    bibliocomment
    if biblioformat is -1 : from verbose import mom_regfor as biblioformat

    to_output = self.execpars['to_output']
    outname = self.execpars['outfile_name']
    
    # END INPUTS
    
    # ordering to_ouput : VERY IMPORTANT!
    # NUMBER must be first (for aesthetics)
    # THE other cols are ordered in order to handle correctly
    # VECTOR_ASSOC columns.
    
    # This is to easen transference of input columns to output catalog...
    
    if not isinstance(biblioorder,list): fbiblioorder = biblioorder.tolist()
    else : fbiblioorder = copy(biblioorder)
    for key in to_output :
        if key not in fbiblioorder : fbiblioorder.append(key)
    biblioorder = num.array(fbiblioorder)
    
    to_output = reg_order(to_output,biblioorder)
    
    # avoid overwriting of file
    
    isthere = access(outname,F_OK)
    
    if isthere == True :
        satisfaction = 0			
        inc = 0
        while satisfaction != 1 :
            incstr = str(inc)
            testname = outname + '.' + incstr
            isalsothere = access(testname,F_OK)
            if isalsothere == True :	    
                inc += 1 
            else : 
                satisfaction = 1
                outname = testname
                file = open(outname,'w')
    else : file = open(outname,'w')
    
    # WRITE HEADER (READY FOR ASSOC PARAMETERS!)
    
    comments = reg_comment(to_output,bibliocomment)
    
    counter = 1
    for key in to_output :
        ext = 0
        check = string.rfind(key,'.')
        if check >= 0 : 
            try : ext = int(key[check+1:])
            except ValueError : ext = key[check+1:]
            keyn = key[0:check]
            if isinstance(ext,str) : keyn = key
        else : keyn = key
     
     	if ext == 0 or check < 0 or isinstance(ext,str):
     	    colnum = string.rjust(str(counter),3)
     	    comment = ' ' + string.ljust(keyn,16)+comments[key]
     	    line = '# ' + colnum + comment
     	    print >> file, line
     	    counter += 1
     	else : counter+=1   
    
    del counter, key, keyn, colnum, line
        
    # WRITE OBJECTS
    
    try: nobj = len(self[self.keys()[0]])
    except TypeError : nobj = 1
    types = par_format(to_output,biblioformat,_verbose)
    
    if bibliolength != -1:
        blanks = reg_length(to_output,bibliolength)
    else:
        blanks = {}
        for key in to_output : blanks[key] = 0
    
    
    for key in to_output:
        typein = types[key]
        if typein == 'Char' and blanks[key] == 0:
            if key+'.1' in self: tofetch = key+'.1'
            else : tofetch = key    
            try: blanks[key] = max([len(reg) for reg in self[tofetch]])
            except ValueError: blanks[key] = 1
            except TypeError : stop()
    
    def dropzeros(number):
       sat = 0
       while sat is not 1:
           point = number.rfind('.')
           e = number.rfind('e')
           if point != -1 and e==-1:
                if number[-1] == '0' and (len(number) - point) > 2:
                    number = number[0:-1]
                else: 
                    sat = 1
                    return number
           else: 
                sat = 1
                return number
    
    
    for counter in range(nobj) :
        
        line = ''
        
        for key in to_output :
	
            typein = types[key]
            
            try : reg = self[key][counter]
            except TypeError : reg = self[key]
            if typein[0:3].lower() == 'int' :
                if reg != 0 : discr = num.log10(abs(reg))
                else : discr = 1
                tmp = '%'+'%i'%discr+'i'
                x = tmp % reg
                if blanks[key] == 0 : space = 13
                else : space = blanks[key]
                
            elif typein[0:5].lower() == 'float' :
                if reg != 0 : 
                   try: discr = abs(num.log10(abs(reg)))
                   except TypeError: stop()
                else : discr = 1
                if discr < 3 :
                    f = '%13.8f'
                    x = dropzeros(f % reg)
                    del f
                else : x = dropzeros('%9.7e' % reg)
                if blanks[key] == 0 : space = 13
                else : space = blanks[key]
            
            elif typein[0:4].lower() == 'char' :
                x = reg
                space = blanks[key]
            # print key, typein, reg, newcol
            
            elif typein[0] == '%':
                x = typein % reg
                if blanks[key] == 0 : space = 13
                else : space = blanks[key]            
            
            newcol = string.rjust(x,space)
            
            line += newcol + '  '
        
        print >> file, line
    
    try: del line, newcol, key, counter, x, reg, discr
    except NameError: pass
    file.close()
    
    return None
    

def read_param(inputf):
    """Reads a param file and returns it as a list
    neglects commented lines
    """

    # ALL IMPORT STUFF
    import string
    from pdb import set_trace as stop
    # END IMPORT  
    
    # Start of program
    
    file = open(inputf,'r')
    all_lines = file.readlines()  # all file as a list of strings
    file.close()

    params = []

    for line in all_lines:
        words = string.split(line)
        if len(words) >0:
            try: 
                comment_index = words.index("#")
                words = words[:comment_index]
            except ValueError : pass
            if len(words) >0:
                params.append(words[0])
    return params

def initialize(self,biblioformat=-1,verbose=True):
        """ """
        _verbose = verbose
        del verbose
        # IMPORT STUFF
        import numpy as num
        from pdb import set_trace as stop
        from CommonTools.IOtools import par_format        
        # END IMPORT
               
        if biblioformat == -1 : from verbose import mom_regfor as biblioformat
        
        to_output = self.execpars['to_output']
        to_output_formats = par_format(to_output,biblioformat,_verbose)
        nobj_out = len(self[self.keys()[0]])
        
        for item in to_output :
            
            if item not in self:
                if to_output_formats[item] == 'Char' :                     
                    self[item] = num.array(['0']*nobj_out,dtype='|S100')
                else :		    
                    self[item] = num.array([0]*nobj_out,dtype=to_output_formats[item])
	
