#! /usr/bin/env python

from pdb import set_trace as stop

def moments_condor(infile):
    """This script divides an execution of moments in smaller bits to 
    feed Condor, and do it in several machines at a time."""
    # IMPORT STUFF
    from read_sex import read_filesex
    from pdb import set_trace as stop
    from momsource import momcat
    import string
    from moments_condor import write_infile_chunk,write_submit
    from moments_condor import write_sexcat_chunk
    import os
    # END IMPORT
    
    print 'THIS SCRIPT WORKS... BUT DOESNT ALLOW TO USE moments.py UNDER CONDOR'
    print 'IT PROVIDES A USELESS SOLUTION... PITY'
    sys.exit()
    
    # INPUTS
    #mode = 'DIVIDE' # DIVIDE, UNITE
    #moments_infile = 'infile_momGOODSz.txt'
    #ntasks = 100
    #condor_guide = 'momGOODSz_condor.txt'
    # END INPUTS
    
    #    parse inputs file
    
    # if mode == DIVIDE
    
    #    divide sextractor catalog in pieces
    
    #    create new inputs files
    
    #    execute condor
    
    #    dump file to use by UNITE to reorder all the information
    
    # if mode == UNITE:
    
    #     check wether all parts of the task where done
    
    #     if last == True : Do final file
    
    prim = read_filesex(infile)
    mode = prim['mode'][0]
    mom_infile = prim['moments_infile'][0]
    ntasks = int(prim['ntasks'][0])
    condor_guide = prim['condor_guide'][0]
    condor_submit = prim['condor_submit'][0]
    executable = prim['executable'][0]
    f_c_n = read_filesex(mom_infile)['father_cat_name'][0]
    m_c =  read_filesex(mom_infile)['outfile_name'][0]
    
    if mode == 'DIVIDE':
        # Load sextractor catalog
        fcat = open(f_c_n,'r')
        allcat = fcat.readlines()
        catheader = [line for line in allcat if line[0] == '#']
        lenheader = len(catheader)
        catdata = [line for line in allcat if line[0] != '#']
        fcat.close()
        
        # Divide sextractor catalog in pieces
        nchunks = len(catdata) / ntasks
        restchunk = len(catdata) % ntasks
        
        f_c_n_root = f_c_n[0:string.rfind(f_c_n,'.')]
        f_c_n_ext = f_c_n[string.rfind(f_c_n,'.'):]
        mom_infile_root = mom_infile[0:string.rfind(mom_infile,'.')]
        mom_infile_ext = mom_infile[string.rfind(mom_infile,'.'):]
        
        m_c_root = m_c[0:string.rfind(m_c,'.')]
        m_c_ext = m_c[string.rfind(m_c,'.'):]
        
        fguide = open(condor_guide,'a')
        print >> fguide,'%s    %s    %s' % \
        ('# sexcat_chunk','infile_chunk','outcat_chunk')
        
        for i in range(nchunks):
            
            f_c_chunk = '%s_condor%i%s' % (f_c_n_root,i,f_c_n_ext)
            mom_infile_chunk = '%s_condor%i%s' % \
            (mom_infile_root,i,mom_infile_ext)
            m_c_chunk = '%s_condor%i%s' % (m_c_root,i,m_c_ext)
            
            # Write guide-file with results from this task    
            print >> fguide,'%s   %s   %s' % \
            (f_c_chunk,mom_infile_chunk,m_c_chunk)
            
            # Creates individual input files for moments.py
            write_infile_chunk(mom_infile_chunk,mom_infile,\
            f_c_chunk,m_c_chunk)
            
            # Creates individual input catalogs for moments.py
            startline = i * ntasks + lenheader
            if i<nchunks-1 : endline = (i+1)*ntasks + lenheader
            else : endline = len(allcat)
            
            write_sexcat_chunk(f_c_chunk,allcat,startline,endline)
        
        fguide.close()
         
        # Create .submit file for condor
        argument = '%s_condor' % (mom_infile_root,)
        write_submit(condor_submit,executable,argument,nchunks)    
        
        # execute Condor
        #os.system('condor_submit %s' % condor_submit)
        
    elif mode == 'UNITE':
        pass
    
def write_infile_chunk(mom_infile_chunk,mom_infile,f_c_chunk,m_c_chunk):
            """Writes an input file for moments 
            (dividing a moments.py-task for condor...)"""
            # IMPORT STUFF
            import numpy as num
            # END IMPORT
            fin = open(mom_infile,'r')
            all_in_lines = fin.readlines()
            fin.close()
            f_c_n_line = num.where(num.array(['father_cat_name' in \
            line for line in all_in_lines]))[0][0]
            m_c_line = num.where(num.array(['outfile_name' in \
            line for line in all_in_lines]))[0][0]
            
            all_in_lines[f_c_n_line] = 'father_cat_name     %s\n' % f_c_chunk
            all_in_lines[m_c_line] = 'outfile_name      %s\n' % m_c_chunk
            
            fout = open(mom_infile_chunk,'w')
            for line in all_in_lines : print >> fout,line[0:-1]
            fout.close()
            
def write_submit(condor_submit,executable,argument,nchunks):    
    """Writes .submit files for condor."""
    
    fsubmit = open(condor_submit,'w')
    submit = """####################################################################
#                                                                  #
# condor.submit                                                    #
#                                                                  #
####################################################################

Executable = %s
Universe   = vanilla
Requirements = Memory >= 256 && OpSys == "Linux" && Arch == "INTEL" && Machine != "temple.ll.iac.es"
Rank       = Memory >= 64
arguments  = -f %s$(Process).txt
output     = out$(Process).txt
error      = err$(Process).txt
Log        = condor.log
getenv     = True
notification = Never
Queue %i
""" % (executable,argument,nchunks)

    print >>fsubmit,submit
    fsubmit.close()
    
    
def write_sexcat_chunk(f_c_chunk,allcat,startline,endline):
    """Writes individual catalogs as inputs for moments.py 
    (under Condor...)"""
    
    header = [line for line in allcat if line[0] == '#']
    data = allcat[startline:endline]
    
    fout = open(f_c_chunk,'a')
    for line in header: print >> fout,line[0:-1]
    for line in data : print >> fout,line[0:-1]
    fout.close()
    
    
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
    
    moments_condor(file)
