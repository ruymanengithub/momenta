#! /usr/bin/env python 

import numpy as num

def read_inputs_T(infile):
    """Reads inputs for truncation.py"""
    # IMPORT STUFF
    from read_sex import read_filesex, read_param
    from pdb import set_trace as stop
    # END IMPORT STUFF   
    
    from verbose_truncation import mom_parfor as trunc_parfor
    
    def updatedict(child,parent):
        for key in parent.keys():
           if key not in child : child[key] = parent[key]
        return child
    
    def updatelist(child,parent):
        for item in parent:
           if item not in child : child.append(item)
        return child
    
    # MANDATORY
    exec_mandatory = [] 
    # always computed, by default
    sexpars_mandatory = ['XMIN_IMAGE','YMIN_IMAGE','XMAX_IMAGE',\
    'YMAX_IMAGE','NUMBER','THETA_IMAGE','ELLIPTICITY','X_IMAGE','Y_IMAGE']
    to_output_mandatory = ['NUMBER'] # always output, by default
    # END MANDATORY
    
    prim = read_filesex(infile,trunc_parfor)
    
    sexpars = read_param(prim['sexpars_f'][0])
    toexec = read_param(prim['toexec_f'][0])
    to_output = read_param(prim['to_output_f'][0])
    
    execpars = read_filesex(prim['execpars_f'][0],trunc_parfor)
    defaults = read_filesex(prim['defaults_f'][0],trunc_parfor)    
    
    execpars = updatedict(execpars,defaults)
    sexpars = updatelist(sexpars,sexpars_mandatory)
    toexec = updatelist(toexec,exec_mandatory)
    to_output = updatelist(to_output,to_output_mandatory)
    
    exec_dict = {}
    
    exec_dict['father_img_name'] = prim['father_img_name'][0]
    exec_dict['father_cat_name'] = prim['father_cat_name'][0]
    exec_dict['father_sex_name'] = prim['father_sex_name'][0]
    exec_dict['father_seg_name'] = prim['father_seg_name'][0]
    if 'father_mask_name' in prim:  
        exec_dict ['father_mask_name'] = prim['father_mask_name'][0]
    exec_dict['outfile_name'] = prim['outfile_name'][0]
    exec_dict['sexpars'] = sexpars
    exec_dict['toexec'] = toexec
    exec_dict['to_output'] = to_output
    for key in execpars : exec_dict[key] = execpars[key]
    
    return exec_dict # dictionary with all the inputs for moments.py
