#! /usr/bin/env python 

from pdb import set_trace as stop
import numpy as num

class momsource(dict):
    """Core class of moments.py"""
    from momlib import make_stamp, make_mask, wrap2mef,getwindow
    from radial import radial_v1,radial_v3,SaveRadial
    from petrosian import petrosian, petromsk, petrosian2
    from Gini import gini
    from Clumpiness import clumpy
    from Concentration import concent
    from algorithms import ellipse, peak, snr, makeskystamp, \
    getmoments, getsky
    #from AxisAsymmetry import axis_asymmetry
    #from Asymmetry import asymmetry
    #from M20 import M20Lotz, M20Azzo
    from development import FFactor, AC, Radii, Radii2, M2, Basics, \
    FindPeaksII, ApFlx
    from Excentricity import Excentricity
    from segmentation import FindClumps
    #from graphs_algorithms import gini_graph, M20_graph, \
    #axis_asymmetry_graph, asymmetry_graph, concent_graph, \
    #clumpy_graph, petropeakcenter_graph, radial_graph, \
    #FindPeaks_graph, FindClumps_graph, stamp_graph, _getGraphId,\
    #Excentricity_graph, AC_graph, FFactor_graph
    # VISITORS:
    from trunc import trunc
    from trunc_man import trunc_man
    
    def __init__(self,execpars):
        """Initializes momsource"""
        self.execpars = execpars
        self['flags'] = 0L
	

class momcat(dict):
    """Catalog-like class of moments.py"""
    from read_sex import read_catsex, dumpcat, initialize
    
    def __init__(self,execpars):
        """Initializes momcat"""
        self.execpars = execpars
        self.flags = {}
        self.ancillary = {}
