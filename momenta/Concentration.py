#! /usr/bin/env python

def concent(self,verbose=False,dograph=False):
    """Updates in self the 'Concentration', acording to the definition given 
    in Bershady et al., 2000 (extracted from Lotz et al., 2004)
    
    C = 5 log (r80/r20)
    where r? is the radius of the aperture that contains ? % of the 
    total flux of the object. The total flux is that contained within 1.5rp of 
    the galaxy's center, as in Lotz et al., 2004 and Conselice 2003... 
    modified: is the asymptotic value of the flux growth curve. How is 
    the center determined? Lotz et al. use that given by the determination 
    of the 'A' parameter.
    """
    # IMPORT STUFF
    import scipy.interpolate as interp
    from flags import addflag, isflagon, allflags
    from pdb import set_trace as stop
    # END IMPORT
    
    # INPUTS
    rprof = self['M_RADIAL']
    r_big_proc = self.execpars['r_big_proc'][0]
    r_small_proc = self.execpars['r_small_proc'][0]
    # END INPUTS
    
    if not isflagon(self['flags'],allflags['NORADIAL']):
        cumulflx = rprof['cumulflx']
        radii = rprof['radii']
        normal = cumulflx/cumulflx[-1] * 100.
        print 'Min radius = %f pixels' % radii[0]
        
        def lazy(normal,radii,proc,addflag):
            to_solve_f = normal - proc
            to_solve_s = interp.splrep(radii,to_solve_f,k=1,s=0)
            roots = interp.sproot(to_solve_s)
            radius = -1
            flag = 0L
            if len(roots) == 0:  flag = addflag(flag,1L) # no solution
            else : 
                radius = max(roots)
                if len(roots) > 1: flag = addflag(flag,2L) # more than one radius.
            
            return radius, flag
        
        r_small, flag_small = lazy(normal,radii,r_small_proc,addflag)
        r_big, flag_big = lazy(normal,radii,r_big_proc,addflag)
        
        cflag = 0L
        if isflagon(flag_small,1L) or isflagon(flag_big,1L) : cflag = \
        addflag(cflag,allflags['NOCONC'])
        if isflagon(flag_small,2L) or isflagon(flag_big,2L):
            cflag = addflag(cflag,allflags['NOCONC'])
            cflag = addflag(cflag,allflags['MANYCONC'])
        C = -99.0
        if not isflagon(cflag,allflags['NOCONC']) : \
        C = 5. * num.log10(r_big / r_small)
        
        print 'C = %5.3f\n' % C
        self['M_C'] = C
        if not isflagon(cflag,allflags['NOCONC']) and\
        not isflagon(cflag,allflags['MANYCONC']) : 
            self['M_RBIG'] = r_big
            self['M_RSMALL'] = r_small
        else:
            self['M_RBIG'] = -99.0
            self['M_RSMALL'] = -99.0
        self['flags'] = addflag(self['flags'],cflag)
    else : 
        self['M_C'] = -99.0
        self['flags'] = addflag(self['flags'],allflags['NOCONC'])
    
    if dograph and not isflagon(self['flags'],allflags['NOCONC']):
        self.concent_graph()
    else : pass
   
    if verbose : return C,r_small,r_big
