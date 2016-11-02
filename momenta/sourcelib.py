#! /usr/bin/env python


def  get_stripe(run,rerun,camcol,field,release,version):
    """Gets StripeNum and mu0 for given run, rerun, camcol and field in an imaging version ('Best' or 'Target') 
       of a Data Release.
    """
    # IMPORT STUFF
    import string
    import re
    from internals import coverage_path
    from pdb import set_trace as stop
    import sys
    # END STUFF

    version = version.lower()
    release = str(release)
    
    file = 'tsChunk.dr%s.%s.par' % (release, version)
    file = coverage_path + file

    f = open(file)
    lines = f.readlines()
    f.close()
    
    pattern = "TSSEGLONG %s %s %s" % (str(run), str(camcol), str(rerun))
        
    # search in lines for the given run, camcol and rerun. Get the stripeNumber and startMu
    
    for line in lines :
        index = string.find(line,pattern)    
        if index != -1 : 
            words = string.split(line,' ')
            stripeNum = words[4]
            mu0 = words[6]       
            break

    try : mu0 = mu0
    except NameError: 
        stop()
        print 'An unexpected error happened'
        sys.exit()

    #columns in tsChunk.dr?.?.par file : 
    #TSSEGLONG, run, camCol, rerun, stripeNumber, strip[2], startMu, endMu, field0, nFields, status[20], 
    #startMuStripe, skyVersion

    return stripeNum, mu0


def getfilenames(self):
    """Gets the paths of relevant file names to 'sloansource' instances."""
    # IMPORT STUFF
    from sourcelib import get_stripe
    # END IMPORT
    
    root = self.execpars['paths']['root']
    filters = self.execpars['filters']
    
    run = self['RUN']
    camcol = self['CAMCOL']
    field = self['FIELD']
    release = self['RELEASE']
    version = self['VERSION']
    
    stripeNum, mu0 = get_stripe(run,rerun,camcol,field,release,version)
    
    run = str(run)
    rerun = str(rerun)
    camcol = str(camcol)
    field = str(field)
    stripeNum = str(tstripeNum)
    mu0 = str(mu0)

    fversion = version.lower()    
    if version.lower() == 'best' : ver_tag = '1'
    else :ver_tag = '0'

    # fpC
    
    frun = string.rjust(run,6).replace(' ','0')
    ffield = string.rjust(field,4).replace(' ','0')  
    fpC = {}
    for filter in filters :
        ffilter = filter+camcol
        interp = (root,run,rerun,camcol,category,frun,ffilter,ffield)
        fpC[filter] = 'imaging/%s/%s/corr/%s/%s-%s-%s-%s.fit.gz'  % interp

	# fpM

    fpM = {}
    for filter_ in filters_k:
        ffilter = filter_+camcol
        interp = (root,run,rerun,camcol,category,frun,ffilter,ffield)
        fpM[filter] = 'imaging/%s/%s/objcs/%s/%s-%s-%s-%s.fit' % interp

	    
	# psField

    interp = (root,run,rerun,camcol,category,frun,camcol,ffield)
    psField = 'imaging/%s/%s/objcs/%s/%s-%s-%s-%s.fit' % interp

	    
	# tsObj

    frerun = string.rjust(rerun,2).replace(' ','0')
    interp = (root,fversion,stripeNum,mu0,ver_tag,camcol,category,frun,camcol,frerun,ffield)
    tsObj = '%simaging/inchunk_%s/stripe%s_mu%s_%s/%s/%s-%s-%s-%s-%s.fit' % interp

    
    # tsField
    
    interp = (root,fversion,stripeNum,mu0,ver_tag,camcol,category,frun,camcol,frerun,ffield)
    tsField = 'imaging/inchunk_%s/stripe%s_mu%s_%s/%s/%s-%s-%s-%s-%s.fit' % interp
    
    # stripeNum, mu0, fpC*, fpM*, tsObj, psField, tsField 
    
    self['STRIPENUM'] = stripeNum
    self['MU0'] = mu0
    self['fpC'] = fpC
    self['fpM'] = fpM
    self['tsObj'] = tsObj
    self['psField'] = psField
    self['tsField'] = tsField
    
    for filter in filters:
        self['fpC%s'%filter] = fpC[filter]
        self['fpM%s'%filter] = fpM[filter]

    return None
    
def getbasicinfo(self):
    """Gets basic info from tsField files for each filter."""
    # IMPORT STUFF
    import pyfits
    from internals import exptime, scale # SDSS values
    # END IMPORT
    
    tsFname = self['tsField']
    filters = self.execpars['filters']
    
    tsField = pyfits.getdata(tsFname)
    
    toretrieve = ('quality','sky','skySig','aa','aaErr','kk','kkErr','airmass','gain','psf_width')
    flags = ('s','f','f','f','f','f','f','f','f','f') # Single value or Filter value
    orderfilters = ('u','g','r','i','z')
    
    # RAW VALUES
    
    retrieve = {}
    for index in range(toretrieve):
        if flags[index] == 's':
            retrieve[toretrieve[index]] = tsField.field(toretrieve[index])
        if flags[index] == 'f':
            for filter in filters:
                order = orderfilters.index(filter)
                retrieve['%s.%s' % (toretrieve[index],filter)] = tsField.field(toretrieve[index])[0][order]
    
    
    # PROCESSED VALUES
    
    sky = tsField.field('sky')
    skySig = tsField.field('skySig')
    aa = tsField.field('aa')
    kk = tsField.field('kk')
    airmass = tsField.field('airmass')
    psf_width = tsField.field('psf_width')
    
    
    zeropoints = - (aa + kk * airmass)
    sky_dn_pix = sky * exptime * (scale**2.0) * 10.0**(0.4*zeropoints)
    psf_pix = psf_width / scale
    
    proc = {'zeropoint':zeropoints,'sky_dn_pix':sky_dn_pix,'psf_pix':psf_pix}
    procflag = {'zeropoint':'f','sky_dn_pix':'f','psf_pix':'f'}
    
    for key in proc:
        if procflag[key] == 's':
            retrieve[key] = proc[key]
        if flags[index] == 'f':
            for filter in filters:
                order = orderfilters.index(filter)
                retrieve['%s.%s' % (key,filter)] = proc[key][0][order] # ?
    
    self['photofield'] = {}
    for key in retrieve: self['photofield'][key] = retrieve[key]
    
    return None
    
def getpsfmodel(self):
    """Gets psf models for each filter from psField files."""
    # IMPORT STUFF
    import os
    from internals import fieldxc, fieldyc
    # END IMPORT
    
    psFname = self['psField']
    filters = self.execpars['filters']
    orderfilters = ('u','g','r','i','z')
    objname = self['NAME']
    pathout = self.execpars['paths']['out']
    
    self['PSFMODS'] = {}
    
    if 'XCENTER' in self:
        xc = self['XCENTER']
    else: xc = fieldxc
    
    if 'YCENTER' in self:
       yc = self['YCENTER']
    else: yc = fieldyc
    
    for filter in filters : 
        outname = '%s%s_%s_psf.fit' % (pathout,objname,filter)
        filtercode = orderfilters.index(filter) + 1
        execline = 'read_PSF %s %s %s %s' % (psFname,filtercode,xc,yc,outname)
        os.system(execline)
        self['psfmods'][filter] = outname
        

def getmask(self):
    """Gets object masks for a given filter from fpM files."""
    # IMPORT STUFF
    from internals import masktype
    from sdsslib import quitpath
    import os
    # END IMPORT 
    
    fpMname = self['fpM'][filter]
    pathout = self.exepcars['paths']['out']
    objname = self['NAME']
    filter = self.execpars['maskfilter']
    
    outname = '%s%s_%s_msk.pol' % (pathout,objname,filter)
    execline= 'read_mask -p %s %s %s' % (fpMname,masktype,outname)
    os.system(execline)
    
    self['polymask'] = quitpath(outname)

def identify(self):
    """Identifies an object in the tsObj file to retrieve information from it."""
    
