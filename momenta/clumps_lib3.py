#! /usr/bin/env python

# IMPORT STUFF
import numpy as num
from pdb import set_trace as stop
import pyfits
from time import time
import os
from momsource import momsource
from time import time
from pylab import figure,close,show
# END IMPORT

isthere = os.path.exists

def read_inputs_CL(infile):
    """Reads inputs for clumps.py"""
    # IMPORT STUFF
    from read_sex import read_filesex, read_param
    from pdb import set_trace as stop
    # END IMPORT STUFF
    
    from clumps_verbose import cl_parfor
    
    def updatedict(child,parent):
        for key in parent.keys():
           if key not in child : child[key] = parent[key]
        return child
    
    def updatelist(child,parent):
        for item in parent:
           if item not in child : child.append(item)
        return child
    
    # MANDATORY
    exec_mandatory = ['M_RADIAL','M_PETRO','M_PETROMSK','M_ELLIP','M_MOM'] 
    # always computed, by default
    sexpars_mandatory = ['XMIN_IMAGE','YMIN_IMAGE','XMAX_IMAGE',\
    'YMAX_IMAGE','NUMBER','THETA_IMAGE','ELLIPTICITY','X_IMAGE','Y_IMAGE']
    to_output_mandatory = ['NUMBER'] # always output, by default
    # END MANDATORY
    
    prim = read_filesex(infile,cl_parfor)
    
    sexpars = read_param(prim['sexpars_f'][0])
    toexec = read_param(prim['toexec_f'][0])
    to_output = read_param(prim['to_output_f'][0])
    
    execpars = read_filesex(prim['execpars_f'][0],cl_parfor)
    defaults = read_filesex(prim['defaults_f'][0],cl_parfor)    
    
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

def ClRadialPlot(self,ptype='intens',eps='',figid=1):
    """Does a Radial Plot... 'very' specific to this task."""
    # IMPORT STUFF
    from pylab import errorbar,xlabel,ylabel,title,axhline,xlim,ylim,savefig,plot
    # END IMPORT
    # INPUTS
    scale = self.execpars['scale'][0]
    magzero = self.execpars['magzero'][0]
    areaobject = float(len(num.where(self['SEXMASK']==0)[0]))
    ellipgal = self['ELLIPTICITY']
    skyradius = (areaobject / (num.pi * (1.-ellipgal)))**0.5
    skyradius *= self.execpars['petrofactor'][0] 
    clsubradial = self.execpars['clsubradial'][0] == 1
    try : doRadialModel = self.execpars['doRadialModel'][0] == 1
    except KeyError : doRadialModel = False
    # skyradius a 50% bigger than object radius, + or -
    # END INPUTS
    
    def domag(value):
        if isinstance(value,float):
            if value > 0.: return magzero - 2.5 * num.log10(value)
            else:  return -99.0
        else:
            pos = num.where(value > 0.)
            mag = value * 0. - 99.0
            try : mag[pos] = magzero - 2.5 * num.log10(value[pos])
            except : pass
            return mag
    
    def domu(value):
        if isinstance(value,float):
            if value > 0. : return magzero - 2.5 * num.log10(value/scale**2.)
            else: return -99.0
        else:
            pos = num.where(value > 0.)
            mu = value * 0. - 99.0
            try: mu[pos] = magzero - 2.5 * num.log10(value[pos]/scale**2.)
            except : pass
            return mu

    def dologarc2(value):
        if isinstance(value,float):
            if value > 0. : return num.log10(value * scale**2.)
            else : return -99.0
        else:
            pos = num.where(value > 0.)
            logarc2 = value * 0. - 99.0
            try: logarc2[pos] = num.log10(value[pos] * scale**2.)
            except : pass
            return logarc2
    
    if ptype in ['clfluxes','clareas']:
        try: test = self['CL_catdata']['RADII'] * scale
        except : return None
    
    if ptype == 'intens':
        def getXYplot(groupkey='M_RADIAL',errkey=''):
            xplot = self[groupkey]['radii'] * scale
            intens  = self[groupkey]['intens'].copy()
            if errkey != '': eintens = self[groupkey][errkey].copy()
            pos = num.where(intens>0.)
            xplot = xplot[pos] ; intens = intens[pos] 
            yplot = domu(intens)
            if errkey != '': 
                eintens = eintens[pos]
                yerr = 2.5 * num.log10(1.+eintens/intens)
            del pos
            if errkey != '': return xplot,yplot,yerr
            else: return xplot,yplot
        
        xplot,yplot,yerr = getXYplot('M_RADIAL','eintens')
        
        if clsubradial and doRadialModel:
            xplotmod,yplotmod = getXYplot('M_MOD_RADIAL')
        
        plotdata = {'xlabel':'arcsecs','ylabel':'mag/arcsec2',\
        'title':'Intensity Profile'}
        ms = 1
        xlimval = (0,skyradius*scale)
        ylimval = (yplot[num.where(xplot<=xlimval[1])].max(),\
        yplot[num.where(yplot != -99.0)].min())
        
    if ptype == 'clfluxes':
        xplot = self['CL_catdata']['RADII'] * scale
        FLUX_ISO = self['CL_catdata']['FLUX_ISO']
        pos = num.where(FLUX_ISO>0.)
        xplot = xplot[pos]
        FLUX_ISO = FLUX_ISO[pos]
        yplot = domag(FLUX_ISO)
        yerr = yplot * 0.
        plotdata = {'xlabel':'arcsecs','ylabel':'mag',\
        'title':'Clumps Fluxes'}
        ms = 8
        if self['meanFcl'] > 0.:
            plotdata['mean'] = domag(self['meanFcl'])
        if self['medianFcl'] > 0.:
            plotdata['median'] = domag(self['medianFcl'])
        if self['p25Fcl'] > 0.:
            plotdata['p25'] = domag(self['p25Fcl'])
        if self['p75Fcl'] > 0.:
            plotdata['p75'] = domag(self['p75Fcl'])
        if len(xplot)>0:
            xlimval = (0,skyradius*scale)
            ylimval = (yplot.max()+0.5,yplot.min()-0.5)
        else : xlimval = (0.,1) ; ylimval = (-1.,1.)
        
    if ptype == 'clareas':
        xplot = self['CL_catdata']['RADII'] * scale
        yplot = dologarc2(self['CL_catdata']['ISOAREA_IMAGE'])
        yerr = yplot * 0.
        plotdata = {'xlabel':'arcsecs','ylabel':'log(arcsecs2)',\
        'title':'Clumps Areas'}
        ms = 8
        if self['meanAcl'] > 0.:
            plotdata['mean'] = dologarc2(self['meanAcl'])
        if self['medianAcl'] > 0.:
            plotdata['median'] = dologarc2(self['medianAcl'])
        if self['p25Acl'] > 0.:
            plotdata['p25'] = dologarc2(self['p25Acl'])
        if self['p75Acl'] > 0.:
            plotdata['p75'] = dologarc2(self['p75Acl'])
        xlimval = (0,skyradius*scale)
        yrange = yplot.max() - yplot.min()
        ylimval = (yplot.min()-0.1*yrange,yplot.max()+0.1*yrange)
    
    figure(figid)
    if ptype == 'intens':
        errorbar(xplot,yplot,yerr=yerr,xerr=None,fmt='b.',ecolor=None)
        if clsubradial and doRadialModel:
            plot(xplotmod,yplotmod,'kx--')
    else: 
        if len(xplot)>0: plot(xplot,yplot,'bo',ms=ms)
    
    if plotdata != None:
        xlabel(plotdata['xlabel'])
        ylabel(plotdata['ylabel'])
        title(plotdata['title'])
        if 'mean' in plotdata:
            axhline(plotdata['mean'],color='k',linestyle='--')
        if 'median' in plotdata :
            axhline(plotdata['median'],color='b',linestyle='-')
        if 'p25' in plotdata:
            axhline(plotdata['p25'],color='r',linestyle='--')
        if 'p75' in plotdata:
            axhline(plotdata['p75'],color='r',linestyle='--')
    
    ylim(ylimval)
    xlim(xlimval)
    
    if eps != '' : savefig(eps)
    else : show()
    
    if eps != '':  close(figid)
    
    return None
    
def ClDisplay(self,img,epsn='',figid=1):
    """Makes .eps with an image. Specific to this task"""
    # IMPORT STUFF
    from Moments.graphs_local import Paint, doStamp, Convert
    from pylab import imshow
    # END IMPORT
    timetag = '%f' % time()
    tmppng = 'tmpCl_%s.png' % timetag
    tmpjpg = 'tmpCl_%s.jpg' % timetag
    if isthere(tmppng) : os.system('rm %s' % tmppng)
    if isthere(tmpjpg) : os.system('rm %s' % tmpjpg)
    
    stampformat='PNG'
    stamppars = '--palette=rainbow -l --min=%.e --max=%.e' % \
    (0.,img.max())
    doStamp(img,tmppng,format=stampformat,
    parameters=  stamppars)
    Convert(tmppng,tmpjpg)
    Painted = Paint(tmpjpg)
    Painted.load()
    
    if epsn == '' and figid != None:
        try: close(figid)
        except : pass
        show._needmain = False
        figure(figid) ; show()
        imshow(Painted.im,origin='lower')
        
        return None
    
    Painted.save(tmppng)
    Convert(tmppng,epsn)
    os.system('rm %s %s' % (tmpjpg,tmppng))
    
    return None
        
        
def ClLatex(self,figlist,plotlist):
    """Does a pdf file with graphical output. Specific to this task."""
    # IMPORT STUFF
    from Moments.latex import LaTeX
    from graphs_local import Paint
    import string
    # END IMPORT
    # INPUTS
    object = self['name']
    try: tag = self.execpars['tag'][0]
    except : tag = ''
    band = '$%s$'  % tag
    band = string.replace(band,'_','\_')
    try: z = self['MCz']
    except : z = -99.0
    # END INPUTS
    
    # Table with Images:
    figtxt = {} 
    for fign in figlist:
        if fign != None:
            P = Paint(fign) ; P.load()
            imsize = P.size ; del P
            if imsize[1]>=imsize[0]: dimkey = 'width'
            else : dimkey = 'height'
            figtxt[fign] = '\includegraphics[%s=5cm]{%s}' % (dimkey,fign)
        else:
            figtxt[fign] = ''
    
    Fig_Latex = ['\\begin{longtable}{|c|c|}',\
    '\hline']
    
    cols = ''
    for i in range(len(figlist)):
        if i % 2 == 0: cols += '%s &' % figtxt[figlist[i]]
        else: cols += ' %s\\\\ \hline' % figtxt[figlist[i]]
    if len(figlist) % 2 != 0: cols += ' \\\\ \hline'
    
    Fig_Latex.append(cols)
    Fig_Latex.append('\end{longtable}')
    
    # Table with Plots:
    plotxt = {} 
    for plotn in plotlist:
        if plotn != None:
            plotxt[plotn] = '\includegraphics[width=8cm]{%s}' % plotn
        else:
            plotxt[plotn] = ''
    
    Plot_Latex = ['\\begin{longtable}{|c|c|}',\
    '\hline']
    
    cols = ''
    for i in range(len(plotlist)):
        if i % 2 == 0: cols += '%s &' % plotxt[plotlist[i]]
        else: cols += ' %s\\\\ \hline' % plotxt[plotlist[i]]
    if len(plotlist) % 2 != 0: cols += ' \\\\ \hline'
    
    Plot_Latex.append(cols)
    Plot_Latex.append('\end{longtable}')
    
    # Numeric Results:
    scale = self.execpars['scale'][0]
    magzero = self.execpars['magzero'][0]
    
    def domag(flux):
        if isinstance(flux,float):
            if flux > 0. : return magzero - 2.5 * num.log10(flux)
            else : return -99.0
        else:
            pos = num.where(flux > 0.)
            mag = flux * 0. - 99.0
            if len(pos[0])>0 : mag[pos] = magzero -2.5 * num.log10(flux[pos])
            return mag
    def doarc2(pixels):
        if isinstance(pixels,float):
            if pixels > 0: return num.log10(pixels * scale**2.)
            else : return -99.0
        else:
            pos = num.where(pixels > 0.)
            arc2 = pixels * 0. - 99.0
            if len(pos[0])>0: arc2[pos] = num.log10(pixels[pos] * scale**2.)
            return arc2
    
    Results_Latex = ['\section{Results}',
    'Ncl = %i\\\\' % self['Ncl'],
    'Fcl(mag) = %.2f phiFcl = %.2f\\\\' % (domag(self['Fcl']),\
    self['phiFcl']),
    'Acl(log(arcs2)) = %.2f phiAcl = %.2f\\\\' % (doarc2(self['Acl']),\
    self['phiAcl']),
    'maxFcl minFcl meanFcl medianFcl sigmaFcl p25Fcl p75Fcl (mag)\\\\',
    ('%.2f '*7+'\\\\') % tuple(domag(num.array([self['maxFcl'],
    self['minFcl'],self['meanFcl'],self['medianFcl'],self['sigmaFcl'],
    self['p25Fcl'],self['p75Fcl']]))),
    'maxAcl minAcl meanAcl medianAcl sigmaAcl p25Acl p75Acl (log(arcs2))\\\\',
    ('%.2f '*7+'\\\\') % tuple(doarc2(num.array([self['maxAcl'],
    self['minAcl'],self['meanAcl'],self['medianAcl'],self['sigmaAcl'],
    self['p25Acl'],self['p75Acl']])))]
    
    # PDF
    
    pdfid = 'CL_%s%s' % (object,tag)
    latexfile = '%s.tex' % pdfid
    pdffile = '%s.pdf' % pdfid
    psfile = '%s.ps' % pdfid
    
    latex = LaTeX()
    header = 'Object: %s z=%s, band: %s\\\\' % (object,z,band)
    latex.body.append(header)
    
    for line in Fig_Latex:
        latex.body.append(line)
    
    for line in Plot_Latex:
        latex.body.append(line)
    
    for line in Results_Latex:
        latex.body.append(line)
        
    figs2erase = {}
    
    for i in range(len(figlist)): 
        if figlist[i] != None: figs2erase[i] = figlist[i]
    for j in range(len(plotlist)) : 
        if plotlist[j] != None: figs2erase[i+j+1] = plotlist[j]
    
    latex.Write(latexfile)
    latex.Compile(latexfile,cleanafter=True,\
    figures=figs2erase)
    latex.Ps2Pdf(psfile,pdffile,cleanafter=True)

    self['pdf'] = pdffile
