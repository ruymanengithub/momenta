#! /usr/bin/env python

from pdb import set_trace as stop
import numpy as num
from pylab import show,plot,imshow,xlabel,ylabel,figure,close,ylim,xlim,\
axhline,axvline,connect,disconnect,text,errorbar,savefig,title
from scipy import interpolate as interpol
from scipy import polyfit, polyval
from time import time
import os

def trunc_man(self):
    """Gets truncation radius in spiral galaxies... In colaboration with
    Nacho Trujillo. Manual version."""
    # IMPORT STUFF
    from deriv import deriv
    from momsource import momsource
    from Moments.graphs_local import Paint, doStamp, Convert
    from Moments.latex import LaTeX
    from copy import copy
    from radial import back_radial
    from CARSMath import polyfitw
    from LadFit import LadFit
    from time import time
    import string, os
    from trunc_man import clicker,getManTrunc
    from trunc_man import TRadialPlot,DrawOnGal,EllipseOnDs9
    from flags_truncation import addflag,allflags
    from trunc_lib import show_on_ds9, save_img_ds9, get_coo_imexam
    # END IMPORT
    
    isthere = os.path.exists
    
    # INPUTS
    img = self['STAMP'].copy()
    objectmask = self['MASK'].copy()
    othermask = self['MASKOTHER'].copy()
    skysigma = self.execpars['sigma_sky'][0]
    object = self['name']
    z = self['MCz']
    scale = self.execpars['scale'][0]
    zeropoint = self.execpars['magzero'][0]
    tag = self.execpars['tag'][0]
    dograph = self.execpars['dograph'][0] == 1
    saveradial = self.execpars['saveradial'][0] == 1
    try : commfile = self.execpars['commfile'][0]
    except NameError : commfile = False
    timetag = '%f' % time()
    # END INPUTS
    
    band = string.replace(tag,'_','')
    imgshape = img.shape
    Title = '%s %s %s' % (object,z,band)
    
    dummy = copy(self) # we need to modify some parameters in self
    
    # get background (1st way: median on sky area)
    
    back1 = self['SKY_MEDIAN']
    dummy['BACKGROUND'] = back1
    
    # get elipsoidal parameters
    
    dummy.getmoments()
    dummy.ellipse()
    
    getmompars = ['M_X','M_Y','M_X2','M_Y2','M_XY']
    for par in getmompars : self[par] = copy(dummy[par])
    
    ellipsepars=['M_A','M_B','M_THETA','M_ELONG','M_ELLIP']
    for par in ellipsepars : self[par] = copy(dummy[par])
    
    # Some inputs are changed... elipsoidal apertures computed by the
    # program itself, not provided by sextractor.
    
    dummy['X_IMAGE'] = dummy['M_X']
    dummy['Y_IMAGE'] = dummy['M_Y']
    dummy['ELLIPTICITY'] = dummy['M_ELLIP']
    dummy['THETA_IMAGE'] = dummy['M_THETA']
    
    newcoo = (-99.0,-99.0)
    Tx,Ty = newcoo
    
    # GET RADIAL PROFILE (INTERACTIVE)
    
    doneRadial = False
    while not doneRadial :
     
        # get radial profile.
     
        dummy.execpars['doLogRadial'] = [0] # Linear increase in radius
        dummy.execpars['RadialStep'] = [1.0] # 1 pix step. 
     
        dummy.radial_v3(MaxBadRatio=0.25,ExtrapolBadPix=False,\
        IgnoreBadValues=True,dograph=False)
     
        # extract "previous" radii, intensities and errors in intensities.
        radialp = dummy['M_RADIAL']
        intensp = radialp['intens']
        radiip = radialp['radii']
        del radialp
        
        ellipgal = dummy['M_ELLIP']
        xgal = dummy['X_IMAGE'] - dummy['MXMIN_IMAGE']
        ygal = dummy['Y_IMAGE'] - dummy['MYMIN_IMAGE']
        thetagal = dummy['M_THETA']
     
        # get background (2nd way: on radial profile)
     
        areaobject = float(len(num.where(objectmask==0)[0]))
        skyradius = (areaobject / (num.pi * (1.-ellipgal)))**0.5
        skyradius *= 1.5 # skyradius a 50% bigger than object radius, + or -
     
        back2 = back_radial(intensp,radiip,skyradius)
     
        # subsctract background again if necessary
    
        if back2 != None:
            backtotal = back1+back2
            dummy['BACKGROUND'] = backtotal
            
            ppcumulflx = dummy['M_RADIAL']['cumulflx'].copy()
            ppnpix = dummy['M_RADIAL']['npix'].copy()
            ppnpixout = dummy['M_RADIAL']['npixout'].copy()
            ppintens = dummy['M_RADIAL']['intens'].copy()
            
            ppcumulflx = ppcumulflx - back2 * (ppnpix-ppnpixout)
            ppintens = ppintens - back2
            
            dummy['M_RADIAL']['intens'] = ppintens.copy()
            dummy['M_RADIAL']['cumulflx'] = ppcumulflx.copy()
            del ppcumulflx,ppintens,ppnpix,ppnpixout
            
            ##dummy.radial_v3(MaxBadRatio=0.25,ExtrapolBadPix=False,\
            ##IgnoreBadValues=True,dograph=False)
        else:
            backtotal = back1
            back2 = 0.
    
    
        # SAVE RADIAL PROFILE
        if saveradial:
            radialf = '%s_RADTRUNC%s.txt' % (object,tag)
            if isthere(radialf) : os.system('rm %s' % radialf)
            dummy.SaveRadial(radialf)
        else : radialf = 'None'
    
        # Extract final radii, intensities and errors in intensities
    
        radial = dummy['M_RADIAL']
        intens = radial['intens']
        eintens = radial['eintens']
        radii = radial['radii']
        del radial
     
        # GET log(I)
        
        positive = num.where(intens > 0.)
        intens = intens[positive]
        eintens = eintens[positive]
        radii = radii[positive]
        
        
     
        logI = num.log10(intens)
        mu = zeropoint - 2.5 * ( logI + num.log10(1./scale**2.) )
        emu = 2.5 * num.log10(1.+eintens/intens)
     
        # Local scale
     
        h =  -1.086 /( deriv(radii,mu)) 
     
        # show image with elliptical aperture equal to sky radius.
     
        tmppng = 'tmpTrunc_%s.png' % timetag
        tmpjpg = 'tmpTrunc_%s.jpg' % timetag
        if isthere(tmppng) : os.system('rm %s' % tmppng)
        if isthere(tmpjpg) : os.system('rm %s' % tmpjpg)
        stampimg = img-backtotal
        stampimg[num.where(othermask == 1)] = 0.
        #doStamp(stampimg,tmppng,palette='rainbow',\
        #minhisto=20,maxhisto=99.5,format='PNG')
        stampformat='PNG'
        stamppars = '--palette=rainbow -l --min=%.e --max=%.e' % \
        (0.,stampimg.max())
        doStamp(stampimg,tmppng,format=stampformat,
        parameters=  stamppars)
        Convert(tmppng,tmpjpg)
        Painted = Paint(tmpjpg)
        qimg = 1 - ellipgal ; paimg = thetagal
        centerimg = (xgal,ygal)
        ellippar = {'qimg':qimg,'paimg':paimg,'centerimg':centerimg}
        
        close(1)
        show._needmain = False
        figure(1) ; show()
        ellippar['radimg'] = skyradius
        DrawOnGal(Painted,ellippar,color='red')
        imshow(Painted.im,origin='lower')
     
        # show intensity profile (radius in arcsecs, Intensity in mag/arcsec**2.)
        close(2)
        show._needmain = False
        figure(2) ; show()
    
        yplot = mu
        yerr = emu
        xplot = radii * scale
        xlimval = (0.,min(skyradius*scale,xplot.max()))
        ylimval = (29.0,min(19.5,yplot.min()))
        
        TRadialPlot(xplot,yplot,yerr,xlimval,ylimval,plotdata=None)
        title(Title)
        
        # Show image and mask on ds9
        imagesds9 = (stampimg,othermask.astype('Float32'))
        imgds9n = 'DS9_%s_img.fits' % object
        mskds9n = 'DS9_%s_msk.fits' % object
        namesds9 = (imgds9n,mskds9n)
        if isthere(imgds9n) : os.system('rm %s %s' % namesds9)
        myds9 = show_on_ds9(imagesds9,namesds9)
        
        # Plot sky radius aperture on ds9
        EllipseOnDs9(myds9,xgal,ygal,skyradius,skyradius*qimg,\
        thetagal,color='white')
        #os.system('echo "ellipse %f %f %f %f %f" | xpaset ds9 regions' %
        #(xgal,ygal,skyradius,skyradius*qimg,thetagal))
        
                # RECENTER?
        
        doCenter = raw_input('Wanna recenter? yes = y ')
        try : doCenter = string.lower(doCenter)[0] == 'y'
        except IndexError: doCenter = True
        
        if doCenter:
            newcoo = get_coo_imexam()
            if newcoo[0] > 0.:
                Tx = newcoo[0] + dummy['MXMIN_IMAGE']
                Ty = newcoo[1] + dummy['MYMIN_IMAGE']
                dummy['X_IMAGE'] = Tx ; dummy['Y_IMAGE'] = Ty
        else:
            dummy['X_IMAGE'] = dummy['M_X']
            dummy['Y_IMAGE'] = dummy['M_Y']
            # Tx = -99.0 ; Ty = -99.0 ; newcoo = Tx,Ty
        
        if not doCenter:
            doneRadial = raw_input('Are you done with the radial profile? yes = y ')
            try: doneRadial = string.lower(doneRadial)[0] == 'y'
            except IndexError: doneRadial = False
        else: doneRadial = False
        
    del dummy
    
    # GET TRUNCATION RADIUS (INTERACTIVE)
    
    TruncTypes = ['I','II','III']
    
    satisfied = False
    while not satisfied:
        
        # is there truncation?
        TruncType = 'I'
        istrunc = raw_input('Is there any truncation? Yes = y ')
        istrunc = string.lower(string.strip(istrunc))[0] == 'y'
        
        if istrunc:
            ttype = raw_input('Which truncation type: 2 or 3? ')
            TruncType = TruncTypes[int(ttype)-1]
            print '\nYou think this is type %s\n' % TruncType
        
        figure(2)
        
        #   Pinpoint Truncation Radius and Truncation mu
        
        if istrunc:
            doneMeasure = False
            while not doneMeasure:
                print '\n Measure of Truncation 1\n'
                Tr, Tmu = getManTrunc(xplot,yplot,xlimval,ylimval,figid=2) # Tr in arcseconds!
                redo = raw_input('if want to mark Tr again, press "y" ')
                if string.lower(redo) != 'y' : doneMeasure = True
            
        # Show Trunc radius on intensity plot.
        
        show._needmain = False
        figure(2) ; show()
        
        if istrunc:
            Trplot = Tr ; Tmuplot = Tmu 
        else:
            Trplot = 0. ; Tmuplot = 0. ;
         
        plotdata = {'Tr':Trplot,'Tmu':Tmuplot,'TruncType':TruncType,
        'Title':Title}
        TRadialPlot(xplot,yplot,yerr,xlimval,ylimval,plotdata)
        
        # Show trunc radius as ellipse on image.
            
        close(1)
        show._needmain = False
        figure(1) ; show()
        
        if istrunc: ellippar['radimg'] = Tr/scale
        else : ellippar['radimg'] = 0.
        
        DrawOnGal(Painted,ellippar,color='red')
        imshow(Painted.im,origin='lower')
        
        # Plot Truncation radius on ds9
        myds9.xpaset('regions delete all')
        EllipseOnDs9(myds9,ellippar['centerimg'][0],\
        ellippar['centerimg'][1],ellippar['radimg'],ellippar['radimg']*qimg,\
        thetagal,color='white')
       
        if not istrunc:
            TruncType= 'I'
            Tmu = -99.0 ; 
            Tr = -99.0 ;
        
        # satisfied? if satisfied quit, if not, repeat.
        ans = raw_input('Are you done with analyzing this object? yes=y ')
        satisfied = string.lower(string.strip(ans))[0] == 'y'
    
    # Believable results?
    Believe = raw_input('Do you deem these results as correct? yes=y ')
    Believe = '%s' % (string.lower(string.strip(Believe))[0] == 'y')
    # Comments?
    comments = ''
    comments = raw_input('Any declaration for the press? ')
    
    close(1) ; close(2)
    
    # Flagging
    # NO FLAGGING FOR MANUAL MODE
    
    if dograph:
     
        # GRAPHIC OUTPUTS (pdffile made with latex):
     
        # image with elipsoidal aperture on Truncation radius : png
        pngGal = '%s_IMGTRUNC%s.png' % (object,tag)
        epsGal = '%s_IMGTRUNC%s.eps' % (object,tag)
        
        readyimg = False
        while not readyimg:
            save_img_ds9(myds9,pngGal,format='jpeg')
            Convert(pngGal,epsGal)
            os.system('gv %s &' % epsGal)
            readyimg = raw_input('Wanna repeat image capture? yes = y ')
            readyimg = string.strip(string.lower(readyimg))[0] != 'y'
            if not readyimg : os.system('rm %s %s' % (pngGal,epsGal))
            else: os.system('rm %s' % pngGal)
    
        # radial profile with truncation radius : eps
        epsRad = '%s_FIGTRUNC%s.eps' % (object,tag)
        figure(4)
        plotdata['Title'] = ''
        TRadialPlot(xplot,yplot,yerr,xlimval,ylimval,plotdata)
        savefig(epsRad)
        close(4)
     
        fig1 = epsGal
        if imgshape[1]>=imgshape[0]: dimkey1 = 'width'
        else : dimkey1 = 'height'
        
        fig2 = epsRad
        dimkey2 = 'width'
        
        # table with figures:
        Fig_Latex = [
        '\\begin{longtable}{|c|}',\
        '\hline',\
        '\includegraphics[%s=5cm]{%s}\\\\' % (dimkey1,fig1),\
        '\hline',\
        '\includegraphics[%s=8cm]{%s}\\\\' % (dimkey2,fig2),\
        '\hline',\
        '\end{longtable}']
        
        # Numeric Results:
        Results_Latex = [
        '\section{Results}',\
        'I trust these Results: %s\\\\' % Believe,\
        'Trunctype : %s\\\\' % TruncType,\
        'Tr = %.2f"\\\\'% (Tr,),\
        'Tmu = %.2f mag/arcsec2\\\\' % (Tmu,),
        'Comments : %s\\' % comments]
        
        # PDF
        
        pdfid = '%s_TRUNCPDF%s' % (object,tag)
        latexfile = '%s.tex' % pdfid
        pdffile = '%s.pdf' % pdfid
        psfile = '%s.ps' % pdfid
        
        latex = LaTeX()
        header = 'Object: %s z=%s, band: %s\\\\' % (object,z,band)
        latex.body.append(header)
        
        for line in Fig_Latex:
            latex.body.append(line)
        
        for line in Results_Latex:
            latex.body.append(line)
        
        figs2erase = {1:fig1,2:fig2}
        latex.Write(latexfile)
        latex.Compile(latexfile,cleanafter=True,\
        figures=figs2erase)
        latex.Ps2Pdf(psfile,pdffile,cleanafter=False)
        
    else:
        pdffile = 'None'
    
    if commfile:
        f = open(commfile,'a')
        commline = '%s   %s' % (object,comments)
        print >> f,commline
        f.close()
    
    # OUTPUTS
    
    skyradius_save = skyradius * scale
    
    outdata = {'TType':TruncType,'Tx':Tx,'Ty':Ty,'Tback1':back1,
    'Tskyradius':skyradius_save,
    'Tback2':back2,'Tmu':Tmu,'Tr':Tr,
    'Tscale':scale,'Tzero':zeropoint,'TBelieve':Believe,
    'Tradialf':radialf,'Tpdf':pdffile}
    
    for key in outdata: self[key] = copy(outdata[key])
    
    if isthere(imgds9n) : os.system('rm %s %s' % namesds9)
    os.system('rm %s %s' % (tmppng,tmpjpg))
    
    try:
        isopen = myds9.isOpen()
        if isopen : myds9.xpaset('exit')
    except NameError: pass
    
    return None

class clicker(dict):
    """Class to be used with pylab.connect to grab coordinates interactively from a display with
    the mouse."""
    
    def __init__(self,nclicks=1):
        self.data = []
        self.nclicks = nclicks
        self.nevents = 0

    def get_clicks(self,event):
        # get the x and y pixel coords
        x, y = event.xdata, event.ydata
        try:
            xlimval = self['xlimval']
            ylimval = self['ylimval']
        except KeyError: pass
        
        if event.inaxes:
            self.nevents += 1
            print '\ncoords = %.1f, %.1f' % (x,y)
            #axhline(y,color='k')
            #axvline(x,color='k')
            plot([x],[y],'kx',ms=20,mew=2)
            try:
                xlim(xlimval)
                ylim(ylimval)
            except NameError: pass
            
            if len(self.data) < self.nclicks :
                self.data.append([x,y])
            elif len(self.data) == nclicks:
                index = self.nevents % self.nclicks - 1 
                self.data[index] = [x,y]

def getManTrunc(xplot,yplot,xlimval,ylimval,figid=1):
    # IMPORT STUFF
    # END IMPORT
    
    click = clicker(2) # 2 is the number of events to register
    click['xlimval'] = xlimval ; click['ylimval'] = ylimval
    event = connect('button_press_event',click.get_clicks)
    ready = raw_input('When you have marked the truncation, kick any key ')
    if ready: disconnect(event)
    try: 
        Tr = click.data[0][0] ; Tmu = click.data[0][1]
    except IndexError:
        Tr = xplot.max() ; Tmu = yplot.max()

    print '\nTr : %.2f' % (Tr,)
    print '\nTmu : %.2f' % (Tmu,)
    
    # Draw selected Truncation
    close(figid)
    show._needmain = False
    figure(figid) ; show()
    plot(xplot,yplot,'b.')
    axvline(Tr,color='k',linestyle='-') # Show Truncation Radius
    axhline(Tmu,color='k',linestyle='-') # Show Truncation mu
    #ylim(num.ceil(yplot.max()),num.floor(yplot.min()))
    xlim(xlimval) ; ylim(ylimval)
    xlabel('arcsec') ; ylabel('mag/arcsec^2')

    return Tr,Tmu

    
    
def TRadialPlot(xplot,yplot,yerr,xlimval,ylimval,plotdata=None):
    """Does a Radial Plot... specific to this task."""
    
    errorbar(xplot,yplot,yerr=yerr,xerr=None,fmt='b.',ecolor=None)
    
    try:
        if plotdata != None:
            Tr = plotdata['Tr']
            Tmu = plotdata['Tmu']
            TruncType = plotdata['TruncType']
            Title = plotdata['Title']
        
        text(Tr*1.2,Tmu-1.5,'Type %s' % TruncType)
        text(Tr*1.2,Tmu-1.,'T. Radius = %.2f"' % Tr)
        text(Tr*1.2,Tmu-0.5,'T. mu = %.2f mag/arcsec^2' % Tmu)
        if Tr > 0.:
            axvline(Tr,color='k',linestyle='-') # Show Truncation Radius
            axhline(Tmu,color='k',linestyle='-')
        if Title != '': title(Title)
    
    except NameError : pass 
    ylim(ylimval)
    xlim(xlimval)
    xlabel('arcsec') ; ylabel('mag/arcsec^2')

def DrawOnGal(Painted,ellippar=None,color='red'):
    """shows an image of the galaxy and draws an ellipse on it. Specific to 
    this task"""
    
    Painted.load()

    if ellippar != None:
        
        try:
            qimg = ellippar['qimg'] ; paimg = ellippar['paimg']
            centerimg = ellippar['centerimg'] ; radimg = ellippar['radimg']
            if radimg > 0.:
                Painted.DrawEllipse(centerimg,radimg,qimg,paimg,\
                color=color,linewidth=2)
        except KeyError,NameError: pass
    
def EllipseOnDs9(myds9,x,y,major,minor,angle,color='white'):
    """Draws an ellipse on a ds9 display."""
    if major == 0. : return None
    
    text = \
"""# Region file format: DS9 version 4.0 
# Filename : ds9.fits
global color=%s font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source 
image
ellipse(%f,%f,%f,%f,%f)""" % (color,x,y,major,minor,angle)
    
    tmpreg = 'tmpDS9ellipse%f.reg' % time()
    f = open(tmpreg,'w')
    f.write(text)
    f.close()
    myds9.xpaset('regions load %s' % tmpreg)
    os.system('rm %s' % tmpreg)
    return None
    
