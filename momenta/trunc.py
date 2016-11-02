#! /usr/bin/env python

from pdb import set_trace as stop
import numpy as num
from pylab import show,plot,imshow,xlabel,ylabel,figure,close,ylim,xlim,\
axhline,axvline,connect,disconnect,text,errorbar,savefig,title
from scipy import interpolate as interpol
from scipy import polyfit, polyval,stats
from time import time
import os

def trunc(self):
    """Gets truncation radius in spiral galaxies... In colaboration with
    Nacho Trujillo."""
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
    from trunc import clicker, selreg, get_cutTr,get_EDTr
    from trunc import TRadialPlot,DrawOnGal,EllipseOnDs9
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
    
    # GET RADIAL PROFILE (NOW INTERACTIVE)
    
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
            
            dummy['M_RADIAL']['cumulflx'] = ppcumulflx.copy()
            dummy['M_RADIAL']['intens'] = ppintens.copy()
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
        xlimval = (0.,min(skyradius*scale*1.2,xplot.max()))
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
        else : doneRadial = False
    
    del dummy
    
    # GET TRUNCATION RADIUS (INTERACTIVE)
    
    satisfied = False
    while not satisfied:
        # use 'Robust' fitting?
        DoRobust = raw_input('Use Robust fitting? no = n ')
        DoRobust = string.lower(string.strip(DoRobust)) != 'n'
        
        # is there truncation?
        TruncType = 'I'
        istrunc = raw_input('Is there any truncation? Yes = y ')
        try: istrunc = string.lower(string.strip(istrunc))[0] == 'y'
        except IndexError: istrunc = False
    
        figure(2)
        
        #   mark inner exponential region (2 events). show selected points.
        
        doneReg1 = False
        while not doneReg1:
             
            print '\n SELECTION OF REGION 1\n'
            b1prior = [xplot.min(),yplot.min()]
            b2prior = [(xplot.max()-xplot.min())/2.,yplot.max()]
            b1x, b1y,b2x,b2y,reg1 = selreg(xplot,yplot,b1prior,b2prior,\
            xlimval,ylimval,figid=2)
            redo1 = raw_input('if want to repeat sel. Reg. 1, press "y" ')
            if string.lower(redo1) != 'y' : doneReg1 = True
            if len(reg1[0]) < 3: 
                print 'At least pick 3 points'
                doneReg1 = False
            
        # mark outer exponential region (2 events). show selected points.
        if istrunc:
         
            doneReg2 = False
            while not doneReg2:
             
                print '\n SELECTION OF REGION 2\n'
                b3prior = [(xplot.max()-xplot.min()) / 2.,yplot.min()]
                b4prior = [xplot.max(),yplot.max()]
                b3x, b3y, b4x, b4y,reg2 = selreg(xplot,yplot,b3prior,b4prior,\
                xlimval,ylimval,figid=2)
                redo2 = raw_input('if want to repeat sel. Reg. 2, press "y" ')
                if string.lower(redo2) != 'y' : doneReg2 = True
                if len(reg1[0]) < 3: 
                    print 'At least pick 3 points'
                    doneReg1 = False
        
        # fit to lines both sets of points: 
        # I use either CARSMath.polyfitw or LadFit 
        # (a translation of idl's ladfit "robust" absolute
        # deviance fitting)
        
        if not DoRobust:
            weights1 = 1./emu[reg1]**2.
            coeffs1,yfit1,yband1,sigma1,a1 = polyfitw(xplot[reg1],\
            yplot[reg1],weights1,ndegree=1,return_fit=1)
            ec1 = a1[0,0]**0.5 ; em1 = a1[1,1]**0.5 #
        else:
            AD1 = [0.] ; sigma1 = [0.]
            coeffs1 = LadFit(xplot[reg1],yplot[reg1],AD1,sigma1)
            sigma1 = sigma1[0]
            ec1 = 0. ; em1 = 0.
        
        c1, m1 = coeffs1
        r1 = stats.pearsonr(xplot[reg1],yplot[reg1])[0]
        h1 = 1.086/m1 ; 
        eh1 = em1 * h1/m1 ; 
        
        print '\nc1=%.2f+-%.2f, m1=%.2f+-%.2f, r1=%.3f' % \
        (c1,ec1,m1,em1,r1)    
        
        if istrunc:
            if not DoRobust:
                weights2 = 1./emu[reg2]**2.
                coeffs2,yfit2,yband2,sigma2,a2 = polyfitw(xplot[reg2],yplot[reg2],\
                weights2,ndegree=1,return_fit=1)
                ec2 = a2[0,0]*0.5 ; em2 = a2[1,1]**0.5 #
            else :
                AD2 = [0.] ; sigma2 = [0.]
                coeffs2 = tuple(LadFit(xplot[reg2],yplot[reg2],AD2,sigma2))
                sigma2 = sigma2[0]
                ec2 = 0. ; em2 = 0.
            
            c2, m2 = coeffs2
            r2 = stats.pearsonr(xplot[reg2],yplot[reg2])[0]
            print 'c2=%.2f+-%.2f, m2=%.2f+-%.2f, r2=%.3f' % \
            (c2,ec2,m2,em2,r2)
            
            # Compute Trunc Radius, Intensity (min abs dev method)
            
            DoEqDev = raw_input('\nDo you want to try EqDev method? yes=y ')
            DoEqDev = string.lower(string.strip(DoEqDev))[0] == 'y'
            
            if DoEqDev:
                doneEDreg = False
                while not doneEDreg:
                    b1EDprior = [b1x,b1y]
                    b2EDprior = [b4x,b4y]
                    b1EDx, b1EDy, b2EDx, b2EDy,EDreg = \
                    selreg(xplot,yplot,b1EDprior,b2EDprior,xlimval,ylimval,figid=2)
                    if len(EDreg[0]) < 3: 
                        print 'At least pick 3 points'
                    else: doneEDreg = True
                    
                Tr,eTr,Tmufit,eTmufit = \
                get_EDTr(xplot,yplot,EDreg,coeffs1,coeffs2,sigma1,sigma2)
                if -99.0 in (Tr,eTr,Tmufit,eTmufit):
                    istrunc= False
                    print """\nSomething went wrong with finding the truncation...\n"""
            else:
                # Compute Trunc Radius, Intensity, h1, h2. (cut method)
                Tr,eTr,Tmufit,eTmufit = get_cutTr(coeffs1,coeffs2,sigma1,sigma2)
            
            h2 = 1.086/m2
            eh2 = em2 * h2/m2
            
            if h1 > h2 : TruncType = 'II'
            else: TruncType= 'III'
            
            # Now the Truncation intensity on the intensity curve by
            # means of interpolation
            
            if istrunc:
                Tspl = interpol.splrep(xplot,yplot,s=0)
                Tmu = interpol.splev(Tr,Tspl)
            else: Tmu = -99.0
            
            print 'Truncation Type : %s' % TruncType
            print 'Radius = %.2f +- %.2f' % (Tr,eTr)
            print 'mufit = %.2f +- %.2f' % (Tmufit,eTmufit)
            print 'mu = %.2f\n' % (Tmu,)
            print 'h1 = %.2f+-%.2f", h2=%.2f+-%.2f"' % (h1,eh1,h2,eh2)
        
        # Show fits and trunc radius on intensity plot.
        
        try: xfit1 = num.array([b1x*0.8,max(Tr*1.2,b2x*1.2)])
        except UnboundLocalError: xfit1 = num.array([b1x*0.8,b2x*1.2])
        yfit1 = polyval(coeffs1[::-1],xfit1)
        
        if istrunc:
            xfit2 = num.array([min(Tr*0.8,b3x*0.8),b4x*1.2])
            yfit2 = polyval(coeffs2[::-1],xfit2)
        else:
            xfit2 = copy(xfit1) ; yfit2 = copy(yfit1)
        
        close(2)
        show._needmain = False
        figure(2) ; show()
        
        if istrunc:
            Trplot = Tr ; Tmuplot = Tmu ; h2plot = h2
            Tmufitplot = Tmufit
        else:
            Trplot = 0. ; Tmuplot = (yfit1.max()-yfit1.min())/2. ;
            Tmufitplot = Tmuplot
            h2plot = -99.0
         
        plotdata = {'xfit1':xfit1,'yfit1':yfit1,'xfit2':xfit2,'yfit2':yfit2,
        'Tr':Trplot,'Tmu':Tmuplot,'Tmufit':Tmufitplot,'TruncType':TruncType,
        'h1':h1,'h2':h2plot,'Title':Title}
        TRadialPlot(xplot,yplot,yerr,xlimval,ylimval,plotdata)
            
        # Show trunc radius as ellipse on image.
            
        close(1)
        show._needmain = False
        figure(1) ; show()
        
        if istrunc: ellippar['radimg'] = Tr/scale
        else : ellippar['radimg'] = h1/scale
        
        DrawOnGal(Painted,ellippar,color='red')
        imshow(Painted.im,origin='lower')
        
        # Plot Truncation radius on ds9
        myds9.xpaset('regions delete all')
        EllipseOnDs9(myds9,ellippar['centerimg'][0],\
        ellippar['centerimg'][1],ellippar['radimg'],ellippar['radimg']*qimg,\
        thetagal,color='white')
        #os.system('echo "ellipse %f %f %f %f %f" | xpaset ds9 regions' %
        #(ellippar['centerimg'][0],ellippar['centerimg'][1],
        #ellippar['radimg'],ellippar['radimg']*qimg,thetagal))
        
        if not istrunc:
            TruncType= 'I'
            b3x = -99.0 ; b4x = -99.0
            b3y = -99.0 ; b4y = -99.0
            m2 = -99.0 ; em2 = -99.0
            c2 = -99.0 ; ec2 = -99.0
            h2 = -99.0 ; eh2 = -99.0
            r2 = -99.0
            Tmu = -99.0 ; Tmufit = -99.0 ; eTmufit = -99.0
            Tr = -99.0 ; eTr = -99.0
        #
        
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
    if DoRobust:
        self['flags'] = addflag(self['flags'],allflags['ROBUST'])
    if istrunc:
        if DoEqDev:
            self['flags'] = addflag(self['flags'],allflags['EDMETHOD'])
        else:
            self['flags'] = addflag(self['flags'],allflags['CUTMETHOD'])
     
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
        
        #figure(3)
        # ellippar['radimg'] = Tr/scale
        #DrawOnGal(Painted,ellippar,color='red')
        #Painted.save(pngGal)
        #close(3)
        #Convert(pngGal,epsGal)
        #os.system('rm %s' % pngGal)
    
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
        'b1x = %.1f  b2x = %.1f\\\\' % (b1x,b2x),\
        'b3x = %.1f  b4x = %.1f\\\\' % (b3x,b4x),\
        'm1 = %.2f +- %.2f   c1 = %.2f +- %.2f r1 = %.3f\\\\' %\
        (m1,em1,c1,ec1,r1),\
        'm2 = %.2f +- %.2f   c2 = %.2f +- %.2f r2 = %.3f\\\\' %\
        (m2,em2,c2,ec2,r2),\
        'h1 = %.2f +- %.2f"   h2 = %.2f +- %.2f"\\\\' %\
        (h1,eh1,h2,eh2),\
        'Tr = %.2f +- %.2f"\\\\'% (Tr,eTr),\
        'Tmufit = %.2f +- %.2f mag/arcsec2\\\\' % (Tmufit,eTmufit),\
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
    'Tback2':back2,'Tb1x':b1x,'Tb1y':b1y,'Tb2x':b2x,'Tb2y':b2y,'Tb3x':b3x,
    'Tb3y':b3y,'Tb4x':b4x,'Tb4y':b4y,'Tm1':m1,'eTm1':em1,'Tm2':m2,
    'eTm2':em2,'Tc1':c1,'eTc1':ec1,'Tc2':c2,'eTc2':ec2,'Th1':h1,
    'eTh1':eh1,'Th2':h2,'eTh2':eh2,'Trpearson1':r1,'Trpearson2':r2,
    'Tmu':Tmu,'Tmufit':Tmufit,'eTmufit':eTmufit,'Tr':Tr,'eTr':eTr,
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

def selreg(xplot,yplot,leftprior,rightprior,xlimval,ylimval,figid=1):
    # IMPORT STUFF
    # END IMPORT
    
    click = clicker(2) # 2 is the number of events to register
    click['xlimval'] = xlimval ; click['ylimval'] = ylimval
    event = connect('button_press_event',click.get_clicks)
    ready = raw_input('when region is selected, press any key ')
    if ready: disconnect(event)
    try: 
        b1x = click.data[0][0] ; b2x = click.data[1][0]
        b1y = click.data[0][1] ; b2y = click.data[1][1]
    except IndexError:
        b1x = leftprior[0] ; b2x = rightprior[0]
        b1y = leftprior[1] ; b2y = rightprior[1]
        
    print '\nRegion : %.2f to %.2f\n' % (b1x,b2x)
    reg = num.where((xplot>=b1x) & (xplot<=b2x) & 
    (yplot>=b1y) & (yplot<=b2y))
    
    # Draw region selected
    close(figid)
    show._needmain = False
    figure(figid) ; show()
    plot(xplot,yplot,'b.')
    plot(xplot[reg],yplot[reg],'r+')
    #ylim(num.ceil(yplot.max()),num.floor(yplot.min()))
    xlim(xlimval) ; ylim(ylimval)
    xlabel('arcsec') ; ylabel('mag/arcsec^2')

    return b1x,b1y,b2x,b2y,reg

def get_cutTr(coeffs1,coeffs2,sigma1=0.,sigma2=0.):
        """Retrieves truncation radius, and intensity, given both fits, as the
        cut point of both lines"""
        def get_cut(m1,m2,c1,c2):
            cutx = - (c1-c2) / (m1-m2)
            cuty = m1 * cutx + c1
            return cutx, cuty
        m1,m2,c1,c2 = coeffs1[1],coeffs2[1],coeffs1[0],coeffs2[0]
        Tr,Tmufit = get_cut(m1,m2,c1,c2)
        if (sigma1 != 0) and (sigma2 != 0.):
            veTr = [] ; veTmufit = []
            for i in range(-1,2,2):
                for j in range(-1,2,2):
                    dTr, dTmufit = get_cut(m1,m2,c1+i*sigma1,c2+j*sigma2)
                    veTr.append(dTr)
                    veTmufit.append(dTmufit)
            veTr = num.array(veTr) ; veTmufit = num.array(veTmufit)
            eTr = (veTr.max()-veTr.min()) / 2.
            eTmufit = (veTmufit.max()-veTmufit.min()) / 2.
        else: 
            eTr = -99.0 ; eTmufit = -99.0
        return Tr,eTr,Tmufit,eTmufit
    
def find_cutsx(x,y):
        """Finds x-points at which a sequence of coordinates cuts 
        the y=0 axis. Used by get_EDTr"""
        xcut = []
        pivots = num.where((y[0:-1]*y[1:])<=0.)
        if len(pivots[0]>0):
            for pivot in pivots[0]:
                if (y[pivot] == 0.) and (y[pivot+1] != 0.): 
                    xcut.append(x[pivot])
                else:
                    y1,y2 = y[pivot],y[pivot+1]
                    x1,x2 = x[pivot],x[pivot+1]
                    m = (y2-y1) / (x2-x1)
                    b = y1 - m * x1
                    xcut.append(-b/m)
            if len(xcut) == 1 : xcut = xcut[0]
        else:
            xcut = None
        return xcut

def get_EDTr(xplot,yplot,EDreg,coeffs1,coeffs2,sigma1=0.,sigma2=0.):
        def get_ED(xplot,yplot,EDreg,coeffs1,coeffs2):   
            xED = xplot[EDreg] ; yED = yplot[EDreg]
            AD1 = num.abs(yED-polyval(coeffs1[::-1],xED))
            AD1spl = interpol.splrep(xED,AD1,s=0)
            AD2 = num.abs(yED-polyval(coeffs2[::-1],xED))
            to_solveED = AD1-AD2
            TrED = find_cutsx(xED,to_solveED)
            #print 'Checkline: TrED ',TrED
            #figure(5) ; show._needmain=False ; show()
            #plot(xED,to_solveED)
            #stop()
            #close(5)
            if TrED != None :
                if isinstance(TrED,list): 
                    AD1 = num.array(interpol.splev(TrED,AD1spl))
                    TrED = TrED[AD1.argmax()]
                Trwithfit1 = polyval(coeffs1[::-1],TrED)
                Trwithfit2 = polyval(coeffs2[::-1],TrED)
                TmuED = (Trwithfit1 + Trwithfit2)/2.
            else: TrED,TmuED = -99.,-99.
            return TrED,TmuED
        TrED,TmuED = get_ED(xplot,yplot,EDreg,coeffs1,coeffs2)
        eTrED = -99.0 ; eTmuED = -99.0
        if (sigma1 != 0) and (sigma2 != 0.):
            veTrED = [] ; veTmuED = []
            NoError = False
            for i in range(-1,2,2):
                for j in range(-1,2,2):
                    coeffs1p = [coeffs1[0]+i*sigma1,coeffs1[1]]
                    coeffs2p = [coeffs2[0]+j*sigma2,coeffs2[1]]
                    dTrED, dTmuED = get_ED(xplot,yplot,EDreg,coeffs1p,coeffs2p)
                    veTrED.append(dTrED)
                    veTmuED.append(dTmuED)
                    if (dTrED == -99.) or (dTmuED == -99.):
                       NoError = True
                       break
            if not NoError:
                veTrED = num.array(veTrED) ; veTmuED = num.array(veTmuED)
                eTrED = (veTrED.max()-veTrED.min()) / 2.
                eTmuED = (veTmuED.max()-veTmuED.min()) / 2.
        return TrED,eTrED,TmuED,eTmuED
    
def TRadialPlot(xplot,yplot,yerr,xlimval,ylimval,plotdata=None):
    """Does a Radial Plot... specific to this task."""
    
    errorbar(xplot,yplot,yerr=yerr,xerr=None,fmt='b.',ecolor=None)
    
    try:
        if plotdata != None:
            xfit1 = plotdata['xfit1']
            yfit1 = plotdata['yfit1']
            xfit2 = plotdata['xfit2']
            yfit2 = plotdata['yfit2']
            Tr = plotdata['Tr']
            Tmufit = plotdata['Tmufit']
            Tmu = plotdata['Tmu']
            TruncType = plotdata['TruncType']
            h1 = plotdata['h1']
            h2 = plotdata['h2']    
            Title = plotdata['Title']
        
        plot(xfit1,yfit1,'k-') # fit line 1
        plot(xfit2,yfit2,'k-') # fit line 2
        text(Tr*1.2,Tmu-1.5,'Type %s' % TruncType)
        text(Tr*1.2,Tmu-1.,'T. Radius = %.2f"' % Tr)
        text(Tr*1.2,Tmu-0.5,'T. mu = %.2f mag/arcsec^2' % Tmu)
        text(Tr*1.2,Tmu,'h1 = %.2f", h2 = %.2f"' % (h1,h2))
        if Tr > 0.:
            axvline(Tr,color='k',linestyle='-') # Show Truncation Radius
            axhline(Tmufit,color='k',linestyle='--') # Show Truncation mu
            axhline(Tmu,color='k',linestyle='-')
        if Title != '': title(Title)
    
    except NameError : pass 
    # ylim(num.ceil(yplot.max()),num.floor(yplot.min()))
    # xlim(num.floor(xplot.min()),num.ceil(xplot.max()))
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
            Painted.DrawEllipse(centerimg,radimg,qimg,paimg,\
            color=color,linewidth=2)
        except KeyError,NameError: pass
    
def EllipseOnDs9(myds9,x,y,major,minor,angle,color='white',width=1):
    """Draws an ellipse on a ds9 display."""
    text = \
"""# Region file format: DS9 version 4.0 
# Filename : ds9.fits
global color=%s font="helvetica 10 normal" select=1 highlite=1 edit=1  move=1 delete=1 include=1 fixed=0 width=%i source 
image
ellipse(%f,%f,%f,%f,%f)""" % (color,width,x,y,major,minor,angle)
    
    tmpreg = 'tmpDS9ellipse%f.reg' % time()
    f = open(tmpreg,'w')
    f.write(text)
    f.close()
    myds9.xpaset('regions load %s' % tmpreg)
    os.system('rm %s' % tmpreg)
    return None
    
def LineOnDs9(myds9,x1,y1,x2,y2,color='white',width=1):
    """Draws a line on a ds9 display."""
    text = \
"""# Region file format: DS9 version 4.0 
# Filename : ds9.fits
global color=%s font="helvetica 10 normal" select=1 highlite=1 edit=1  move=1 delete=1 include=1 fixed=0 width= %i source 
image
line(%f,%f,%f,%f) # line=0 0""" % (color,width,x1,y1,x2,y2)
    
    tmpreg = 'tmpDS9line%f.reg' % time()
    f = open(tmpreg,'w')
    f.write(text)
    f.close()
    myds9.xpaset('regions load %s' % tmpreg)
    os.system('rm %s' % tmpreg)
    return None
