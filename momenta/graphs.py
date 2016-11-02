#! /usr/bin/env python

# IMPORT STUFF
from pdb import set_trace as stop
import pyfits
import os
# END IMPORT

isthere = os.path.exists

fontsdir = '/net/cipres/scratch1/MYFONTS/'

def doStamp(img,imgname,palette='heat',minhisto=20,maxhisto=99.5,format='PNG'):
    """Creates a scaled stamp"""
    # IMPORT STUFF
    from time import time
    # END IMPORT
    tempfits = 'graphs_doSTAMP_%f.fits' % time()
    if isthere(tempfits) : os.system('rm %s' % tempfits)
    pyfits.writeto(tempfits,img)
    if format == 'JPEG' : format_key = 'j'
    if format == 'PNG' : format_key = 'p'
    execline = 'fitscut -%s --all --autoscale=%f,%f --palette=%s %s %s' % \
    (format_key,minhisto,maxhisto,palette,tempfits,imgname)
    os.system(execline)
    os.system('rm %s' % tempfits)
    
    return None

def Convert(imgname1,imgname2,cleanafter=False):
    """Converts formats of images"""
    os.system('convert %s %s' % (imgname1,imgname2))
    if cleanafter : os.system('rm %s' % imgname1)
    return None

def Plot(datacube,figname,sxlabel='',sylabel='',stitle='',fontsize=20,doOrigin=False):
    """Does a simple plot to an eps file"""
    # IMPORT STUFF
    from pylab import figure,close,plot,savefig,xlabel,ylabel,title
    from pylab import axhline,axvline
    import numpy as num
    import sys
    # END IMPORT
    
    if isinstance(datacube[0],tuple):
        xs = []
        ys = []
        for i in range(len(datacube)):
            xs.append(datacube[i][0])
            ys.append(datacube[i][1])
    elif isinstance(datacube[0],num.ndarray):
        xs.append(datacube[0])
        ys.append(datacube[1])
    else : stop('Execution halted at Moments.graphs')
    
    colors = ['b','g','r','c','m','k']
    vfontsize=20
    i = 0
    figure(1)
    for i in range(len(xs)):
        x = xs[i] ; y = ys[i]
        plot(x,y,'%s-.'% colors[i%len(colors)])
        i+=1
    if doOrigin:
        axhline(y=0,linewidth=1,color='k')
        axvline(x=0,linewidth=1,color='k')
    title(stitle,fontsize=vfontsize) ; xlabel(sxlabel,fontsize=vfontsize) 
    ylabel(sylabel,fontsize=vfontsize)
    savefig(figname)
    close(1)
    return None

class Paint(dict):
    """Class to draw on images"""
    import Image, ImageDraw, ImageFont
    
    def __init__(self,imgname=""):
        self.imgname = imgname
        self.im = None
        self.size = (0,0)
        self.draw = None
    
    def load(self):
        """Loads image"""
        # IMPORT STUFF
        import Image
        # END IMPORT
        self.im = Image.open(self.imgname)
        self.size = self.im.size
        self.InitDraw()
    
    def save(self,imgname,format="PNG"):
        """Saves image"""
        self.im.save(imgname,format)
    
    def InitDraw(self):
        """Initializes 'draw' attribute"""
        import ImageDraw
        self.draw = ImageDraw.Draw(self.im)
        
    def DrawEllipse(self,center,sma,q,angle,color='red',linewidth=2):
        """Draws an ellipse"""
        import numpy as num
        from copy import copy
        
        angle = angle * num.pi / 180. # to radians
        center = (center[0],self.size[1]-center[1]) # origin in upper left corner
        rads = num.arange(360) * 2. * num.pi /359.
        xp = sma * num.cos(rads)
        yp = sma * q * num.sin(rads)
        cosang = num.cos(angle)
        sinang = num.sin(angle)
        x = xp * cosang + yp * sinang
        y = -xp * sinang + yp * cosang
        x += center[0] ; y += center[1]
        xy = []
        for i in range(len(x)) : xy.append((x[i],y[i]))
        self.draw.line(xy,fill=color,width=linewidth)
        xyreverse = copy(xy)
        xyreverse.reverse()
        self.draw.line(xyreverse,fill=color,width=linewidth)
    
    def Graffiti(self,text,center,fontsize=11,color='green',\
     fontname='ActionManExtended'):
        """Writes text on image"""
        # IMPORT STUFF
        from graphs import fontsdir
        import ImageFont
        # END IMPORT
        
        size = self.size
        fullfontname = fontsdir + fontname + '.ttf'
        font = ImageFont.truetype(fullfontname,fontsize)
        textsize = self.draw.textsize(text,font)
        center = (center[0],self.size[1]-center[1]) # origin in upper left corner
        
        # TRY TO MAKE TEXT APPEAR INSIDE IMAGE
        if textsize[0] + center[0] > size[0] : 
            xcenter = max(0,center[0] - textsize[0])
        else : xcenter = center[0]
        if center[1] - textsize[1] < 0 : 
            ycenter = max(0,center[1] + textsize[1])
        else : ycenter = center[1]
        center = (xcenter,ycenter)
        
        self.draw.text(center,text,font=font,fill=color)
    
    def DrawCross(self,icenter,length=10,color='red',linewidth=2):
        """Draws a Cross on an image"""
        from numpy import around
        import numpy as num
        from copy import copy
        import sys
        
        try:
            ncross = len(icenter[0])
        except TypeError : ncross = 1
        
        if ncross == 1 : centers = (copy(icenter),)
        else:
            # Assumption: when there are multiple centers, the argument
            # icenter is the product of a num.where over a 2D array.
            if isinstance(icenter,tuple) and \
            isinstance(icenter[0],num.ndarray) :
                pass
            else : sys.exit('icenter argument is not as expected.\nHalting at graphs.FindPeaks_graph')
            
            centers = []
            for i in range(ncross):
                centers.append((icenter[1][i],icenter[0][i]))
        
        for center in centers:
            center = (center[0],self.size[1]-center[1]) # origin in upper left corner
            xleft = center[0] - int(around(length/2.))
            xright = center[0] + int(around(length/2.))
            yleft = center[1] ; yright = yleft
            xdown = center[0] ; xup = xdown
            ydown = center[1] + int(around(length/2.))
            yup = center[1] - int(around(length/2.))
            hor = [(xleft,yleft),(xright,yright)]
            ver = [(xdown,ydown),(xup,yup)]
            
            self.draw.line(hor,fill=color,width=linewidth)
            self.draw.line(ver,fill=color,width=linewidth)

