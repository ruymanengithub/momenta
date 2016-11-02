#! /usr/bin/env python

from time import time
import os
import string
from pdb import set_trace as stop

def save_img_ds9(myds9,imgname,format='jpeg'):
    """Creates a .png/.jpg with the image currently displayed on ds9."""
    # IMPORT STUFF
    
    # END IMPORT
    
    ready = raw_input('\nKick any key when ds9 is firmly on top\n')
    myds9.xpaset('saveimage %s %s' % (format,imgname))
    
    return None

def show_on_ds9(imgtuple,namestuple):
    """Shows several images on ds9. First make .fits."""
    # IMPORT STUFF
    from CommonTools.ds9 import ds9class
    import pyfits
    # END IMPORT
    
    if len(imgtuple) != len(namestuple) : return False
    myds9 = ds9class()
    
    isopen = myds9.isOpen()
    if isopen : myds9.xpaset('exit')
    myds9.launch()

    for i in range(len(imgtuple)):
        imgname = namestuple[i]
        pyfits.writeto(imgname,imgtuple[i])
    
    for i in range(len(namestuple)):
        myds9.xpaset('frame %i' % (i+1))
        myds9.xpaset('file %s' % namestuple[i])
        myds9.xpaset('scale mode 99.5')
        myds9.xpaset('zoom to 1.75')
        myds9.xpaset('cmap rainbow')
    myds9.xpaset('frame 1')
    
    return myds9


def get_coo_peak():
    """Function to get coordinates on an image displayed on ds9, selecting
    brightest pixel inside a 'box'."""
    # IMPORT STUFF
    
    # END IMPORT
    
def get_coo_imexam():
    """Function to get coordinates on an image displayed on ds9 by means
    of iraf.imexam."""
    # IMPORT STUFF
    from pyraf import iraf
    import pyraf
    import string
    # END IMPORT
    isthere = os.path.exists
    
    tmpcursor = 'tmpCursor_%f.txt' % time()
    
    def parse(infile):
        f = open(infile,'r') ; lines = f.readlines() ; f.close()
        pointers = [i for i in range(len(lines)) if 'COORDINATES' in lines[i]]
        try:
            pointer = pointers[0]
            coo = [float(item) for item in string.split(lines[pointer+2])[0:2]]
            return tuple(coo)
        except IndexError: return (-99.0,-99.0)
    
    iraf.rimexam.setParam('fittype','gaussian')
    iraf.rimexam.setParam('radius','5')
    
    satisf = False
    while not satisf:
        print 'Get over position with cursor and press "a"'
        print 'If no Error shows on screen, then click "q" to go on'
        try: 
            iraf.imexam(Stdout=tmpcursor)
            coo = parse(tmpcursor)
            os.system('rm %s' % tmpcursor)
        except pyraf.irafglobals.IrafError : 
            print '\nIraf Error raised!\n'
            coo = (-99.0,-99.0)
        print 'x=%.3f   y=%.3f' % coo
        satisf = raw_input('are you satisfied with resulting coordinates? yes=y ')
        satisf = string.lower(satisf) == 'y'
        
    if isthere(tmpcursor) : os.system('rm %s' % tmpcursor)
    
    iraf.rimexam.unlearn()
    iraf.imexam.unlearn()
    
    return coo
