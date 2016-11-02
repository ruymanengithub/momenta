#! /usr/bin/env python

"""Script to measure masses with spatial resolution on galaxies, using
images in differents bands. Aimed for GOODS data only."""

# IMPORT STUFF
import numpy as num
from pdb import set_trace as stop
# END IMPORT

# DATA

# FJy = 3631. * 10.**(mag / (-2.5))

# GOODS-ACS

goods_bands = num.array(['B435','V606','i775','z850'])
goods_cw = num.array([4317.40 ,5917.68,7693.03,9054.77])
goods_bw = num.array([293.47,672.31,434.60,593.43])

zeroGOODS = {'B435':25.65288,'V606':26.49341,'i775':25.64053,
'z850':24.84315}

sloan_bands = num.array(['u','g','r','i','z'])
sloan_cw = num.array([3551.,4686.,6165.,7481.,8931.])
ABminusSDSS = num.array([0.04,0.,0.,0.,0.02]) #
sloan_sun = num.array([6.80,5.15,4.67,4.56,4.53])
sloan_sun = sloan_sun + ABminusSDSS

AllMLcoeffs = {
'ug':{'g':(-0.221,0.485),'r':(-0.099,0.345),'i':(-0.053,0.268),
'z':(-0.105,0.226)},
'ur':{'g':(-0.390,0.417),'r':(-0.223,0.299),'i':(-0.151,0.233),
'z':(-0.178,0.192)},
'ui':{'g':(-0.375,0.359),'r':(-0.212,0.257),'i':(-0.144,0.201),
'z':(-0.171,0.165)},
'uz':{'g':(-0.400,0.332),'r':(-0.232,0.239),'i':(-0.161,0.187),
'z':(-0.179,0.151)},
'gr':{'g':(-0.499,1.519),'r':(-0.306,1.097),'i':(-0.222,0.864),
'z':(-0.223,0.689)},
'gi':{'g':(-0.379,0.914),'r':(-0.220,0.661),'i':(-0.152,0.518),
'z':(-0.175,0.421)},
'gz':{'g':(-0.367,0.698),'r':(-0.215,0.508),'i':(-0.153,0.402),
'z':(-0.171,0.322)},
'ri':{'g':(-0.106,1.982),'r':(-0.022,1.431),'i':(0.006,1.114),
'z':(-0.052,0.923)},
'rz':{'g':(-0.124,1.067),'r':(-0.041,0.780),'i':(-0.018,0.623),
'z':(-0.041,0.463)}}

def MassMap(goods_img,redshift):
    """Returns a mass map (on Msun units), given redshift and stamps 
    from GOODS-ACS on all four bands. If not possible, returns 'None'.
    
    Bibliography : 
    Bell, E. F., McIntosh, D. H., Katz N. & Weinberg M. 2003.
    'The K correction', Hogg, Baldry, Blanton & Eisenstein 2002.
   'Cocina fresca', Karlos Arguiñano.
    """
    # IMPORT STUFF
    from Zscope.zlib1 import cosmocalc
    from Mass2D import interpol_img, get_colmap, get_Lmap
    # END IMPORT
   
    drop,drop,DMOD = cosmocalc(redshift)

    goods_img_muJy = goods_img * 0.
    
    # Converting images to micro-Jy (muJy)
    
    for i in range(len(goods_bands)):
        zero = zeroGOODS[goods_bands[i]]
        imgmuJy = 1.e6 * 3631. * goods_img[i] * 10.**(zero/(-2.5))
        goods_img_muJy[i] = imgmuJy.copy()
    del imgmuJy
    
    # Rest-frame central wavelengths
    
    goods_cw_r = goods_cw / (1.+redshift)
    
    # Which of the observed GOODS-bands bracket rest-frame sloan bands?
    
    inter_cw_ix = num.where((goods_cw_r <= sloan_cw.max()) & \
    (goods_cw_r>=sloan_cw.min()))

    if len(inter_cw_ix[0]) != 0:
        inter_cw = sloan_cw[inter_cw_ix]
        inter_bands = sloan_bands[inter_cw_ix]
    else:
        pass

    if len(inter_bands) < 2 : return None    
    
    # Which sloan colors could be built given GOODS-data and redshift?
    fetch_colors = []
    for i in range(len(inter_bands)-1,-1,-1):
        for j in range(i-1,-1,-1):
            addcolor = inter_bands[j]+inter_bands[i]
            if addcolor not in fetch_colors : fetch_colors.append(addcolor)

    possiblecolors = [color for color in fetch_colors if color in MLcoeffs]
    if len(possiblecolors)>0: masscolor = possiblecolors[0]
    else: return None
    
    # Selection of 'reddest' color possible
    
    MLcoeffs = AllMLcoeffs[masscolor][masscolor[-1]]
    interimgs = []
    
    # Now interpolation of GOODS images to the wavelengths which match
    # sloan cw's at z=0.

    for i in range(2):
        ix = num.where(sloan_bands == masscolor[i])[0]
        intercw = sloan_cw[ix] * ( 1.+ z )
        cwbefix = num.where(intercw>=goods_cw)[0][-1]
        cwbef = goods_cw[cwbefix] ; cwaft = goods_cw[cwbefix+1]
        bandbef = goods_bands[cwbefix]
        bandaft = goods_bands[cwbefix+1]
        imgbef = goods_img_muJy[cwbefix]
        imgaft = goods_img_muJy[cwbefix+1]
        interimgs.append(interpol_img(intercw,imgbef,imgaft,cwbef,cwaft))
    
    if None in interimgs : return None
    
    # Luminosity map
    
    zeros = 0.,0. # as interpolated images are already in muJy
    colmap = get_colmap(interimgs,zeros)
    Magsun = sloan_sun(bandaft) # absolute magnitude
    Limg = interimgs[1]
    Lresun = get_Lmap(Limg,Magsun,Dmod,redshift)
    goodmask = ((Lresun > 0.) & (colmap>-99.))
    good = num.where(goodmask)
    
    # And Mass Map

    ax = MLcoeffs[0] ; bx = MLcoeffs[1]
    Mmap = Lresun * 0. - 99.
    Mmap[good] = Lresun[good] * 10.**(ax+bx * colmap[good])

    return Mmap,Lresun,Colmap,goodmask
    
    # Dev Notes:
    #How Lresun is computed:
    #Dmod = 5 * log10(Dl) - 5.
    #Dl = 10.**((Dmod +5)/5.)

    #ABimg = zero -2.5 * log10(img)
    #Jy = 3631. * 10.**(ABimg / (-2.5)) = 3631. * img * 10**(zero/(-2.5))
    #L = Jy * 4 * pi * Dl**2./(1.+z)
    #Lresun = L /Lsun
    #Jysun = 3631. * 10.**(Msun/(-2.5))
    #Lsun = Jysun 4 * pi * 10.**2.
    #Lresun = img * (10.**((zero-Msun)/-2.5)) * (Dl/10.)**2. / (1.+z)


def interpol_img(intercw,imgbef,imgaft,cwbef,cwaft):
    """Interpolates (strictly linear) between 2 images in wavelength."""
    if (cwbef > intercw) or (cwaft < intercw) : return None
    if (cwbef == intercw): return imgbef
    if (cwaft == intercw): return imgaft
    return intercw*(imgaft-imgbef)/(cwaft-cwbef) + imgbef
    
def get_colmap(imgs,zeros):
    """Returns a color map."""
    img1 = imgs[0] ; img2 = imgs[1]
    zero1 = zeros[0] ; zero2 = zeros[1]
    positive = num.where((img1>0.) & (img2>0.))
    colmap = imgs * 0. -99.0
    colmap[positive] = zero1-zero2 - 2.5 * num.log10(img1/img2)
    return colmap

def get_Lmap(img,Msun,Dmod,redshift):
    """Returns Luminosity map, pixel by pixel, in solar luminosities in
    such band. 'img' in muJy."""
    positive = num.where(img>0.)
    Lmap = img * 0. - 99.0
    Dl = 10.**((Dmod +5)/5.)
    zero = -2.5*log10(1.e6)
    Lresun[positive] = img[positive] * \
    (10.**((zero-Msun)/-2.5)) * (Dl/10.)**2. / (1.+z)
    return Lresun
