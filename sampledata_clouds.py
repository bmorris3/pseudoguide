# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 11:36:49 2014

@author: bmorris3
"""

from glob import glob
import pyfits 
import numpy as np

platescale = 0.656 # arcseconds per pixel

path = '/astro/users/bmmorris/ARCSATsampledata/'
files = sorted(glob(path+'GJ*sdss_g*.fits'))
# Took sample images up until this one, then began run
firstgoodfile = path+'GJ1042_sdss_g_20141022_060257.fits'
goodfiles = files[files.index(firstgoodfile):]
flat = pyfits.getdata(path+'/domeflat_sdss_g_004.fits')

firstimage = pyfits.getdata(goodfiles[0])
binning = pyfits.getheader(goodfiles[0])['XBINNING']
goodimages = np.zeros((len(goodfiles), firstimage.shape[0], firstimage.shape[1]))
for i, goodfile in enumerate(goodfiles):
    a = 1 - float(i+1)/len(goodfiles)
    b = float(i+1)/len(goodfiles)
    goodimages[i,:,:] = a*pyfits.getdata(goodfile) + b*flat

import paopg as p 
estimatedPSF = 1.5
init_sourcemask = p.gensourcemask(goodimages[0,:,:], estimatedPSF)

for j, goodimage in enumerate(goodimages[1:,:,:]):
    latest_sourcemask = p.gensourcemask(goodimage, estimatedPSF)
    print 'len(mask) = %d' % (np.sum(latest_sourcemask))
    xcorr_sky, ycorr_sky = p.findpixelshifts(init_sourcemask, latest_sourcemask, binning, platescale)
    #xcorr_sky, ycorrection_sky = np.array([xcorr_pixel, ycorrection_pixel])
    print 'X correction = %.4f arcseconds, Y correction = %.4f arcseconds\n' % (xcorr_sky, ycorr_sky)

from matplotlib import pyplot as plt
fig, ax = plt.subplots(1,2)
ax[0].imshow(np.log(goodimages[0,:,:]), interpolation='nearest', origin='lower')
#ax[1].imshow(np.log(np.sum(goodimages,axis=0)), interpolation='nearest', origin='lower')
ax[1].imshow(np.log(goodimages[len(goodfiles)-1,:,:]), interpolation='nearest', origin='lower')
plt.show()


