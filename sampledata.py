# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 11:36:49 2014

@author: bmorris3
"""

from glob import glob
import pyfits 
import numpy as np

platescale = 0.656 # arcseconds per pixel

files = sorted(glob('/local/tmp/ARCSAT/20141021/GJ*sdss_g*.fits'))
# Took sample images up until this one, then began run
firstgoodfile = '/local/tmp/ARCSAT/20141021/GJ1042_sdss_g_20141022_060257.fits'
goodfiles = files[files.index(firstgoodfile):]

firstimage = pyfits.getdata(goodfiles[0])
binning = pyfits.getheader(goodfiles[0])['XBINNING']
goodimages = np.zeros((len(goodfiles), firstimage.shape[0], firstimage.shape[1]))
for i, goodfile in enumerate(goodfiles):
    goodimages[i,:,:] = pyfits.getdata(goodfile)

import pseudoguide as pg
estimatedPSF = 1.5
init_sourcemask = pg.gensourcemask(goodimages[0,:,:], estimatedPSF)

for j, goodimage in enumerate(goodimages[1:,:,:]):
    latest_sourcemask = pg.gensourcemask(goodimage, estimatedPSF)
    xcorr_pixel, ycorrection_pixel = pg.findpixelshifts(init_sourcemask, latest_sourcemask)
    xcorr_sky, ycorrection_sky = np.array([xcorr_pixel, ycorrection_pixel])*2*platescale
    print 'X correction = %.2f arcseconds, Y correction = %.2f arcseconds\n' % (xcorr_sky, ycorrection_sky)

from matplotlib import pyplot as plt
fig, ax = plt.subplots(1,2)
ax[0].imshow(np.log(goodimages[0,:,:]), interpolation='nearest', origin='lower')
ax[1].imshow(np.log(np.sum(goodimages,axis=0)), interpolation='nearest', origin='lower')
plt.show()


