'''
Execute the script datagenerator.py first, then run guidermaster.py. 
'''
import pyfits 
import numpy as np
from time import sleep
import sys
sys.path.append('../')
import pseudoguide as pg

dimensions = 500
binning = 2
platescale = 0.656 # arcseconds per pixel
xpositions = np.array([100, 300, 400, 250])
ypositions = np.array([100, 300, 250, 400])
driftperexp_0, driftperexp_1 = 1, 2.5

Nimages = 20
testimages = [pg.createimage(xpositions+i*driftperexp_0, ypositions+i*driftperexp_1, dimensions) for i in range(Nimages)]

## Clean up data directory first
import os 
folder = 'data/'
for the_file in os.listdir(folder):
    file_path = os.path.join(folder, the_file)
    try:
        if os.path.isfile(file_path) and '__init__' not in file_path:
            os.unlink(file_path)
    except Exception, e:
        print e

for i, testimage in enumerate(testimages):

    print 'Writing file data/%03d.fits' % i
    hdr = pyfits.Header()
    hdr['XBINNING'] = 2
    pyfits.writeto('data/%03d.fits' % i, testimage, header=hdr)
    print 'Pausing 20 seconds...'
    sleep(20)

