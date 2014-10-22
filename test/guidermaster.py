'''
Execute the script datagenerator.py first, then run guidermaster.py. 
'''
from glob import glob
import pyfits 
import sys
sys.path.append('../')
import os
import pseudoguide as pg
from time import sleep

# Potentially tweakable parameters
estimatedPSF = 1.5      # very rough estimate of PSF width in pixel units
platescale = 0.656      # arcseconds per pixel
guideron = True         # Leaving this switch in for later
rootdir = 'data'        # Directory where to search for the files to guide on
searchstring = "0*.fi*" # If you were to `ls` this string, you should get the files
waitbetweenchecks = 5   # seconds between checks for new files (minimizes CPU usage)
logfilepath = 'guideroffsetlog.txt' # where to save the guider offset log file

# Initializing some variables
firstimage = None
logfile = open(logfilepath,'w')
oldfilelist = []

while guideron:
    sleep(waitbetweenchecks)
    print 'Checking for new files (reminder: this is an infinite loop)'
    # Make list of files currently in output directory
    search_dir = os.path.join(rootdir,searchstring)
    currentfilelist = filter(os.path.isfile, glob(search_dir))
    currentfilelist.sort(key=lambda x: os.path.getmtime(x))
    
    # If a new file has been found:
    if currentfilelist != oldfilelist:
        # If it's the first file, make a source mask from it
        if firstimage is None:
            print 'Caught first file'
            firstimage = pyfits.getdata(currentfilelist[0])
            binning = pyfits.getheader(currentfilelist[0])['XBINNING']
            init_sourcemask = pg.gensourcemask(firstimage, estimatedPSF)
        # If its not first, compare it to the first and calculate a shift
        else: 
            newimagepath = list(set(currentfilelist) - set(oldfilelist))[0]
            print 'new file detected = %s' % newimagepath
            newimage = pyfits.getdata(newimagepath)
            latest_sourcemask = pg.gensourcemask(newimage, estimatedPSF)
            xcorr_sky, ycorr_sky = pg.findshifts(init_sourcemask, latest_sourcemask, binning, platescale)
            print 'X correction = %.4f arcseconds, Y correction = %.4f arcseconds\n' % (xcorr_sky, ycorr_sky)
            
            # Write out the update to the log file
            with open('guideroffsetlog.txt', 'a') as logfile:
                logfile.write('%.4f\t%.4f\n' % (xcorr_sky, ycorr_sky))
        # Update the old file list with the new file
        oldfilelist = currentfilelist
