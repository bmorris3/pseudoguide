# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 20:47:52 2014

@author: bmorris3
"""

import numpy as np
from scipy.signal import correlate2d

def gauss2d(x_centroid, y_centroid, amplitude, width, dimensions):
    '''
    Generate a 2D gaussian.
    
    Parameters
    ----------
    x_centroid : float
        centroid in x coordinates
    y_centroid : float
        centroid in y coordinates
    amplitude : float
        peak value of gaussian at centroid
    width : float
        gaussian sigma parameter
    dimensions : int
        number of pixels along each side of square (to generate)
    '''
    x, y = np.meshgrid(range(dimensions),range(dimensions))
    return amplitude*np.exp(-1*( (x-x_centroid)**2/(2*width**2) + (y-y_centroid)**2/(2*width**2) ))

def createimage(x_centroids, y_centroids, dimensions):
    '''
    Create a simulated image with stars at centroids given by `x_centroids` and `y_centroids`.
    
    Parameters
    ----------
    x_centroids : array-like of dtype=float
        centroids in x coordinates
    y_centroids: array-like of dtype=float
        centroids in y coordinates
    '''
    backgroundlevel = 100
    backgroundimage = np.sqrt(backgroundlevel)*np.random.randn(dimensions,dimensions) + backgroundlevel
    PSFwidth = 1.5
    amplitude = 1*backgroundlevel
    simulateddata = backgroundimage
    for i in range(len(x_centroids)):
        simulateddata += gauss2d(x_centroids[i], y_centroids[i], amplitude, PSFwidth, dimensions)
    return simulateddata

def gensourcemask(image, PSFwidth):
    '''
    Generate a source mask from the image `image`, assuming PSF width `PSFwidth`.
    
    Parameters
    ----------
    image : array-like of dtype=int
        centroids in x coordinates
    PSFwidth: float
        estimated with of the PSF in pixels

    Returns
    -------
    sourcemask: array-like of dtype=int
        Matrix where each element with a source is equal to unity, each element
        without a source is equal to zero.

    '''
    stampdim = 10*int(PSFwidth)
    crosscorrelation = correlate2d(image, gauss2d(stampdim/2, stampdim/2, 10, PSFwidth, stampdim), mode='same')
    crosscorrelation /= np.median(crosscorrelation) # Normalize in place
    sourcemask = np.zeros_like(crosscorrelation) # Intialize output matrix
    sourcemask[crosscorrelation > 1.2] = 1
    return sourcemask

def safetylimits(shift_0, shift_1, binning, platescale):
    '''
    Given a candidate shift suggestion, prevent the telescope from moving too far
    if the candidate shift is large.
    
    Parameters
    ----------
    shift_0 : int
        number of pixels to move in the zero axis
    shift_1 : int
        number of pixels to move in the one axis
        
    Returns
    -------
    safeshift_0: int
        shift along zero axis after limits implemented
    safeshift_1: int
        shift along first axis after limits implemented
    '''
    if -60 <= shift_0*binning*platescale <= 60 and -60 <= shift_1*binning*platescale <= 60: # limit to within 1 arcminute shifts
        return shift_0*binning*platescale, shift_1*binning*platescale
    else:
        return 0, 0
        

def findshifts(image1, image2, binning, platescale):
    '''
    Find the best offset in units of arcseconds along each axis to shift 
    `image2` into the position of `image1`. Also check that the shift
    is small using safetylimits(), if too large an offset, force shift to zero.
    
    Parameters
    ----------
    image1 : array-like of dtype=int
        Source mask of the earlier image
    image2: array-like of dtype=int
        Source mask of the later image
        
    Returns
    -------
    safeshift_0 : int
        Best number of arcseconds to shift along axis 0 of `image2` to overlay it on `image1`
    safeshift_1 : int
        Best number of arcseconds to shift along axis 1 of `image2` to overlay it on `image1`
    '''
    dimensions = np.shape(image1)[0]
    mask1_0 = np.concatenate([np.zeros(dimensions/2.),np.sum(image1,axis=0),np.zeros(dimensions/2.)])
    mask1_1 = np.concatenate([np.zeros(dimensions/2.),np.sum(image1,axis=1),np.zeros(dimensions/2.)])
    mask2_0 = np.concatenate([np.zeros(dimensions/2.),np.sum(image2,axis=0),np.zeros(dimensions/2.)])
    mask2_1 = np.concatenate([np.zeros(dimensions/2.),np.sum(image2,axis=1),np.zeros(dimensions/2.)])
    rollrange = range(-1*dimensions/2, dimensions/2)

    chi2_0 = np.zeros_like(rollrange)
    for i, rollint in enumerate(rollrange):
        chi2_0[i] = np.sum((mask1_0 - np.roll(mask2_0, rollint))**2)    
    bestroll_0 = rollrange[np.argmin(chi2_0)]

    chi2_1 = np.zeros_like(rollrange)
    for i, rollint in enumerate(rollrange):
        chi2_1[i] = np.sum((mask1_1 - np.roll(mask2_1, rollint))**2)    
    bestroll_1 = rollrange[np.argmin(chi2_1)]
    
    # Check to make sure they're within the limits   
    safeshift_0, safeshift_1 = safetylimits(bestroll_0, bestroll_1, binning, platescale)
    
    return safeshift_0, safeshift_1