#!/usr/bin/env python

"""
    This example illustrates how to obtain an MGE fit from a galaxy image
    using the mge_fit_sectors package and how to verify the results.

    V1.0.0: Translated NGC4342 example from the corresponding IDL version.
        Michele Cappellari, Aspen Airport, 8 February 2014
    V1.0.1: Fixed incorrect offset in high-res contour plot.
        Use arcsec pixel scale. MC, Oxford, 18 September 2014
    V1.1.0: Translated M32 example of fitting to two images simultaneously.
        MC, Oxford, 18 June 2015
    V1.1.1: Support both Pyfits and Astropy to read FITS files.
        MC, Oxford, 23 October 2015
    V1.1.2: Make files paths relative to this file, to run the example
        from any directory. MC, Oxford, 20 January 2017
    V1.1.3: Included fit_1d() example. Added important note about
        mge_fit_sectors_regularized(). MC, Oxford, 14 February 2017
    V1.2.0: Illustrates how to efficiently mask an object before MGE fitting.
        Included dist_circle() function. MC, Oxford, 17 March 2017
    V1.3.0: Modified comment about import of mge_fit_sectors_regularized.
        Thanks to Evgeny Vasilyev (Oxford) for the suggestion.
        Included fit_ngc5831_twist() example and corresponding imports.
        MC, Oxford, 27 July 2017
    V1.3.1: Make path relative to package to run the example from any directory.
        MC, Oxford, 17 April 2018


"""

import numpy as np
import math as math
import matplotlib.pyplot as plt
from astropy.io import fits
from os import path

import mgefit
from mgefit.find_galaxy import find_galaxy
from mgefit.mge_fit_1d import mge_fit_1d
from mgefit.sectors_photometry import sectors_photometry
#from mgefit.mge_fit_sectors import mge_fit_sectors
from mgefit.mge_fit_sectors_regularized import mge_fit_sectors_regularized
from mgefit.mge_print_contours import mge_print_contours
from mgefit.mge_fit_sectors_twist import mge_fit_sectors_twist
from mgefit.sectors_photometry_twist import sectors_photometry_twist
from mgefit.mge_print_contours_twist import mge_print_contours_twist

####---------------------------------------------------------------------------

def fit_ngc2824():
    """
    MGE fit to the 2MASS image for NGC2824.  
    """
    #### Set File:

    file_dir = path.dirname(path.realpath(mgefit.__file__))  # path of mgefit

    file = '/home/danielle/mge_fitting/mgefit/images/ngc_2824_2mass_el.fits'

    #### SFind Galaxy outputs and Constants:

    eps    = 0.334
    theta  = 113.8
    xpeak  = 43
    ypeak  = 44

    skylev = 104.307          # Counts/pixel

    seesh  = 1  	      # Seeing shape parameter from the telescope.

    qBound = 0.48     	      # Lower q-bound

    MAGZP  = 20.8861          # This is the magnitude zero point listed in the fits header of the 
    Ak     = 0.01             # Galactic extinction parameter.  Find from NED

    D      = 40.7*10**6       # Distance to NGC 4476 in parsecs

    RCO    = 3600.00          # From the table in the overleaf document, this number needs to be in
                              # parsec.

    Mdyn   = 3.4*10**10       # From the table in the overleaf document


    #### Other Constants:

    ang       = 90-theta
    scale     = 1.0                 # Arcsec/pixel
    minlevel  = 2
    n         = 15

    d         = D*4.8481*10**(-6)   # Parsecs per pixel width so parsecs/arcsec, d=D*tan(theta),
                                    # where theta=1 arcsec, theta needs to be in radians.  

    Msk       = 3.28                # Units are magnitudes

    xlim      = d*40                # CO radius for plot

    ngauss = 15

    hdu = fits.open(file)
    img = hdu[0].data   # So you need to be careful here.  It is not always the hdu[0] that
                        #contains the pixel values.  In the case of this ukidss image it is hdu[1]

    img -= skylev       # subtract sky

##### Find Galaxy:
    plt.clf()
    f = find_galaxy(img, nblob=1, plot=1)
    plt.pause(1)  # Allow plot to appear on the screen
    plt.savefig("n2824_2mass.png")

    # Here we use a single gaussian MGE PSF 
    FWHM=3.13*seesh-0.46 
    # https://irsa.ipac.caltech.edu/data/2MASS/docs/supplementary/seeing/seesum.html
    sigmapsf=FWHM/2.355  # Must divide by the gaussian factor of 2.355 to get the sigma parameter
                         # that is wanted in the MGE models 

    # Perform galaxy photometry
    plt.clf()
    s = sectors_photometry(img, f.eps, f.theta, f.xpeak, f.ypeak,
                           minlevel=minlevel, plot=1)
    plt.pause(1)  # Allow plot to appear on the screen
    plt.savefig("n2824_pgp_2mass.png")

##### Do the actual MGE fit
    qbounds=[qBound,1]
    plt.clf()
    m = mge_fit_sectors_regularized(s.radius, s.angle, s.counts, f.eps,
                        ngauss=ngauss, sigmapsf=sigmapsf,
                        scale=scale, plot=1, qbounds=qbounds)

    plt.pause(1)                      # Allow plot to appear on the screen
    plt.savefig("n2824_mfs_2mass.png") # Show contour plots of the results
    plt.clf()
    plt.subplot(121)
   
    # You may have to play around with the binning and magrange parameter for each galaxy.  
    mge_print_contours(img.clip(minlevel), f.theta, f.xpeak, f.ypeak, m.sol, scale=scale,
                       binning=4, sigmapsf=sigmapsf, magrange=5.5)

    # Extract the central part of the image to plot at high resolution.
    # The MGE is centered to fractional pixel accuracy to ease visual comparson.

    img = img[f.xpeak-n:f.xpeak+n, f.ypeak-n:f.ypeak+n]
    xc, yc = n - f.xpeak + f.xmed, n - f.ypeak + f.ymed
    plt.subplot(122)
    mge_print_contours(img, f.theta, xc, yc, m.sol,
                       sigmapsf=sigmapsf, scale=scale)
    plt.pause(1)  # Allow plot to appear on the screen
    plt.savefig("n2824_cp_2mass.png")
    plt.clf()
    print('  Counts for each Gaussian: ', m.sol[0])
    print('  Sigma_Pixels for each Gaussian: ', m.sol[1])
    print('  qobs for each Gaussian: ', m.sol[2])


    #### MGE Outputs:

    c1 = m.sol[0]
    c2 = m.sol[1]
    c3 = m.sol[2]


####---------------------------------------------------------------------------------------------
# MGE Analysis: 
    Cf = []
    uk = []
    for c in range(0, len(c1)):
        ###### We begin by converting the Total Counts to a Peak Surface brightness.  There is
             # probably a much more eloquent way to code this.
        Cfc = c1[c]/(2*3.14159*c3[c]*(c2[c])**2)
        Cf.append(Cfc)
        ###### Now we use the peak Surface brightness to a Johnson-Cousins K-band surface      
             # brightness in magnitudes/arcsec^2
        ukc = MAGZP-2.5*math.log10(Cf[c])-Ak  #MAGZP-2.5*math.log10(C1)-Ak  Whats C1?
        uk.append(ukc)


    Ip = []
    sigp = [] # sig in arcsec
    sigdp = [] #sig in parsec
    for c in range(0, len(c1)):
        ###### Now we use the K-band surface brightness in magnitudes/arcsec^2 to calculate the a
             # surface density in Solar Luminosities per parsec^2. 
        Ipc = (64800/3.14159)**2*(10**(0.4*(Msk-uk[c])))
        Ip.append(Ipc)
        ###### Now we use the K-band surface density in Solar Luminosities per parsec^2 to
             # calculate the total luminosity of each Gaussian. 
        ###### Step 1.  In the axes semetric case sigmaj=sigmaj', but we want to convert the
             # SigmaPixels from pixels to arcseconds by multipling be the scale
        ###### Which has units of arsec/pixel.  For the final Gaussian inputs we need sigp to nbe
             # in units of parsecs
        sigpc = c2[c]*scale
        sigp.append(sigpc)
        sigd = c2[c]*scale*d
        sigdp.append(sigd)

    Ltot = []
    for c in range(0, len(c1)):
        ###### Step 2.  Calculate Lj for each Gaussian, we need Lj to be in units of Lsol so we
             # need the distance.  
        Ltotc = 2*3.1415*Ip[c]*sigp[c]**2*c3[c]*d**2
        Ltot.append(Ltotc)


    int = []
    for c in range(0 ,len(c1)):
        ###### Step 3.  We now need to evalutate the total light interior to the RCO.  We can
             # solve the tripple integral analytically to yield:
        intc = Ltot[c]*(1-np.exp(-1*RCO**2/(2*sigdp[c]**2)))
        int.append(intc)
        
    inttot = 0
    for c in range(0, len(c1)):
        inttot += int[c]

    ###### Step 4.  Calculate the M/LK,MGE ratio, this is the Mdyn/L ratio
    M2L = Mdyn/inttot

    ###### Step 5.  Now the luminosity density of each Gaussian is given by:
    z1 = 0.0     #Zero parsecs above the midplane. z itself is in arcsec so we really need to use
                 #the d from above.  
    r = np.arange(0, 1000, 0.1) # So this is r measures in arcsec
    z2 = 150/d   #150 pc above the midplane at the distance to the galaxy  
    r1 = np.arange(0, 100, 0.1) # So this is r measures in arcsec
    z3 = 300/d   #300 pc above the midplane at the distance to the galaxy.   
    r2 = np.arange(0, 100, 0.1) # So this is r measures in arcsec


    LD1 = []
    LD2 = []
    LD3 = []
    for c in range(0 ,len(c1)):
        LD1c = (Ltot[c]/(sigp[c]*2.506628275*c3[c]*d)**3) * np.exp(-1/(2*(sigp[c]*sigp[c]))*(r**2+(z1/c3[c])**2))
        LD1.append(LD1c)
        LD11c = (Ltot[c]/(sigp[c]*2.506628275*c3[c]*d)**3) * np.exp(-1/(2*(sigp[c]*sigp[c]))*(r1**2+(z2/c3[c])**2))
        LD2.append(LD11c)
        LD12c = (Ltot[c]/(sigp[c]*2.506628275*c3[c]*d)**3) * np.exp(-1/(2*(sigp[c]*sigp[c]))*(r2**2+(z3/c3[c])**2))
        LD3.append(LD12c)

    pd1 = 0     #Expression for the total stellar volume density Msol/pc^3
    pd2 = 0     #Expression for the total stellar volume density Msol/pc^3
    pd3 = 0     #Expression for the total stellar volume density Msol/pc^3
    for c in range(0, len(LD1)):
        pd1 += LD1[c]*M2L
        pd2 += LD2[c]*M2L
        pd3 += LD3[c]*M2L

##### All printing to terminal:
    print('  Co for each Gaussian: ', Cf)
    print('  uk for each Gaussian: ', uk)  
    print('  Distance to NGC 2824 in parsecs: ', D)
    print('  parsec/pixel for NGC 2824: ', d)
    print('  Ip for each Gaussian: ', Ip)
    print('  sig for each Gaussian in arcsec: ', sigp)
    print('  sig for each Gaussian in parsecs: ', sigdp)
    print('  The Peak Luminosity, Lj, for each Gaussian: ', Ltot)  
    print('  Ltot @RCO and z=inf of each Gaussian: ', int)
    print('  Ltot @RCO and z=inf: ', inttot)   
    print('  Dynamical K-band mass to light ratio: ', M2L)

######  Make a plot of the pD versus radius for z=0, z=3.5.  
    plt.plot(r*d, np.log10(pd1), 'r--',r1*d, np.log10(pd2), 'b-.',r2*d, np.log10(pd3), 'g:')
    plt.xlabel('Radius (parsecs)')
    plt.ylabel(r'(Log$_{10}$($\rho$) ($L_{\odot}pc^{-2}$)')
    plt.title('NGC 2824: 2MASS')
    plt.grid(False)
    plt.xlim(0,xlim) 
    plt.ylim(-3,2)  # These will have to be scaled for each galaxy. Start by running the code
                    # unset and then refine.  
    plt.savefig("LD_2mass_N2824.png")
    plt.pause(1)

####----------------------------------------------------------------------------

if __name__ == '__main__':

    print("\nFitting mge NGC2824-----------------------------------------\n")
    fit_ngc2824()
