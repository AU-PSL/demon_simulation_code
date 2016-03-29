#-------------read_data.py-----------------------------------------------------#
#
# Purpose: This file is an example that explains how to use the astropy package
#          to read data from the data.fits file
#
#   Notes: Requires the astropy package
#          Please refer to the following website for more info:
#              http://docs.astropy.org/en/stable/io/fits/#module-astropy.io.fits
#------------------------------------------------------------------------------#

from astropy.io import fits

# first, we need to get the hdulist
# with DEMON, there are three hdu's: PRIMARY, CLOUD, and TIME_STEP
hdulist = fits.open('data.fits')

print("The structure of your fits file is:")
print(hdulist.info())

# Now we need to grab the data from the TIME_STEP hdu
data = hdulist[2].data

# Now we have all the data for the 100 particles worth of data in the following
# format:
#    data[dt][0]     -- time at dt timestep
#    data[dt][1][n]  -- x position of particle n at dt timestep
#    data[dt][2][n]  -- y position of particle n at dt timestep
#    data[dt][3][n]  -- x velocity of particle n at dt timestep
#    data[dt][4][n]  -- y velocity of particle n at dt timestep

# prints the x position of the 2nd timestep for the 4th particle
print("The x position of the 2nd timestep for the 4th particle is:")
print(data[1][2][3])
