import everest
import matplotlib.pyplot as plt
import os
import sys
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from astropy.timeseries import TimeSeries
from astropy.utils.data import get_pkg_data_filename
import numpy as np
import pyke
import logging
#%matplotlib inline

logger = logging.getLogger()
logger.setLevel(logging.ERROR)                                         #Prevents spam


if os.path.exists("EverestFits-KepFiltered.fits"):
  os.remove("EverestFits-KepFiltered.fits")
  print("EverestFits-KepFiltered.fits clobbered")
else:
  print("The Filtered does not exist")

star = everest.Everest(201128338, season = 102, clobber = True)          #downloads the fitsfile to variable star

fits_file = get_pkg_data_filename(star.fitsfile)                       #star.fitsfile is the path to the star

print("Extension 0:")                                                  # header
print(repr(fits.getheader(fits_file, 0)))
print()
print("Extension 1:")                                                  # data
print(repr(fits.getheader(fits_file, 1)))

#fits.setval(fits_file, 'TTYPE1', value='TIME', ext=1)
fits.setval(fits_file, 'TTYPE2', value='SAP_FLUX', ext=1)              #likely that TTYPE6 should be SAP_FLUX

#possible code edit: prompt the user to change

pyke.kepfilter(star.fitsfile,                                          #filters the fits file and saves it as EverestFits-KepFiltered.fits
               passband='high',
               outfile="EverestFits-KepFiltered.fits",
               overwrite=True,
               datacol='SAP_FLUX',
               function='boxcar',
               logfile='kepfilter.log')


ts=TimeSeries.read(star.fitsfile,format='kepler.fits')                 #creates 2 AstroPY timeseries that can be used
ts1=TimeSeries.read("EverestFits-KepFiltered.fits",format='kepler.fits')

image = fits.open(star.fitsfile)
plt.imshow(image[5].data)                                              #shows the star

print("Pre KepFilter")
plt.plot(ts.time.jd, ts['sap_flux'], 'k.', markersize=1)
plt.xlabel('Julian Date')
plt.ylabel('Raw Flux (e-/s)')
plt.show()

print("Post KepFilter")
plt.plot(ts1.time.jd, ts1['sap_flux'], 'k.', markersize=1)
plt.xlabel('Julian Date')
plt.ylabel('Raw Flux (e-/s)')
plt.show()

model = BoxLeastSquares(ts1.time, ts1["sap_flux"])                     #creates a box least squares periodogram on the filtered timeseries
john1 = model.autopower(.2)

plt.plot(john1.period, john1.power)                                    #plots the periods vs the power
plt.show

n = john1.period[np.argmax(john1.power)]                               #finds the maximum power in the timeseries->n

print("The period is " + str(n) + " days")                             #prints out the period



