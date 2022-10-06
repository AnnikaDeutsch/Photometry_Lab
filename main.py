from astroquery.astrometry_net import AstrometryNet
from astroquery import jplhorizons
from astropy.io import fits
from astropy import units as u
import photutils
from glob import glob
import astropy
import astropy.config as config
import numpy
import matplotlib.pyplot as plt

#(image-bias) / [(flat-bias: individual images) / (flat: average of all images)]

biasimages = glob('./2022-09-01/Flat_V-*bias.fit')
flatimages = glob('./2022-09-01/Flat_V-*flat.fit')
biaslists = []
flatlists = []
for i in biasimages:
    biaslists.append(fits.open(i))
for i in flatimages:
    flatlists.append(fits.open(i))

biasaverage = numpy.zeros((512, 2048))
for k in biaslists:
    biasaverage += k[0].data
biasaverage = biasaverage / len(biaslists)

flataverage = numpy.zeros((512, 2048))
for k in flatlists:
    flataverage += k[0].data
flataverage = flataverage / len(flatlists)

flatcorrect = numpy.zeros((512, 2048))
for j in flatlists:
    flatcorrect += (j[0].data - biasaverage)
flatcorrect = flatcorrect / flataverage

jdates = []
images = glob('./New FITS Headers/WCS CCD Image *.fits')
imagelists = []
for i in images:
    imagelists.append(fits.open(i))
for j in imagelists:
    jdates.append(j[0].header["JD"])
print(jdates)
data = []
for k in jdates:
    bourgeois = jplhorizons.Horizons(id='Bourgeois', location='H81', epochs=k, id_type='smallbody')
    loc = bourgeois.ephemerides(quantities=1)
    datatemp = loc.iterrows('RA', 'DEC')
    for i in datatemp:
        data.append(i)
print(data)
ascensions = []
for i in data:
    ascensions.append(i[0])
decs = []
for i in data:
    decs.append(i[1])
positions = astropy.coordinates.SkyCoord(ascensions*u.deg, decs*u.deg)
apertur = photutils.aperture.SkyCircularAperture(positions, 5.*u.arcsec)
for t in imagelists:
    phot_table = photutils.aperture_photometry(t[0].data / flatcorrect, apertur, wcs=astropy.wcs.WCS(imagelists[0][0].header))
phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
print(phot_table)
plt.plot(phot_table['aperture_sum'])
plt.show()


