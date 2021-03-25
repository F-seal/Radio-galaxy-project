import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.coordinates import SkyCoord
import sys

NSIDE = int(sys.argv[1])

NPIX = hp.nside2npix(NSIDE)
filename = './mask/mask_' + str(NSIDE) + '.fits'

mask = np.ones(NPIX,dtype=int)

f = open('./data/blake_mask22.cat')

bs = np.ones([22,3])

for i in range(22):
    temp = f.readline()
    bs[i,0] = float(temp[0:6])
    bs[i,1] = float(temp[8:13])
    bs[i,2] = float(temp[17:20])

# for ra,dec in bs:
#     idx = hp.ang2pix(NSIDE,ra,dec,lonlat=True)
#     mask[idx] = 0

for i in range(NPIX):
    RA,DEC = hp.pix2ang(NSIDE,i,lonlat=True)
    if DEC <= -40:
        mask[i] = 0
    else:
        eq = SkyCoord(RA,DEC, frame='icrs', unit='deg')
        l, b = eq.galactic.l.value, eq.galactic.b.value
        if np.abs(b) < 10:
            mask[i] = 0
    for ra,dec,r in bs:
        print(r)
        if (RA-ra)**2+(DEC-dec)**2 < r**2:
            mask[i] = 0


hp.write_map(filename,mask,coord='E',overwrite=True,dtype=int)
hp.mollview(mask)
plt.show()