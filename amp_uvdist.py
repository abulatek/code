from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

def amp_uvdist(sbfile,lbfile,binsize):
    plt.clf()
    filenames = [sbfile,lbfile]

    for name in filenames:
        image = fits.open(name)
        u = image[0].data['UU']
        v = image[0].data['VV']
        u *= image[0].header['crval4']/1e3
        v *= image[0].header['crval4']/1e3
#        plt.plot(u,v,'.')

        vis = image[0].data['data']
        real = (vis[:,0,0,0,:,0,0] + vis[:,0,0,0,:,1,0])/2.
        real = np.mean(real,axis=1)
        imag = (vis[:,0,0,0,:,0,1] + vis[:,0,0,0,:,1,1])/2.
        imag = np.mean(imag,axis=1)

        amp = np.sqrt(real**2 + imag**2)
        amp = amp.squeeze()
        uvdist = np.sqrt(u**2 + v**2)
#        plt.plot(uvdist,amp,'.')

###### for binning!

        newuvdist = []
        newreal = []
        newimag = []
        newamp = []
        rstdevs = []
        istdevs = []
        rNs = []
        iNs = []
        amperrs = []
        rangemin = np.amin(uvdist)
        rangemax = np.amax(uvdist)

        for minimum in range(rangemin,rangemax,binsize):
            maximum = minimum + binsize
            uvsubarray = (uvdist > minimum) & (uvdist < maximum)
            uvdist2 = uvdist[uvsubarray]
            uvrep = np.median(uvdist2)
            newuvdist.append(uvrep)

            realsubarray = (real[uvsubarray] != 0)
            real2 = real[uvsubarray][realsubarray]
            realav = np.mean(real2)
            newreal.append(realav)

            imagsubarray = (imag[uvsubarray] != 0)
            imag2 = imag[uvsubarray][imagsubarray]
            imagav = np.mean(imag2)
            newimag.append(imagav)

            rstdev = np.std(real2)
            istdev = np.std(imag2)

            rN = len(real2)
            iN = len(imag2)

            amperr = np.sqrt((rstdev**2/rN)*(realav**2/(realav**2+imagav**2))+(istdev**2/iN)*(imagav**2/(realav**2+imagav**2)))
            amperrs.append(amperr)

        newuvdist = np.asarray(newuvdist)
        newuvdist = newuvdist[~np.isnan(newuvdist)]
        r = np.asarray(newreal)
        r = r[~np.isnan(r)]
        i = np.asarray(newimag)
        i = i[~np.isnan(i)]
        newamp = np.sqrt(r**2 + i**2)
        amperrs = np.asarray(amperrs)
        amperrs = amperrs[~np.isnan(amperrs)]

        plt.plot(newuvdist,newamp,'.')
        plt.errorbar(newuvdist,newamp,yerr=amperrs,fmt='.')
        plt.xlabel('Distance from center of uv-plane (klambda)')
        plt.ylabel('Amplitude (Jy)')
        plt.title('Amplitude versus uv-distance')
        plt.grid(True)
        plt.show()

