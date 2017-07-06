from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
###### In order to determine values to use for rangemin and rangemax,
###### use amp_uvdist.py to plot amp vs. uvdist and bound overlap.

def quotient(sbfile,lbfile,binsize,rangemin,rangemax):
    plt.clf()

###### short baselines

    image = fits.open(sbfile)
    u = image[0].data['UU']
    v = image[0].data['VV']
    u *= image[0].header['crval4']/1e3
    v *= image[0].header['crval4']/1e3

    vis = image[0].data['data']
    real = (vis[:,0,0,0,:,0,0] + vis[:,0,0,0,:,1,0])/2.
    real = np.mean(real,axis=1)
    imag = (vis[:,0,0,0,:,0,1] + vis[:,0,0,0,:,1,1])/2.
    imag = np.mean(imag,axis=1)

    amp = np.sqrt(real**2 + imag**2)
    ampsb = amp.squeeze()
    uvdist = np.sqrt(u**2 + v**2)

###### for binning!

    newuvdist = []
    newreal = []
    newimag = []
    newamp = []
    rstdevs = []
    istdevs = []
    rNs = []
    iNs = []
    amperrssb = []

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

        amperrsb = np.sqrt((rstdev**2/rN)*(realav**2/(realav**2+imagav**2))+(istdev**2/iN)*(imagav**2/(realav**2+imagav**2)))
        amperrssb.append(amperrsb)

    newuvdist = np.asarray(newuvdist)
    newuvdist = newuvdist[~np.isnan(newuvdist)]
    r = np.asarray(newreal)
    r = r[~np.isnan(r)]
    i = np.asarray(newimag)
    i = i[~np.isnan(i)]
    newampsb = np.sqrt(r**2 + i**2)
    amperrssb = np.asarray(amperrssb)
    amperrssb = amperrssb[~np.isnan(amperrssb)]

###### long baselines

    image = fits.open(lbfile)
    u = image[0].data['UU']
    v = image[0].data['VV']
    u *= image[0].header['crval4']/1e3
    v *= image[0].header['crval4']/1e3

    vis = image[0].data['data']
    real = (vis[:,0,0,0,:,0,0] + vis[:,0,0,0,:,1,0])/2.
    real = np.mean(real,axis=1)
    imag = (vis[:,0,0,0,:,0,1] + vis[:,0,0,0,:,1,1])/2.
    imag = np.mean(imag,axis=1)

    amp = np.sqrt(real**2 + imag**2)
    amplb = amp.squeeze()
    uvdist = np.sqrt(u**2 + v**2)

###### for binning!

    newuvdist = []
    newreal = []
    newimag = []
    newamp = []
    rstdevs = []
    istdevs = []
    rNs = []
    iNs = []
    amperrslb = []

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

        amperrlb = np.sqrt((rstdev**2/rN)*(realav**2/(realav**2+imagav**2))+(istdev**2/iN)*(imagav**2/(realav**2+imagav**2)))
        amperrslb.append(amperrlb)

    newuvdist = np.asarray(newuvdist)
    newuvdist = newuvdist[~np.isnan(newuvdist)]
    r = np.asarray(newreal)
    r = r[~np.isnan(r)]
    i = np.asarray(newimag)
    i = i[~np.isnan(i)]
    newamplb = np.sqrt(r**2 + i**2)
    amperrslb = np.asarray(amperrslb)
    amperrslb = amperrslb[~np.isnan(amperrslb)]

###### ratio stuff

    print newampsb.shape
    print newamplb.shape
    ampquot = newampsb/newamplb
    ampquoterr = ampquot*np.sqrt((amperrssb**2/newampsb**2)+(amperrslb**2/newamplb**2))
    plt.errorbar(newuvdist,ampquot,yerr=ampquoterr,fmt='.')
    plt.xlabel('Distance from center of uv-plane (klambda)')
    plt.axhline(1)
    plt.ylabel('Short-baseline amplitude divided by long-baseline amplitude')
    plt.title('Ratio of short/long baseline amplitudes versus uv-distance')
    plt.grid(True)
    plt.show()

    x = newuvdist
    y = ampquot
    sigma = ampquoterr
    ymodel = 1

    def chisqr(y,ymodel,sigma):
        chisqr = np.sum((y - ymodel)**2/sigma**2)
        return chisqr

    def redchisqr(y,ymodel,sigma):
        chisqr = np.sum((y - ymodel)**2/sigma**2)
        nu = len(sigma)
###### nu should be the number of data points over the number of parameters.
        return chisqr/nu

    print "chi squared =",chisqr(y,ymodel,sigma)
    print "reduced chi squared =",redchisqr(y,ymodel,sigma)
