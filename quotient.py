###### In order to determine values to use for rangemin and rangemax,
###### use amp_uvdist.py to plot amp vs. uvdist and bound overlap.

def quotient(sbfile,lbfile,binsize,rangemin,rangemax):
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import numpy as np
    plt.clf()

###### short baselines

    image = fits.open(sbfile)
    u = image[0].data['UU']
    v = image[0].data['VV']
    u *= image[0].header['crval4']/1e3
    v *= image[0].header['crval4']/1e3

    vis = image[0].data['data']
    real = (vis[:,0,0,0,:,0,0] + vis[:,0,0,0,:,1,0])/2.
    imag = (vis[:,0,0,0,:,0,1] + vis[:,0,0,0,:,1,1])/2.

    amp = np.sqrt(real**2 + imag**2)
    ampsb = amp.squeeze()
    uvdist = np.sqrt(u**2 + v**2)

###### for binning!

    newuvdist = []
    newreal = []
    newimag = []
    newamp = []

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

    newuvdist = np.asarray(newuvdist)
    newuvdist = newuvdist[~np.isnan(newuvdist)]
    r = np.asarray(newreal)
    r = r[~np.isnan(r)]
    i = np.asarray(newimag)
    i = i[~np.isnan(i)]
    newampsb = np.sqrt(r**2 + i**2)
    rstdev = np.std(r)
    istdev = np.std(i)
    rN = len(r)
    iN = len(i)

    amperrsb = np.sqrt((rstdev**2/rN)*(r**2/(r**2+i**2))+(istdev**2/iN)*(i**2/(r**2+i**2)))

###### long baselines

    image = fits.open(lbfile)
    u = image[0].data['UU']
    v = image[0].data['VV']
    u *= image[0].header['crval4']/1e3
    v *= image[0].header['crval4']/1e3

    vis = image[0].data['data']
    real = (vis[:,0,0,0,:,0,0] + vis[:,0,0,0,:,1,0])/2.
    imag = (vis[:,0,0,0,:,0,1] + vis[:,0,0,0,:,1,1])/2.

    amp = np.sqrt(real**2 + imag**2)
    amplb = amp.squeeze()
    uvdist = np.sqrt(u**2 + v**2)

###### for binning!

    newuvdist = []
    newreal = []
    newimag = []
    newamp = []

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

    newuvdist = np.asarray(newuvdist)
    newuvdist = newuvdist[~np.isnan(newuvdist)]
    r = np.asarray(newreal)
    r = r[~np.isnan(r)]
    i = np.asarray(newimag)
    i = i[~np.isnan(i)]
    newamplb = np.sqrt(r**2 + i**2)
    rstdev = np.std(r)
    istdev = np.std(i)
    rN = len(r)
    iN = len(i)

    amperrlb = np.sqrt((rstdev**2/rN)*(r**2/(r**2+i**2))+(istdev**2/iN)*(i**2/(r**2+i**2)))

###### ratio stuff

    ampquot = newampsb/newamplb
    ampquoterr = ampquot*np.sqrt((amperrsb**2/newampsb**2)+(amperrlb**2/newamplb**2))
    plt.errorbar(newuvdist,ampquot,yerr=ampquoterr,fmt='.')
    plt.xlabel('Distance from center of uv-plane (klambda)')
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

    print chisqr(y,ymodel,sigma)
