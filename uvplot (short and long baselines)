def amp_uvdist(sbfile,lbfile,binsize):
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
#    plt.plot(u,v,'.')

    vis = image[0].data['data']
    real = (vis[:,0,0,0,:,0,0] + vis[:,0,0,0,:,1,0])/2.
    imag = (vis[:,0,0,0,:,0,1] + vis[:,0,0,0,:,1,1])/2.

    amp = np.sqrt(real**2 + imag**2)
    ampsb = amp.squeeze()
    uvdist = np.sqrt(u**2 + v**2)
#    plt.plot(uvdist,amp,'.')

###### for binning! 

    newuvdist = []
    newreal = []
    newimag = []
    newamp = []
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

    newuvdist = np.asarray(newuvdist)
    newuvdist = newuvdist[~np.isnan(newuvdist)]
    r = np.asarray(newreal)
    r = r[~np.isnan(r)]
    i = np.asarray(newimag)
    i = i[~np.isnan(i)]
    newamp = np.sqrt(r**2 + i**2)
    rstdev = np.std(r)
    istdev = np.std(i)
    rN = len(r)
    iN = len(i)

    amperrsb = np.sqrt((rstdev**2/rN)*(r**2/(r**2+i**2))+(istdev**2/iN)*(i**2/(r**2+i**2)))
#    plt.plot(newuvdist,newamp,'.')
    plt.errorbar(newuvdist,newamp,yerr=amperrsb,fmt='.')

###### long baselines

    image = fits.open(lbfile)
    u = image[0].data['UU']
    v = image[0].data['VV']
    u *= image[0].header['crval4']/1e3
    v *= image[0].header['crval4']/1e3
#    plt.plot(u,v,'.')

    vis = image[0].data['data']
    real = (vis[:,0,0,0,:,0,0] + vis[:,0,0,0,:,1,0])/2.
    imag = (vis[:,0,0,0,:,0,1] + vis[:,0,0,0,:,1,1])/2.

    amp = np.sqrt(real**2 + imag**2)
    amplb = amp.squeeze()
    uvdist = np.sqrt(u**2 + v**2)
#    plt.plot(uvdist,amp,'.')

###### for binning!

    newuvdist = []
    newreal = []
    newimag = []
    newamp = []
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

    newuvdist = np.asarray(newuvdist)
    newuvdist = newuvdist[~np.isnan(newuvdist)]
    r = np.asarray(newreal)
    r = r[~np.isnan(r)]
    i = np.asarray(newimag)
    i = i[~np.isnan(i)]
    newamp = np.sqrt(r**2 + i**2)
    rstdev = np.std(r)
    istdev = np.std(i)
    rN = len(r)
    iN = len(i)

    amperrlb = np.sqrt((rstdev**2/rN)*(r**2/(r**2+i**2))+(istdev**2/iN)*(i**2/(r**2+i**2)))
#    plt.plot(newuvdist,newamp,'.')
    plt.errorbar(newuvdist,newamp,yerr=amperrlb,fmt='.')

###### error things

    quotient = ampsb/amplb
    print quotient
