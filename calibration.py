from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
###### Both angle inputs should be in degrees.
###### My values: V4046     PA = 76     incl = 33
######            MWC480    PA = -34    incl = 36

def calibration(sbfile,lbfile,pos_angle,incl,binsize):
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

###### deprojected uv-distances
        phi = np.radians(np.degrees(np.arctan2(v,u)) - pos_angle)
        da = uvdist*np.sin(phi)
        db = uvdist*np.cos(phi)*np.cos(np.radians(incl))
        deprojuvdist = np.sqrt(da**2 + db**2)
#        plt.plot(da,db,'.')
#        plt.plot(deprojuvdist,amp,'.')

###### for binning
        newuvdist = []
        newreal = []
        newimag = []
        newamp = []
        rstdevs = []
        istdevs = []
        rNs = []
        iNs = []
        amperrs = []
        rangemin = 10
        rangemax = 900

        for minimum in range(rangemin,rangemax,binsize):
            maximum = minimum + binsize
            uvsubarray = (deprojuvdist > minimum) & (deprojuvdist < maximum)
            addition = (maximum-minimum)/2.
            median = minimum + addition
            newuvdist.append(median)

            realsubarray = (real[uvsubarray] != 0)
            real2 = real[uvsubarray][realsubarray]
            realav = np.nanmean(real2)
            newreal.append(realav)

            imagsubarray = (imag[uvsubarray] != 0)
            imag2 = imag[uvsubarray][imagsubarray]
            imagav = np.nanmean(imag2)
            newimag.append(imagav)

            rstdev = np.nanstd(real2)
            istdev = np.nanstd(imag2)
            rlist = list(real2[~np.isnan(real2)])
            ilist = list(imag2[~np.isnan(imag2)])
            rN = len(rlist)
            iN = len(ilist)
            denom = realav**2+imagav**2
            amperr = np.sqrt((rstdev**2/rN)*(realav**2/denom)+(istdev**2/iN)*(imagav**2/denom))
            amperrs.append(amperr)

        newuvdist = np.asarray(newuvdist)
        newuvdist = newuvdist[~np.isnan(newuvdist)]
        newuvdist = [i for i in newuvdist if i <= max(deprojuvdist)]
        rilength = len(newuvdist)
        newuvdist = np.asarray(newuvdist)

        r = newreal[:rilength]
        r = np.asarray(r)
        i = newimag[:rilength]
        i = np.asarray(i)
        remove = (np.isnan(r)) & (np.isnan(i))
        newuvdist = newuvdist[~remove]
        newamp = np.sqrt(r**2 + i**2)
        newamp = newamp[~remove]
        amperrs = amperrs[:rilength]
        amperrs = np.asarray(amperrs)
        amperrs = amperrs[~remove]

        plt.errorbar(newuvdist,newamp,yerr=amperrs,fmt='.')
        plt.xlabel('Deprojected uv-distance (klambda)')
        plt.ylabel('Amplitude (Jy)')
        plt.title('Amplitude versus deprojected uv-distance')
        plt.grid(True)
        plt.show()
