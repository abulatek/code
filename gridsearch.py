import numpy as np

def grid_search(fname):
    fakedata = np.recfromtxt(fname)
    x = fakedata[:,0]
    y = fakedata[:,1]
    sigma = fakedata[:,2]

    mrange = np.arange(-2,2,0.01)
    brange = np.arange(-2,2,0.01)

    def chi2(y,ymodel,sigma):
        chi2 = np.sum((y - ymodel)**2/sigma**2)
        return chi2

    def redchi2(y,ymodel,sigma):
        chi2 = np.sum((y - ymodel)**2/sigma**2)
        nu = 2
        return chi2/nu

    values = []
    for m in mrange:
        for b in brange:
            ymodel = m * x + b
            k = chi2(y,ymodel,sigma)
            values.append(k)
