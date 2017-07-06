import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def best_fit(fname):
    plt.clf()

    fakedata = np.recfromtxt(fname)
    x = fakedata[:,0]
    y = fakedata[:,1]
    uncertainty = fakedata[:,2]

    def line(x,m,b):
        return m * x + b

    params, errors = curve_fit(line,x,y,sigma=uncertainty)

    print "m =",params[0],"+/-",errors[0,0]**0.5
    print "b =",params[1],"+/-",errors[1,1]**0.5
    
    plt.errorbar(x,y,yerr=uncertainty,fmt=None)
    xfine = np.linspace(0.,10.,10)
    plt.plot(xfine,line(xfine,params[0],params[1]),'r-')
# Above, instead of 'params[0],params[1]', you can do *params.
    plt.show()
