import numpy as np
import matplotlib.pyplot as plt

def chisqr(y,ymodel,sigma):
    chisqr = np.sum((y - ymodel)**2/sigma**2)
    return chisqr
def mcmc(fname,iterations):
    fakedata = np.recfromtxt(fname)
    x = fakedata[:,0]
    y = fakedata[:,1]
    sigmay = fakedata[:,2]
    bchain = []
    mchain = []
    baccepted = []
    maccepted = []

    b = np.random.uniform(-3,3)
    m = np.random.uniform(-3,3)
    sigmab = .85
    sigmam = .15
    n = 0
    for i in range(iterations):
        bprime = np.random.normal(loc=b,scale=sigmab)
        ymodel = m * x + bprime
        chi2 = chisqr(y,ymodel,sigmay)
        ymodeln = m * x + b
        chi2n = chisqr(y,ymodeln,sigmay)
        exponent = -(chi2-chi2n)/2.
        ratio = np.exp(exponent)
        u = np.random.rand()
        if u <= ratio:
            b = bprime
            bchain.append(bprime)
            baccepted.append(bprime)
        elif u > ratio:
            b = b
            bchain.append(b)
        mprime = np.random.normal(loc=m,scale=sigmam)
        ymodel = mprime * x + b
        chi2 = chisqr(y,ymodel,sigmay)
        ymodeln = m * x + b
        chi2n = chisqr(y,ymodeln,sigmay)
        exponent = -(chi2-chi2n)/2.
        ratio = np.exp(exponent)
        u = np.random.rand()
        if u <= ratio:
            m = mprime
            mchain.append(mprime)
            maccepted.append(mprime)
        elif u > ratio:
            m = m
            mchain.append(m)
        n += 1

    baccept = (float(len(baccepted))/float(len(bchain)))*100.
    maccept = (float(len(maccepted))/float(len(mchain)))*100.
    print "The acceptance percentage for b is {}%.".format(baccept)
    print "The acceptance percentage for m is {}%.".format(maccept)

    nchain = range(0,n)
    ax1 = plt.subplot(211)
    plt.plot(nchain,bchain,linewidth=1.0)
    plt.setp(ax1.get_xticklabels(),visible=False)
    ax1.set_title('Chains over time')
    ax1.set_ylabel('b value')
    ax1.grid(True)
    ax2 = plt.subplot(212,sharex=ax1)
    plt.plot(nchain,mchain,linewidth=1.0)
    plt.setp(ax2.get_xticklabels(),fontsize=8)
    ax2.set_ylabel('m value')
    ax2.set_xlabel('Iteration')
    ax2.grid(True)
    plt.show()
