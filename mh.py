import numpy as np
import matplotlib.pyplot as plt

def chisqr(y,ymodel,sigma):
    chisqr = np.sum((y - ymodel)**2/sigma**2)
    return chisqr
### The input 'iterations' is how many iterations you want for each parameter (total is then 2*iterations).
def mcmc(fname,iterations):
### Import the data.
    fakedata = np.recfromtxt(fname)
    x = fakedata[:,0]
    y = fakedata[:,1]
    sigmay = fakedata[:,2]
    mchain = []
    bchain = []
    maccepted = []
    baccepted = []

### Choose a random float to start.
    m = np.random.uniform(-3,3)
    b = np.random.uniform(-3,3)
### Choose step sizes for each parameter--a downfall of Metropolis-Hastings.
    sigmam = .06
    sigmab = .3
    n = 0
    for i in range(iterations):
### Generate a test value for m.
        mprime = np.random.normal(loc=m,scale=sigmam)
### Calculate chi2 values for the trial state and the current state.
        mymodel = mprime * x + b
        mchi2 = chisqr(y,mymodel,sigmay)
        mymodeln = m * x + b
        mchi2n = chisqr(y,mymodeln,sigmay)
### Determine the ratio given by equation 11 in Ford 2005.
        mratio = np.exp(-(mchi2-mchi2n)/2.)
### Repeat for the other parameter.
        bprime = np.random.normal(loc=b,scale=sigmab)
        bymodel = m * x + bprime
        bchi2 = chisqr(y,bymodel,sigmay)
        bymodeln = m * x + b
        bchi2n = chisqr(y,bymodeln,sigmay)
        bratio = np.exp(-(bchi2-bchi2n)/2.)
### Draw a random number, u, from a uniform distribution between 0 and 1.
        u = np.random.rand()
        if u <= mratio:
### If u <= the ratio, set the next parameter value to the test value.
            mchain.append(mprime)
            bchain.append(b)
            maccepted.append(mprime)
            n += 1
        elif u > mratio:
### If u > the ratio, keep the same start values.
            mchain.append(m)
            bchain.append(b)
            n += 1
### Repeat for the other parameter.
        v = np.random.rand()
        if v <= bratio:
            bchain.append(bprime)
            mchain.append(m)
            baccepted.append(bprime)
            n += 1
        elif v > bratio:
            mchain.append(m)
            bchain.append(b)
            n += 1
### Move onto the next value if the value is accepted.
        if u <= mratio:
            m = mprime
        if v <= bratio:
            b = bprime

### Calculate acceptance fractions.
    maccept = (float(len(maccepted))/float(len(mchain)))*100.
    baccept = (float(len(baccepted))/float(len(bchain)))*100.
    print "The acceptance percentage for m is {}%.".format(maccept)
    print "The acceptance percentage for b is {}%.".format(baccept)

### Get rid of the burn-in.
    nchain = range(0,n)
    mhist = mchain[1000:]
    bhist = bchain[1000:]
    print mhist[120]
    print bhist[120]

### Print mean and standard deviation of the two-dimensional posterior distributions.
    mmean = np.mean(mhist)
    mstddev = np.std(mhist,ddof=1)
    bmean = np.mean(bhist)
    bstddev = np.std(bhist,ddof=1)
    print "The mean of the m value distribution is {}.".format(mmean)
    print "The standard deviation of the m value distribution is {}.".format(mstddev)
    print "The mean of the b value distribution is {}.".format(bmean)
    print "The standard deviation of the b value distribution is {}.".format(bstddev)

### This plot shows the entire chains.
    plt.figure(1)
    ax1 = plt.subplot(211)
    plt.plot(nchain,mchain,linewidth=1.0)
    plt.setp(ax1.get_xticklabels(),visible=False)
    ax1.set_title('Chains over time')
    ax1.set_ylabel('m value')
    ax1.grid(True)
    ax2 = plt.subplot(212,sharex=ax1)
    plt.plot(nchain,bchain,linewidth=1.0)
    plt.setp(ax2.get_xticklabels(),fontsize=8)
    ax2.set_ylabel('b value')
    ax2.set_xlabel('Iteration')
    ax2.grid(True)
    plt.show()

### This is a plot which only samples the posterior distribution (no burn-in).
    plt.figure(2)
    plt.plot(mhist,bhist)
    plt.grid(True)
    plt.ylabel("b")
    plt.xlabel("m")
    plt.title("b and m in Parameter Space (after burn-in)")
    plt.show()

### This plot is the two-dimensional posterior distributions.
    plt.figure(3)
    ax1 = plt.subplot(211)
    plt.hist(mhist,bins=50)
    plt.setp(ax1.get_xticklabels(),fontsize=8)
    ax1.set_title('Two-dimensional Posterior Distributions')
    ax1.set_xlabel('m')
    ax1.set_ylabel('Frequency')
    ax1.grid(True)
    ax2 = plt.subplot(212,sharey=ax1)
    plt.hist(bhist,bins=50)
    plt.setp(ax2.get_xticklabels(),fontsize=8)
    ax2.set_xlabel('b')
    ax2.set_ylabel('Frequency')
    ax2.grid(True)
    plt.tight_layout(pad=2.0, w_pad=1.0, h_pad=1.0)
    plt.show()
