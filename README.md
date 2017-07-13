# code
This repository is evidence of my attempt to learn how Python works. I have written my own code to analyze ALMA data, and I wrote my own MCMC Metropolis-Hastings algorithm in order to be able to eventually use emcee (The MCMC Hammer) to find best-fit parameters for protoplanetary disks. This work was done at Wesleyan University in an REU through the KNAC program. My advisor, Dr. Kevin Flaherty, deserves much credit for the success of this code.

The main features of this repository include
  -   calcheck.py, which calculates the amplitude offset between two ALMA visibilities, and
  -   mcmc.py, which is an MCMC Metropolis-Hastings algorithm to calculate (linear) best-fit parameters for a data set.

Thanks much, and I hope that this might help someone at some point (but there are programs built-in to CASA and Miriad that do what these codes do!).
