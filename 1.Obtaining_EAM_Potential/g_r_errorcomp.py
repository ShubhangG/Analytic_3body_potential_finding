#!/usr/bin/env python
#########################################
##Code authors
##Paul Yang, yyang173@illinois.edu
#Shubhang Goswami, sgoswam3@illinois.edu 
#########################################

import pandas as pd 
import numpy as np
import pickle as pk 
import matplotlib.pyplot as plt
import collections
    
def distance(posi,posj,box_length):
    """ return the distance between particle i and j according to the minimum image convention. """

    dist = 0.0
    if len(posi) != len(posj):
        print("Incorrect dimensions!")
        return
    ndim = len(posi)
    A = np.zeros(ndim)
    for i in range(ndim):
        d = abs(posi[i]-posj[i])
        A[i] = min(d,box_length-d)

    dist = np.linalg.norm(A)
    return dist

def g(r):
    bins = collections.defaultdict(list)
    for pos in np.array_split(r, tsteps):
        for i in range(len(pos)):
            for j in range(i,len(pos)):
                dist = round(distance(pos[i],pos[j],lattice_dim), 1)    
                if dist>=lattice_dim/2 or dist==0:      #binning same distance objects
                    continue
                if dist in bins:
                    bins[dist][i]+=1
                else:
                    bins[dist]= np.zeros(tsteps)
    return bins
def corr(trace):
    """ calculate the autocorrelation of a trace of scalar data
    pre:  trace should be a 1D iterable array of floating point numbers
    post: return the autocorrelation of this trace of scalars 
    """

    correlation = 1.0
    su =0.0
    n = len(trace)
    for k in range(1,n):
        cc = compcorr(trace,k)
        if cc<=0:
            break
        su = su + cc

    correlation = 1 + 2*su
    return correlation

def compcorr(trace, k):
    R = 0.0
    n= len(trace)
    mu = np.mean(trace)
    sigma = np.std(trace)
    sigma = sigma**2
    for t in range(n-k):
        R = R + (trace[t] - mu)*(trace[t+k] - mu)
    
    R = R/((n-k)*sigma)
    return R

# end def corr

def error(trace):
    """ calculate the standard error of a trace of scalar data
    for uncorrelated data, this should match np.std(trace)/np.sqrt(len(trace))
    pre:  trace should be a 1D iterable array of floating point numbers
    post: return the standard error of this trace of scalars 
    """

    #stderr = 0.0
    N_eff = len(trace)/corr(trace)
    pk.dump(N_eff, open("N_eff.pk","wb"))
    sigma = np.std(trace)
    stderr = sigma/np.sqrt(N_eff)
    return stderr


if __name__ == '__main__':
    cdf = pd.read_table("positions.dat", delim_whitespace=True, header=None, usecols=[1,2,3]) #Read positions
    cdf.columns = ['X', 'Y', 'Z']
    lattice_dim = 7.4152                                        #box length L
    tsteps=len(cdf)/64
    cf=cdf.as_matrix()
    rho= 64/(lattice_dim**3)
    factor =(4/3)*np.pi*rho*64*1250
    b = {k: v.sum()/factor for k,v in g(cf).iteritems()}        
    e = {k: error(v)/factor for k,v in g(cf).iteritems()}       
    err = sorted(e.items())
    lists = sorted(b.items())

    x,y = zip(*lists)
    a, yerr = zip(*err)
    z=np.zeros(len(y))
    zerr = np.zeros(len(yerr))
    for i in range(len(y)):
        if i==0:
            z[i] = y[i]/(x[i]**3)
            zerr[i] = yerr[i]/(a[i]**3) 
        else:
            z[i] = y[i]/(x[i]**3-x[i-1]**3)
            zerr[i] = yerr[i]/(a[i]**3-a[i-1]**3)

    plt.errorbar(x,z, yerr=zerr)
    plt.savefig("g(r)_errorbars.png")
    plt.show()