#!/usr/bin/env python
#########################################
##Code authors
#Shubhang Goswami, sgoswam3@illinois.edu 
#########################################

import pandas as pd 
import numpy as np 
from itertools import product
import matplotlib.pyplot as plt
import pickle as pk 

			

def distance(posi,posj,box_length):
    """ return the distance between particle i and j according to the minimum image convention. """

    dist = 0.0
    if len(posi) != len(posj):		
        print("Incorrect dimensions!")
        return
    ndim = len(posi)	#Number of dimensions, usually 3- x,y,z
    A = np.zeros(ndim)
    for i in range(ndim):
        d = abs(posi[i]-posj[i])	#Consider periodic boundary conditions
        A[i] = min(d,box_length-d)	#Choose the shortest distance

    dist = np.linalg.norm(A)
    return dist

def legal_kvecs(maxk):
    kvecs = []
    p=np.arange(0,maxk)
    kvecs = [x for x in product(p,repeat=3)]	#Permutations of k
    # calculate a list of legal k vectors
    return np.array(kvecs)*2*np.pi/lattice_dim	

def rhok(kvec, pset):
    value = 0.0
    #computes \sum_j \exp(i * k \dot r_j) 
    for ptcl in pset:
        value=value + np.exp(1j*np.dot(kvec,ptcl))
    return value
# end def

def Sk(kvecs, pset):
    """ computes structure factor for all k vectors in kList
     and returns a list of them """
    sk_list=[]
    for kvec in kvecs:
        sk_list.append(np.real((1.0/64)*rhok(kvec,pset)*rhok(-1*kvec,pset)).tolist())
    return [i/sum(sk_list) for i in sk_list]

#################################################################################################
if __name__ == '__main__':
    #READ positions input, it would be tstep*Natoms long and 3 columns wide
    cdf = pd.read_table("positions.dat", delim_whitespace=True, header=None, usecols=[1,2,3])
    cdf.columns = ['X', 'Y', 'Z']   #Set the name for each column
    lattice_dim = 7.41516299        #The lattice constant given in bohr
    tsteps=len(cdf)/64              #Time steps, we used 64 atoms so /64

    #Convert to matrix form for easy manipulaiton
    cf=cdf.as_matrix()
    kvecs=legal_kvecs(5)	
    pos = cf
    var = np.split(pos,tsteps)	#Split into different configurations
    SK_list = []
    for row in var:
        SK_list.append(Sk(kvecs,row))

    #Calculate s(k)
    Sk = np.array(SK_list)	
    sk_list = np.mean(Sk, axis=0)
    sk_list[0]=0

    #Normalizing
    # fact = 1.0/64
    # sk_list = [v*fact for v in sk_list]

    kmags  = [np.linalg.norm(kvec) for kvec in kvecs]
    sk_arr = np.array(sk_list)


    unique_kmags = np.unique(kmags)
    unique_sk    = np.zeros(len(unique_kmags))
    for iukmag in range(len(unique_kmags)):
        kmag    = unique_kmags[iukmag]
        idx2avg = np.where(kmags==kmag)
        unique_sk[iukmag] = np.mean(sk_arr[idx2avg])

    #Pickle the file, so you don't need to calculate again and again
    pk.dump((unique_kmags, unique_sk) , open("sok.p","wb"))
    plt.plot(unique_kmags,unique_sk)
    plt.xlabel("k values")
    plt.ylabel("S(k)")
    plt.title("Plotting S(K)")
    plt.savefig('Skplot.png')
    plt.close()

################################################################################################
'''
#Time to compute g(r)
bins = {}
for pos in np.array_split(cf, tsteps):
	for i in range(len(pos)):
		for j in range(i,len(pos)):
			dist = round(distance(pos[i],pos[j],lattice_dim), 1)	
			if dist>=lattice_dim/2 or dist==0:		#binning same distance objects
				continue
			if dist in bins:
				bins[dist]+=1
		 	else:
		 		bins[dist]= 1


#Normalizing g(r)
delg=lattice_dim/(2*len(bins))
#vb=((0.2)**3-(0.1)**3)*delg**3
rho= 64/(lattice_dim**3)
factor =(4/3)*np.pi*rho*tsteps*64/2
#print lists
#bkeys = bins.keys().sort
normalised_bin = {k: v/factor for k, v in bins.iteritems()}
lists= sorted(normalised_bin.items())


x,y = zip(*lists)
z=np.zeros(len(y))
for i in range(len(y)):
	if i==0:
		z[i] = y[i]/(x[i]**3)
	else:
		z[i] = y[i]/(x[i]**3-x[i-1]**3)


pk.dump((x,z), open("bins.p","wb"))

#Plotting g(r)
plt.plot(x,z)
plt.title("Plotting g(r)")
plt.xlabel("The distance")
plt.ylabel("g(r)")
plt.savefig("gorplot.png")
plt.close()
'''