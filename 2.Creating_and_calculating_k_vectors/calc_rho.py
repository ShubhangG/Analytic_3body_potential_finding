#!/usr/bin/env python
#########################################
##Code authors
#Chong Han, chan104@illinois.edu 
#########################################

import numpy as np
import math

if __name__ == '__main__':
    
    cutoff = 4
    total_configurations = 49

    #make k vectors, since rho_k and rho_(-k) are complex conjugate, only one need to be calculated
    k_vectors=[]
    for i in range(0,cutoff+1):
        for j in range(-cutoff,cutoff+1):
            for k in range(-cutoff,cutoff+1):
                if i*i+j*j+k*k<=cutoff*cutoff:
                    if (i>0) or (i==0 and j>0) or (i==0 and j==0 and k>=0): #make sure just count the positive k
                        k_vectors.append([i,j,k])

    data = np.loadtxt('data/SGpos.dat')

    file = open('k_vectors.dat','w')
    for i in range(len(k_vectors)):
        print >>file, k_vectors[i][0],k_vectors[i][1],k_vectors[i][2]
    file.close()


    #config_index is line 
    file2 = open('rho.dat','w')
    def calc(config_index):
        L=data[config_index,0] # box length
        pos = data[config_index,1:64*3+1].reshape(-1,3)
        k_cos=[0.0 for i in range(len(k_vectors))] #real part
        k_sin=[0.0 for i in range(len(k_vectors))] #imaginary part
        for i in range(len(k_vectors)):
            k=[math.pi*2/L*k_vectors[i][0],math.pi*2/L*k_vectors[i][1],math.pi*2/L*k_vectors[i][2]]
            for j in range(64):
                k_cos[i]+=math.cos(k[0]*pos[j,0]+k[1]*pos[j,1]+k[2]*pos[j,2])
                k_sin[i]+=math.sin(k[0]*pos[j,0]+k[1]*pos[j,1]+k[2]*pos[j,2])
            print >>file2, k_cos[i],k_sin[i],
        print >>file2, ''

    for i in range(total_configurations):
        calc(i)
