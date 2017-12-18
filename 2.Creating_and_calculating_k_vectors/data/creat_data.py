#!/usr/bin/env python
#########################################
##Code authors
#Chong Han, chan104@illinois.edu 
#########################################

import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    
    #force difference from dummy.force
    f_error =[]
    file=open('SG7k_r1.15_eam1.force','r')
    line='1'
    i=0
    ave=0
    while line != '':
        if (str(i)+': ') in line:
            #print line.split()
            for j in range(64*3):
                f=float(line.split()[4])
                f0=float(line.split()[5])
                f_error.append(f0-f)
                line=file.readline()
            #print i,ave
            i+=1
            ave=0
        line=file.readline()
    file.close()


    #Compute energy difference between predicted and real values
    e_error = np.loadtxt('SG7k_r1.15_eam1.energy')[:,6]
    e_error = -1*e_error # e0-e instead of e-e0
    e0= np.loadtxt('SG7k_r1.15_eam1.energy')[:,4]

    configs = 49
    #positions and box length
    #too lazy to write readlines, just np.loadtxt to get positions then assign box lengths
    pos=np.loadtxt('SG_potfit_temp7000K_rs1.15.dat') #Reads from potfit input file

    box=[7.4152 for i in range(configs)]

    file=open('SGpos.dat','w')

    for i in range(configs):
        print >> file, box[i],
        for j in range(64):
            print >> file,pos[i*64+j,1],pos[i*64+j,2],pos[i*64+j,3],
        print >> file, ''
    file.close()
    file=open('SGdft_eam_difference.dat','w')
    for i in range(configs):
    	print >>file, e_error[i],
    	for j in range(64*3):
    		print >>file, f_error[i*64*3+j],
    	print >>file,''
    file.close()

    file = open('SG_dft7Krs1.15.dat','w')
    for i in range(configs):
        print >>file,e0[i],
        for j in range(64):
            print >>file,pos[i*64+j,4],pos[i*64+j,5],pos[i*64+j,6],
        print >>file,''
    file.close()
