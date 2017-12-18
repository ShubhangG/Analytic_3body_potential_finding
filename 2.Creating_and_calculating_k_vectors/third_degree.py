# python third_degree.py > 
#!/usr/bin/env python
#########################################
##Code authors
#Chong Han, chan104@illinois.edu 
#########################################


import numpy as np
import math
import cmath

cutoff = 4
total_configurations = 41
data = np.loadtxt('data/SGpos.dat')
k_vectors=np.loadtxt('k_vectors.dat')
rho_k_data=np.loadtxt('rho.dat')

output = open('third_degree.dat','w')

def run(config_index):
    rho_k_tripled={}
    L=data[config_index,0] # box length
    pos = data[config_index,1:64*3+1].reshape(-1,3)

    def get_k(k_index):
        if k_index<0:
            return -k_vectors[-k_index,:]
        else:
            return k_vectors[k_index,:]

    def get_rho_k(k_index):
        if k_index<0:
            return complex(rho_k_data[config_index,-k_index*2],-rho_k_data[config_index,-k_index*2+1])
        else:
            return complex(rho_k_data[config_index,k_index*2],rho_k_data[config_index,k_index*2+1])

    def rho_k_derivative(k_index,iatom,dimension):
        k=get_k(k_index)    
        k_scaled=k*math.pi*2/L
        #i*k_i*e^(i*vec_k*vec_r)
        return complex(0,1)*k_scaled[dimension]*cmath.exp(1j*(k_scaled[0]*pos[iatom,0]+k_scaled[1]*pos[iatom,1]+k_scaled[2]*pos[iatom,2]))

    def calc_one_k(i,j,m):
        ret=np.array([0.0 for t in range(64*3+1)])
        rho_k1=get_rho_k(i)
        rho_k2=get_rho_k(j)
        rho_k3=get_rho_k(m)
        ret[0]=(rho_k1*rho_k2*rho_k3).real
        #taking derivatives
        for iatom in range(64):
            for dimension in range(3):
                drho_k1=rho_k_derivative(i,iatom,dimension)
                drho_k2=rho_k_derivative(j,iatom,dimension)
                drho_k3=rho_k_derivative(m,iatom,dimension)
                ret[iatom*3+dimension+1]=(drho_k1*rho_k2*rho_k3+rho_k1*drho_k2*rho_k3+rho_k1*rho_k2*drho_k3).real
        return ret

    def calc():
        total=0
        l=len(k_vectors)
        #make sure k1<=k2<=k3, so 1/6 
        for i in range(-l+1,l):
            k1=get_k(i)
            for j in range(i,l):
                k2=get_k(j)
                for m in range(j,l):
                    k3=get_k(m)
                    #satisfy k1+k2+k3==0
                    if ((k1+k2+k3)[0]==0 and (k1+k2+k3)[1]==0 and(k1+k2+k3)[2]==0):
                        k1_squared=int(k1[0]*k1[0]+k1[1]*k1[1]+k1[2]*k1[2])
                        k2_squared=int(k2[0]*k2[0]+k2[1]*k2[1]+k2[2]*k2[2])
                        k3_squared=int(k3[0]*k3[0]+k3[1]*k3[1]+k3[2]*k3[2])
                        key=tuple(sorted([k1_squared,k2_squared,k3_squared]))
                        #using index instead of actual values
                        value=calc_one_k(i,j,m)
                        if key in rho_k_tripled:
                            rho_k_tripled[key]+=value
                        else:
                            rho_k_tripled[key]=value

    def printout():
        #each rho_k_tripled has 1(energy)+64*3(forces)
        keys=sorted(rho_k_tripled)
        for i in range(64*3+1):
            for j in keys:
                print >>output, rho_k_tripled[j][i],
        print >>output,''
        # for j in keys:
        #     print j
    calc()
    printout()

for i in range(total_configurations):
    print "running configuration ",i
    run(i)