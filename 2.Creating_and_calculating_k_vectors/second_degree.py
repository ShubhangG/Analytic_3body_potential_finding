import numpy as np
import math
import cmath

cutoff = 4
total_configurations = 41
data = np.loadtxt('data/SGpos.dat')
k_vectors=np.loadtxt('k_vectors.dat')
rho_k_data=np.loadtxt('rho.dat')
output = open('second_degree.dat','w')
mapping_k=[0 for i in range(cutoff*cutoff+1)]
n_k_squared=0

def get_rho_k(config_index,k_index):
    if k_index<0:
        return -k_vectors[-k_index,:],complex(rho_k_data[config_index,-k_index*2],-rho_k_data[config_index,(-k_index)*2+1])
    else:
        return k_vectors[k_index,:],complex(rho_k_data[config_index,k_index*2],rho_k_data[config_index,k_index*2+1])

def rho_k_derivative(config_index,k_index,iatom,dimension):
    L=data[config_index,0] # box length
    pos = data[config_index,1:64*3+1].reshape(-1,3)
    k,rho_k=get_rho_k(config_index,k_index)    
    k_scaled=k*math.pi*2/L
    #i*k_i*e^(i*vec_k*vec_r)
    return complex(0,1)*k_scaled[dimension]*cmath.exp(1j*(k_scaled[0]*pos[iatom,0]+k_scaled[1]*pos[iatom,1]+k_scaled[2]*pos[iatom,2]))


def create_k_map():
    global n_k_squared
    flag=[0 for i in range(cutoff*cutoff+1)]
    for i in range(len(k_vectors)):
        k=k_vectors[i,:]
        k_squared=int(k[0]*k[0]+k[1]*k[1]+k[2]*k[2])
        flag[k_squared]+=1
    for i in range(cutoff*cutoff+1):
        if flag[i]>0:
            mapping_k[i]=n_k_squared
            n_k_squared+=1

create_k_map()

def map_k(k):
    k_squared=int(k[0]*k[0]+k[1]*k[1]+k[2]*k[2])
    return mapping_k[k_squared]


def calc(config_index):
    L=data[config_index,0] # box length
    rho_k_squared=[0.0 for i in range(n_k_squared)]

    pos = data[config_index,1:64*3+1].reshape(-1,3)
    for j in range(len(k_vectors)):
        k,rho_k=get_rho_k(config_index,j)
        mapped_k=map_k(k)
        rho_k_squared[mapped_k]+=abs(rho_k)*abs(rho_k)

    print >>output, rho_k_squared[0],
    #for k!=0 the rho_k should be multiply by 2
    for i in range(1,n_k_squared):
        print >>output, rho_k_squared[i]*2,

    #taking derivative
    for iatom in range(64):
        for dimension in range(3):
            rho_k_squared=[0.0 for i in range(n_k_squared)] #derivative  storage
            for j in range(len(k_vectors)):
                k,rho_k=get_rho_k(config_index,j)
                mapped_k=map_k(k)
                '''
                old method
                k_scaled=[math.pi*2/L*k_vectors[j][0],math.pi*2/L*k_vectors[j][1],math.pi*2/L*k_vectors[j][2]]
                #real part of -i*2*k_i*e^(-i*vec_k*vec_r)*rho_k
                rho_k_squared[mapped_k]+=2*k_scaled[dimension]*(math.cos(
                    -k_scaled[0]*pos[iatom,0]-k_scaled[1]*pos[iatom,1]-k_scaled[2]*pos[iatom,2])*\
                rho_k[config_index,2*j+1]+math.sin(
                    -k_scaled[0]*pos[iatom,0]-k_scaled[1]*pos[iatom,1]-k_scaled[2]*pos[iatom,2])*\
                rho_k[config_index,2*j])
                '''
                #real part of -i*2*k_i*e^(-i*vec_k*vec_r)*rho_k
                rho_k_squared[mapped_k]+=2*(rho_k_derivative(config_index,-j,iatom,dimension)*rho_k).real

            for i in range(n_k_squared):
                print >>output, rho_k_squared[i],
    print >>output, ''

#print mapping_k
for  i in range(total_configurations):
    print "running configuration ",i
    calc(i)
