import numpy as np
import math

cutoff = 4
k_vectors=np.loadtxt('k_vectors.dat')
data = np.loadtxt('data/pos.dat')

mapping_k=[0 for i in range(cutoff*cutoff+1)]
n_k_squared=0
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



def calc(shift):
    L=data[0,0] # box length
    pos = data[0,1:64*3+1].reshape(-1,3)
    pos[0][0]+=shift
    k_cos=[0.0 for i in range(len(k_vectors))] #real part
    k_sin=[0.0 for i in range(len(k_vectors))] #imaginary part
    for i in range(len(k_vectors)):
        k=[math.pi*2/L*k_vectors[i][0],math.pi*2/L*k_vectors[i][1],math.pi*2/L*k_vectors[i][2]]
        for j in range(64):
            k_cos[i]+=math.cos(k[0]*pos[j,0]+k[1]*pos[j,1]+k[2]*pos[j,2])
            k_sin[i]+=math.sin(k[0]*pos[j,0]+k[1]*pos[j,1]+k[2]*pos[j,2])

    rho_k_squared=[0.0 for i in range(n_k_squared)]

    for j in range(len(k_vectors)):
        k=k_vectors[j,:]
        mapped_k=map_k(k)
        if j==0:
            rho_k_squared[mapped_k]+=k_cos[j]*k_cos[j]+k_sin[j]*k_sin[j]
        else:
            rho_k_squared[mapped_k]+=2*(k_cos[j]*k_cos[j]+k_sin[j]*k_sin[j])
    return np.array(rho_k_squared)

delta=1e-6

print (calc(delta)-calc(-delta))/(2*delta)
print calc(0)
