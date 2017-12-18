import math
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from itertools import product


energy_scale=30

n_parameters=15
def create_r(maxr):
    rvec = []
    dr = 0.2
    p = np.arange(0,maxr,dr) - 1.0*maxr/2
    rvec = [x for x in product(p,repeat=3)]
    return np.asarray(rvec)

def avg_r(r,V):
    #print r
    rmag = np.linalg.norm(r,axis=1)
    rids = map(int,map(round,rmag*100))
    unique_rid = np.unique(rids)
    unique_rmags = np.zeros(len(unique_rid))
    unique_V    = np.zeros(len(unique_rid),dtype=complex)
    for iurmag in range(len(unique_rid)):
        rid = unique_rid[iurmag] # shell ID
        sel = rids == rid # select this shell
        unique_V[iurmag] = np.mean(V[sel])
        unique_rmags[iurmag] = float(rid)/100

    return unique_rmags, unique_V

def run(start,N,flag):
    print 'start:',start,'N:',N,'flag:',flag
    #since there are 15 parameters to fit
    a=np.loadtxt('second_degree.dat')
    dat=np.loadtxt('data/SGdft_eam_difference.dat')
    E0 = np.loadtxt('data/SG7k_r1.15_eam1.energy',usecols=[4])
    F0 = np.loadtxt('data/SG7k_r1.15_eam1.force',usecols=[5])
    #print dat[0]
    #dat=np.loadtxt('data/dft.dat')
    if flag:
        select_forces=True
        select_energy=True
        neglect_constant=False
    else:
        select_forces=True
        select_energy=False
        neglect_constant=True

    #select forces
    if not select_energy:
        a=a[:,n_parameters:]
        dat=dat[:,1:]
    if not select_forces:
        a=a[:,:n_parameters]
        dat=dat[:,:1]
    if select_energy:
        a[:,:n_parameters]=energy_scale*a[:,:n_parameters]
        dat[:,:1]=energy_scale*dat[:,:1]

    a=a[start:start+N].reshape(-1,n_parameters)
    b=dat[start:start+N].reshape(-1,1)
    if neglect_constant:
        a=a[:,1:]
    print '#matrix A size:',len(a),len(a[0])
    print '#b size:',len(b)
    print a

    x, resid, rank, sigma = linalg.lstsq(a, b) #a is the A_alphabeta

    print '#rank:',rank

    for i in range(len(x)):
        print x[i,0]
    #plt.plot(x,'bo')
    #plt.show()
    #print '#residue',resid
    s=0
    for i in range(len(b)):
        s+=b[i]*b[i]
    print '#original ',s[0]
    # print b
    #print 'rank',rank
    #print 'sigma',sigma
    sum=0
    print len(b)

    residue= b-np.dot(a,x)
    #print residue
    
    EminEo = np.square(residue[::193]/energy_scale)
    FminF0 = np.square(residue.reshape(41,193)[:,1:])
    #FminF0 = residue
    #print residue
    #print np.sum(FminF0)
    #print FminF0
    print "Energy residue error: ", np.sum(EminEo)/np.sum(E0**2)
    print "Energy RMSE error: ", np.sqrt(np.sum(EminEo))
    print "Force residue error: ", np.sum(FminF0)/np.sum(F0**2)
    print "Force RMSE error: ", np.sqrt(np.sum(FminF0))

    return x    
    #print x

def expandingx(unik,x,k):
    #print dic
    p = dict(zip(unik,x))
    #print p
    vx = np.zeros(len(k),dtype=complex)
    #idx = np.arange(len(k))
    for i in range(len(k)):
        magk = np.linalg.norm(k[i])
        if magk in unik:
            vx[i] = p[magk]
    #print vx
    return vx

def plotting_k(x):

    k = np.loadtxt("k_vectors.dat")
    #print len(k)
    uniqK = np.sort(np.unique(np.linalg.norm(k,axis=1)))
    vx = expandingx(uniqK,x,k)

    # plt.scatter(uniqK,x*1000)
    # plt.ylabel("Coefficient of k (e-03)")
    # plt.ylim([-1e-01,0.2])
    # plt.xlabel("k vectors")
    # plt.title("Plotting V_coeff for rho")
    # plt.savefig("Plotting_V2body.png")
    # plt.show()
    # plt.close()

    posits=create_r(8)

    Vlong = []
    
    for r in posits:
        # kdotr = dict()
        # for l in k:
        #     if np.linalg.norm(l) in kdotr:
        #         continue
        #     else:
        #         if np.linalg.norm(l) in uniqK:
        #             kdotr[np.linalg.norm(l)]= np.dot(l,r).tolist()
        #print kdotr.values()
        #print len(kr)
        #assert 1==0
        #print "length of kr", len(kr)
        #print "BLAAAHHH\n"
        #assert 1==0
        #print np.exp(1j*kr)
        #assert 1==0
        kr = np.asarray([np.dot(l,r) for l in k])
        #print len(kr)
        Vlong.append(np.sum(vx*np.exp(1j*kr)))

    Vlong = np.asarray(Vlong)
    #print len(Vlong)
    #assert 1==0
    #print Vlong
    #assert 1==0
    Vpfit = np.loadtxt("dummy.plot")[:1000]
    a,y = avg_r(posits,Vlong)
    #print y
    #assert 1==0
    #Vhat = savgol_filter(f, 31, 3)
    plt.plot(a,(7.4152**3)*np.absolute(y),c="r", label="2 body long range")
    plt.plot(Vpfit[:,0],Vpfit[:,1],label="Potfit short range")
    #plt.ylim([-5e-03,5e-03])
    plt.xlabel("r (bohr)")
    plt.ylabel("V (a.u)")
    plt.title("Adding 2 bdy term to potential")
    plt.legend()
    plt.savefig("2bodylrV2.png")
    plt.show()
    plt.close()
        # print abs(b.reshape(48,193)).mean()
    # print abs(FminF0).mean()
    # print np.sqrt((FminF0**2.).mean())
    #print "Chi of E: ", np.sqrt(np.mean(EminEo))/(energy_scale)
#     def calc_error(m):

#         resi=m.reshape(N,-1)
#         f_error=[]
#         e_error=[]

#         for i in range(N):
#             e_error.append(resi[i,0]/energy_scale)
# #            print resi[i,0]/energy_scale
#             force=resi[i,-192:]
#             ave=np.average(abs(force))
#             f_error.append(ave)

#         error0 = np.loadtxt('error.dat')
#         error0=error0[start:start+N]
#         e_xi_squared=0
#         f_xi_squared=0
#         for i in range(len(e_error)):
# #            print e_error[i],error0[i][0]
#             e_xi_squared+=e_error[i]*e_error[i]/(error0[i][0]*error0[i][0])
#             f_xi_squared+=f_error[i]*f_error[i]/(error0[i][1]*error0[i][1])*64*3
#         print '#energy xi^2:',e_xi_squared,'force xi^2:',f_xi_squared

#         print "#total xi^2:",e_xi_squared+f_xi_squared

#     calc_error(residue)
#     calc_error(b)

#run(0,16,True)
x = run(0,50,True)
plotting_k(x)
# run(16,15,True)
# run(31,8,True)

# run(0,16,False)
# run(16,15,False)
# run(31,8,False)
