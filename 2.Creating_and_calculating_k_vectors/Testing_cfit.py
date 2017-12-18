#!/usr/bin/env python
import numpy as np
import pandas as pd
import lmfit
from scipy.optimize import curve_fit
from scipy.optimize import basinhopping
import matplotlib.pyplot as plt
import scipy.interpolate
from sklearn.metrics import mean_squared_error

def fit_method0(X,*e):
  k1,k2,k3 = X
  #print k1,k2,k3
  e=np.array(e)
  return e[17]*1.0/((k1)*(k2)*(k3))*(1.0/e[k1]-1)*(1.0/e[k2]-1)*(1.0/e[k3]-1) + 0*np.linalg.norm(e)

def fit_method1(X,*e):
  k1,k2,k3 = X
  e=np.array(e)
  return e[k1]*e[k2]*e[k3] + 0*np.linalg.norm(e)

def fit_method2(X,*e):
  k1,k2,k3 = X
  return e[0]*np.cos((k1-e[1])*e[2])*np.exp(-k1/e[3])*np.cos((k2-e[4])*e[5])*np.exp(-k2/e[6])*np.cos((k3-e[7])*e[8])*np.exp(-k3/e[9]) + 0*np.linalg.norm(e)

def fit_method3(X,*e):
  k1,k2,k3 = X

  return e[0]*np.sin(k1*e[1])/(k1*e[1])*np.sin(k2*e[2])/(k2*e[2])*np.sin(k3*e[3])/(k3*e[3])

def fit_method4(X,*e):
  k1,k2,k3 = X
  e=np.array(e)
  return e[17]*np.cos((k1)*e[k1])*np.exp(-k1/e[k1])*np.cos((k2)*e[k2])*np.exp(-k2/e[k2])*np.cos((k3)*e[k3])*np.exp(-k3/e[k3]) + 0*np.linalg.norm(e)

def fit_method5(X,*e):
  k1,k2,k3 = X
  return e[0]*np.cos((k1-e[1])*e[2])*np.exp(-k1/e[3])*np.cos((k2-e[4])*e[5])*np.exp(-k2/e[6])*np.cos((k3-e[7])*e[8])*np.exp(-k3/e[9]) + e[10]/((k1+0.0001)*(k2+0.0001)*(k3+0.0001))
#def diff_nfit(k1,k2,k3,Vk):

def fit_method6(X,*e):
  k1,k2,k3 = X
  e=np.array(e)
  return e[17]*np.cos((k1-e[18])*e[19])*np.cos((k2-e[20])*e[21])*np.cos((k3-e[22])*e[23]) + e[k1]*e[k2]*e[k3]/((k1+0.0001)*(k2+0.0001)*(k3+0.0001))

def fit_method7(X,*e):
  k1,k2,k3 = X
  e=np.array(e)
  return e[17]*np.cos((k1-e[18])*e[19])*np.exp(-k1/e[24])*np.cos((k2-e[20])*e[21])*np.exp(-k2/e[25])*np.cos((k3-e[22])*e[23])*np.exp(-k3/e[26]) + e[k1]*e[k2]*e[k3]/((k1+0.0001)*(k2+0.0001)*(k3+0.0001))

def genparamfrom2(params):
  e = np.arange(17)
  nprm = params[0]*np.cos((e-params[1])*params[2])*np.exp(-e/params[3])
  return nprm;

def plot_stuff(X,val_arr,popt1,model1,popt2,model2):
   fig,ax = plt.subplots(1,1)

   prediction = model1(X,*popt1) #+ model2(X,*popt2)
   print np.sqrt(mean_squared_error(val_arr,prediction))
   #print prediction
   plt.plot(val_arr,prediction,'o')

   #plt.xlim([-0.00004,0.00002])
   #plt.ylim([-0.00004,0.00002])
   plt.xlabel("Real value")
   plt.ylabel("Predicted value")
   plt.plot(val_arr,val_arr)
   plt.show()

def simanneal(model,guess,niter,step,t):
    finalPams =  basinhopping(model,guess,niter=niter,stepsize=step,T=t)
    f=open('annealing_step'+repr(step)+"_iter"+repr(niter)+"_T"+repr(t)+".dat","w")
    f.write(repr(finalPams))
    f.close()
    return finalPams.x, finalPams.fun




if __name__ == '__main__':
   df = pd.read_json('purekvalues.json')
   sel = (df[['k1','k2','k3']].apply(lambda x: x==0)).any(axis=1)
   selmin = (np.isclose(df['x'],-5.72637e-05,rtol=1e-06))
   df1 = df#[~selmin]
   pos  = df1.loc[:,['k1','k2','k3']].values

   val_arr  = df1.loc[:,['x']].values.flatten()
   k1 = pos[:,0]
   k2= pos[:,1]
   k3 = pos[:,2]

   val_arr_new = val_arr*1e05
   #x = pos[:,0]*17**2 + pos[:,1]*17 + pos[:,2]
   #e = np.random.rand(18)
   #fV = e[17]*1.0/((k1)*(k2)*(k3))*(1.0/e[k1]-1)*(1.0/e[k2]-1)*(1.0/e[k3]-1)


   mlist  = [fit_method0,fit_method1,fit_method2,fit_method3, fit_method4, fit_method5, fit_method6, fit_method7]
   nparams= [18,17,10,10,18,11,24,27]

   imodel = 2
   model  = mlist[imodel]
   nparam = nparams[imodel]

   #print val_arr.min(),val_arr.max()
   e0s= np.zeros(nparam)
   e0s=np.array([17 for i in range(nparam)]) #1.7607 for method2 and 1.6 if you remove smallest
   #e0s[17:nparam]= 1.7607
   #e0s[nparam-1] = 2e-06
   #e0s[1:11] = [ -1.97064512e+00,-2.81365470e+00,5.56149770e-01,-1.22714787e+03,1.50077927e+00,3.14573643e+00,2.92169920e+00,1.69773224e+00,9.12906113e-01,1.03941233e+00]
   #e0s[17] = 1.0
   def func4anneal1(e):
     return np.sqrt(mean_squared_error(val_arr_new,e[k1]*e[k2]*e[k3]))
  

   def cnvrttoks(optParam):
     return np.asarray([optParam[0]*np.cos((k-optParam[1])*optParam[2])*np.exp(-k/optParam[3]) for k in range(0,17)])


   X = (pos[:,0],pos[:,1],pos[:,2])

   popt, pcov = curve_fit(model,X,val_arr_new,p0=e0s)
   ne0s = genparamfrom2(popt)
   #print ne0s
   npopt, npcov = curve_fit(fit_method1,X,val_arr_new,p0=ne0s)
   #print npopt

   optPm = cnvrttoks(popt)
   stepsze = 0.5
   niter=200

   pams, rmse =  simanneal(func4anneal1,optPm,niter,stepsze,10)
   print "RMSE", rmse
   pams = np.asarray(pams)
   pams, rmse = simanneal(func4anneal1,optPm,niter,stepsze,1)
   print "RMSE", rmse

   predic = pams[k1]*pams[k2]*pams[k3]

   plt.scatter(val_arr_new,predic)
   plt.xlabel("True Value")
   plt.ylabel("Predicted Value")
   plt.plot(val_arr_new,val_arr_new,c="r")
   plt.savefig("Annealed_plot.png")
   plt.show()
   #plot_stuff(X,val_arr, npopt,fit_method1, popt, model)
  # 