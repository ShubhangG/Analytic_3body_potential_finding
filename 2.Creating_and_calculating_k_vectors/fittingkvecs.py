#!/usr/bin/env python
import numpy as np
import pandas as pd
import lmfit
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#from lmfit import minimize, Parameters

def unfold(kvecs,vals):
  from sympy.utilities.iterables import multiset_permutations

  k1 = kvecs.tolist()
  v1 = vals.tolist()

  for ik in range(len(kvecs)):
    kvec = kvecs[ik]
    val  = vals[ik]

    for perm in multiset_permutations(kvec):
      if np.allclose(perm,kvec):
        continue
      # end if
      k1.append(perm)
      v1.append(val)
    # end for
  # end for ik

  k2 = np.array(k1)
  v2 = np.array(v1)
  
  #return np.append(-k2,k2).reshape(-1,3),np.append(v2,v2)
  return k2,v2
# end def unfold

def reciprocal(k):
	return np.asarray([1.0/0.0001 if x==0 else 1.0/x for x in k])
#end def reciprocal

def residual(params,x,data):
	epsil=params['epsil']
	model = 1.0/epsil*x
	return (model-data)
#end def residual

def fit_lm(df):
	k1s = reciprocal(df['k1'].values)
	k2s = reciprocal(df['k2'].values)
	k3s = reciprocal(df['k3'].values)

	Vx = df['x'].values
	eps = np.zeros(len(Vx))
	eps = np.asarray([l+0.00042 for l in eps])
	model= lmfit.models.ExpressionModel("epsil1 * epsil2 * epsil3 * x")
	params = model.make_params(epsil1=0.00014,epsil2=0.00014,epsil3=0.00014)
	#params = model.guess(Vx,x=k1s*k2s*k3s)
	#params.add('epsil',value=0.00042)
	fit = model.fit(Vx,params,x=k1s*k2s*k3s)
	print fit.fit_report()
	#fit.plot_fit()
#end def fit_lm
def fit_method2(k,v):
	rK = k[:,0]*k[:,1]*k[:,2]
	epsK=v*rK
	#print Eps_in
	model = lmfit.models.PolynomialModel(2, independent_vars=['x'], missing='drop')
	#params = model.guess(epsK,x=rK)
	fit = model.fit(epsK,x=rK)
	print fit.fit_report()
#end def fit_method2

def fit_method3(x,*e):
	k1 = x/17**2
	k2 = (x-k1*17*17)/17
	k3 = x-(x/17)*17
	#print k1,k2,k3
	e=np.array(e)
	return e[17]*1.0/((k1)*(k2)*(k3))*(1.0/e[k1]-1)*(1.0/e[k2]-1)*(1.0/e[k3]-1)

def fit_method4(x,*e):
	k1 = x/17**2
	k2 = (x-k1*17*17)/17
	k3 = x-(x/17)*17
	#print k1,k2,k3
	e=np.array(e[0])
	return e[17]*1.0/((k1)*(k2)*(k3))*(1.0/e[k1]-1)*(1.0/e[k2]-1)*(1.0/e[k3]-1)

fl = 0
def fit_method5(X,*e):
	k1,k2,k3 = X
	if fl==1:
		e = np.array(e[0])
	return e[0]*np.cos((k1-e[1])*e[2])*np.exp(-k1/e[3])*np.cos((k2-e[4])*e[5])*np.exp(-k2/e[6])*np.cos((k1-e[7])*e[8])*np.exp(-k3/e[9])

def fit_method6(X,*e):
	k1,k2,k3 = X
	if fl==1:
		e= np.array(e[0])
	return e[0]*np.sin(k1*e[1])/(k1*e[1])*np.sin(k2*e[2])/(k2*e[2])*np.sin(k3*e[3])/(k3*e[3])

def create_dat4fitting(k1,k2,k3):
	e1 = np.linspace(0.1,0.4,len(k1))
	e2 = np.linspace(0.5,0.8,len(k2))
	e3= np.linspace(1.0,1.4,len(k3))
	e4 = 0.5
	return e4*1.0/(k1*k2*k3)*(1.0/e1-1)*(1.0/e2-1)*(1.0/e3-1)

def testfit(fX):
	fV = create_dat4fitting
	fe0s = np.array([0.1 for i in range(18)])
	fe0s[17] = 1.0
	print fV
	def fake_fitmeth(x,*e):
		k1 = x/5**2
		k2 = (x-k1*5*5)/5
		k3 = x-(x/5)*5
		#print k1,k2,k3
		e=np.array(e)
		return e[5]*1.0/((k1)*(k2)*(k3))*(1.0/e[k1]-1)*(1.0/e[k2]-1)*(1.0/e[k3]-1)
	
	popt, pcov = curve_fit(fake_fitmeth,fX,fV,p0=fe0s)
	print popt
	

if __name__ == '__main__':
	df = pd.read_json('Mykvalues.json')
	sel = (df[['k1','k2','k3']].apply(lambda x: x==0)).any(axis=1)
	#print sel
	df1 = df[~sel]
	#print df1
	#df1=df
	pos  = df1.loc[:,['k1','k2','k3']].values
	val_arr  = df1.loc[:,['x']].values.flatten()
	# print max(val_arr)
	# assert 0
 #  	k1,v1 = unfold(pos,val_arr)
 	k1 = pos[:,0]
 	k2= pos[:,1]
 	k3 = pos[:,2]

 	x = pos[:,0]*17**2 + pos[:,1]*17 + pos[:,2]
 	#print pos[5]
 # 	rK = pos[:,0]*pos[:,1]*pos[:,2]
	# epsK=val_arr*rK

	e = np.linspace(0.1,1.0,18)
	fV = e[17]*1.0/((k1)*(k2)*(k3))*(1.0/e[k1]-1)*(1.0/e[k2]-1)*(1.0/e[k3]-1)
 	#print fV
 	#assert 1==0
 	e0s= np.zeros(17)
 	e0s=np.array([0 for i in range(18)])
 	e0s[17] = 1.0

 	##e0s[[3,7,12]] = -0.00042
 	#print fit_method3(x,e0s)
	 
 	X = (pos[:,0],pos[:,1],pos[:,2])
 	#e0s = np.asarray([1e-05,1,np.pi/3,0.5,0.1,np.pi/3,0.5,0.1,np.pi/3,0.5])
 	#print e0s
 	#e0s = np.asarray([0.1,0.1,0.1,0.1])
 	popt, pcov = curve_fit(fit_method3,x,fV,p0=e0s)
 	#print x.dtype
 	print popt
 	fl=1
 	plt.plot(fV,fit_method4(x,popt),'o')
 	#plt.plot([0,6e-06],[0,6e-06])
 	plt.show()
 	#print po
  	#print pd.DataFrame(k1,v1)
  	#fit_method2(k1,v1)


	# plt.show()
	# x = k1s*k2s*k3s
	# print x
	# #print Vx**2
	# #print residual(params,x,Vx)
	# out= minimize(residual,params,method='leastsq',args=(x,),kws={'data':Vx})

	# print out.params['epsil