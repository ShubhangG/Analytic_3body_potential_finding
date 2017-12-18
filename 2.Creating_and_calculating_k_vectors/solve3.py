#!/usr/bin/env python
#########################################
##Code authors
#Chong Han, chan104@illinois.edu 
#Shubhang Goswami, sgoswam3@illinois.edu
#########################################

import math
import numpy as np
import pandas as pd
import itertools
from scipy import linalg
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import colors
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import griddata
import matplotlib as mpl 

energy_scale=10

n_parameters=154

def run(start,N,flag):
	print 'start:',start,'N:',N,'flag:',flag
	if flag:
		select_forces=True
		select_energy=True
		neglect_constant=False
	else:
		select_forces=True
		select_energy=False
		neglect_constant=True
	#since there are 15 parameters to fit
	a=np.loadtxt('third_degree.dat')
	dat=np.loadtxt('data/SG_dft7Krs1.15.dat')
	E0 = np.loadtxt('data/SG7k_r1.15_eam1.energy',usecols=[4])
	F0 = np.loadtxt('data/SG7k_r1.15_eam1.force',usecols=[5])
	#dat=np.loadtxt('data/qmc.dat')
	#print len(a)
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

	x, resid, rank, sigma = linalg.lstsq(a, b)
	xx=x.reshape(1,-1)[0]
	k=np.loadtxt('third_degree_k.dat')
	
	#print '#residue',resid
	s=0
	for i in range(len(b)):
		s+=b[i]*b[i]
	print '#original:',s

	residue= b-np.dot(a,x)
	#print b.shape,a.shape,x.shape
	#print residue.shape
	EminEo = np.square(residue[::193]/energy_scale)
	FminF0 = np.square(residue.reshape(41,193)[:,1:])
	print "Energy residue error: ", np.sum(EminEo)/np.sum(E0**2)
	print "Energy RMSE error: ", np.sqrt(np.sum(EminEo))
	print "Force residue error: ", np.sum(FminF0)/np.sum(F0**2)
	print "Force RMSE error: ", np.sqrt(np.sum(FminF0))
	#dfdat = pd.read_json("Mykvalues.json")
	col_dim = x[:,0]
	print "max V value :", max(col_dim)
	print "min V value :", min(col_dim)
	#print vecs
	tdegk = np.loadtxt("third_degree_k.dat")
	contourplot(tdegk,col_dim)
	datdf = pd.DataFrame(np.column_stack((tdegk,col_dim)),columns=["k1","k2","k3","x"])
	datdf.to_json("purekvalues.json")
	
	#plotlines2(datdf)
	#binplot(datdf)
	#print tdegk[:,0]
	#thr3Dplot(tdegk[:,0],tdegk[:,1],tdegk[:,2],col_dim)

def thr3Dplot(k1,k2,k3,col_dim):
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	k1,k2=np.meshgrid(k1,k2)
	k3=np.outer(k3.T,k3)
	col_dim = np.outer(col_dim.T,col_dim)
	minn,maxx =col_dim.min(), col_dim.max()
	normi = colors.Normalize(minn,maxx)
	m=plt.cm.ScalarMappable(norm=normi,cmap='jet')
	m.set_array([])
	fcolors=m.to_rgba(col_dim)
	print fcolors.shape
	#print fcolors
	print k1.shape, k2.shape, k3.shape
	#assert 0
	surf = ax.plot_surface(k1,k2,k3,cmap=cm.coolwarm,rstride=1,cstride=1,facecolors=fcolors,vmin=minn,vmax=maxx,shade=False)
	ax.set_xlabel('|k1|')
	ax.set_ylabel('|k2|')
	ax.set_zlabel('|k3|')
	ax.zaxis.set_major_locator(LinearLocator(10))
	cbar= fig.colorbar(m,shrink=0.5,aspect=5)
	cbar.set_label("The coefficient of threebody term")
	plt.savefig("3dksandv.png")
	fig.canvas.show()
	plt.show()

def plotlines2(datdf):
	marklist = ["o","s","v","^","<",">","H","8","p","P","*","h","D","+","x"]
	col = np.linspace(0,1,len(np.unique(datdf["k1"]))+1)
	mrkctr=0
	for term in np.unique(datdf["k1"]):
		sel = (datdf["k1"]==term)
		k3_norm = normalize(datdf[sel]["k3"],'max')
		plt.scatter(datdf[sel]["k2"].values,datdf[sel]["x"].values,c=k3_norm,s=80,marker=marklist[int(term)],cmap=cm.coolwarm,label="|k1|= "+repr(term))
		mrkctr+=1
	cbar = plt.colorbar()
	cbar.set_label("|k3| value normalized",fontsize=22)
	#mpl.rcParams.update({'font.size':22})
	plt.xlabel("|k2|",fontsize=22)
	plt.ylabel("Potential",fontsize=22)
	plt.title("Scatter color plot to visualize potential",fontsize=22)
	plt.legend(loc='upper center',bbox_to_anchor=(0.5,-0.05),ncol=8,fancybox=True,shadow=True)
	plt.ylim([-0.0002,0.0002])
	plt.savefig("scatter3body.png")
	plt.show()

def contourplot(tdegk, col_dim):
	k1 = tdegk[:,0]
	k2 = tdegk[:,1]
	k3 = tdegk[:,2]
	X,Y = np.meshgrid(k1,k2)
	Z = col_dim
	xi = np.linspace(min(k1),max(k1))
	yi = np.linspace(min(k2),max(k2))
	#griddata(k1, k2, grid_x, grid_y, method='cubic')
	zi = griddata((k1,k2),Z,(xi[None,:],yi[:,None]), method='cubic')
	zj = griddata((k2,k1),Z,(xi[None,:],yi[:,None]),method='cubic')
	#levels = np.linspace(-1,1,40)
	cs = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
	cs = plt.contourf(xi,yi,zi,15,cmap=cm.jet)# levels=levels)
	# cs1 = plt.contour(yi,xi,zj,15,linewidths=0.5,colors='k')
	# cs1 = plt.contourf(yi,xi,zj,15,cmap=cm.jet)
	cbar = plt.colorbar()
	cbar.set_label("The coefficient of threebody term")
	plt.grid(linewidth=0.2)
	#plt.scatter(k1,k2,marker='o',c='b',s=5)
	plt.xlabel("|k1|")
	plt.ylabel("|k2|")
	plt.title("Plotting three body on paper")
	plt.savefig("contour_plot3body.png")
	plt.show()

def binplot(datadf):
	col = cm.rainbow(np.linspace(0,1,len(np.unique(datadf["k1"]))))
	ctr=0
	for term,c in zip(np.unique(datadf["k1"]),col):
		sel = (datadf["k1"]==term) & (datadf["k2"]>4)
		if ctr<5:
			ctr+=1
			continue

		x = datadf[sel]["k3"]
		y = datadf[sel]["x"]
		if x.size >1:
			lists = sorted(itertools.izip(*[x,y]))
			new_x, new_y = list(itertools.izip(*lists))
		else:
			new_x = x
			new_y = y
		plt.plot(new_x,new_y,"-o",c=c,label="|k1|^2="+repr(term)+" , |k2|^2>"+str(4))
		ctr+=1

	plt.xlabel("|k3|^2")
	plt.ylabel("Vpot")
	#plt.ylim([-0.00006,0.00006])
	plt.title("Potential vs 3rd k term")
	plt.legend()
	plt.savefig("binningplotk3vsV_5_10k2.png")
	plt.show()


run(0,41,True)
# run(16,15,True)
# run(31,8,True)

#########################
    #print residue
    #print np.sum(FminF0)
    #print FminF0
    #print "Energy residue error: ", np.sum(EminEo)/np.sum(E0)

    	# def calc_error(m):
	# 	resi=m.reshape(N,-1)
	# 	f_error=[]
	# 	e_error=[]

	# 	for i in range(N):
	# 		e_error.append(resi[i,0]/energy_scale)
	# 		force=resi[i,-192:]
	# 		ave=np.average(abs(force))
	# 		f_error.append(ave)

	# 	error0 = np.loadtxt('error.dat')
	# 	error0=error0[start:start+

	# 	e_xi_squared=0
	# 	f_xi_squared=0
	# 	for i in range(len(e_error)):
	# 	    e_xi_squared+=e_error[i]*e_error[i]/(error0[i][0]*error0[i][0])
	# 	    f_xi_squared+=f_error[i]*f_error[i]/(error0[i][1]*error0[i][1])*64*3
	# 	#print '#energy xi^2:',e_xi_squared,'force xi^2:',f_xi_squared
	# 	#print "#total xi^2:",e_xi_squared+f_xi_squared
	# 	print e_xi_squared,f_xi_squared,e_xi_squared+f_xi_squared
		

	# calc_error(b)
	# calc_error(residue)

	# seq=np.argsort(-1*abs(xx))
	# print seq
	# '''
	# for i in seq:
	# 	print k[i]
	# 	print xx[i]
	# '''
	# newa=a[:,seq]

	# for i in range(len(newa[0])):
	# 	x, resid, rank, sigma = linalg.lstsq(newa[:,:i+1], b)
	# 	residue= b-np.dot(newa[:,:i+1],x)
	# 	calc_error(residue)
####3-D Plotting####
	# fig = plt.figure()
	# ax = fig.gca(projection='3d')
	# k1,k2=np.meshgrid(k1,k2)
	# k3=np.outer(k3.T,k3)
	# minn,maxx =col_dim.min(), col_dim.max()
	# normi = colors.Normalize(minn,maxx)
	# m=plt.cm.ScalarMappable(norm=normi,cmap='jet')
	# m.set_array([])
	# fcolors=m.to_rgba(col_dim)
	# #print len(fcolors)
	# #print fcolors
	# #print len(k1), len(k2), len(k3)
	# surf = ax.plot_surface(k1,k2,k3,rstride=1,cstride=1,facecolors=fcolors,vmin=minn,vmax=maxx,shade=False)
	# ax.set_xlabel('|k1|')
	# ax.set_ylabel('|k2|')
	# ax.set_zlabel('|k3|')
	# #ax.zaxis.set_major_locator(LinearLocator(10))
	# fig.colorbar(surf,shrink=0.5,aspect=5)
	# fig.colorbar.set_label("The coefficient of threebody term")
	# plt.savefig("3dksandv.png")
	# # fig.canvas.show()
	# fig = plt.figure()
	# ax = fig.gca(projection='3d')
	# k1, k2, k3, colors = tdegk[:,0], tdegk[:,1], tdegk[:,2], col_dim
	# grid_x, grid_y = np.mgrid[min(k1):max(k1):50j, min(k2):max(k2):50j]
	# z = griddata((k1, k2), k3, (grid_x, grid_y), method='cubic')
	# colors = griddata((k1, k2), colors, (grid_x, grid_y), method='cubic')
	# k2, k1 = grid_y.ravel(), grid_x.ravel()

	# y_coords, x_coords = np.unique(k2), np.unique(k1)
	# y_centers, x_centers = [ arr[:-1] + np.diff(arr)/2 for arr in (y_coords, x_coords)]

	# # Convert back to a 2D grid, required for plot_surface:
	# Y = k2.reshape(y_coords.size, -1)
	# X = k1.reshape(-1, x_coords.size)
	# Z = z.reshape(X.shape)
	# C = colors.reshape(X.shape)

	# #Normalize the colors to fit in the range 0-1, ready for using in the colormap:
	# Z -= Z.min()
	# Z /= Z.max()

	# C -= C.min()
	# C /= C.max()

	# interp_func = RectBivariateSpline(x_coords, y_coords, C.T, kx=1, ky=1) # the kx, ky define the order of interpolation. Keep it simple, use linear interpolation.
	# ax.set_xlabel('|k1|')
	# ax.set_ylabel('|k2|')
	# ax.set_zlabel('|k3|')

	# surf =ax.plot_surface(X,Y,Z, facecolors=cm.hot(interp_func(x_centers, y_centers).T), rstride=1,  cstride=1, linewidth=1, antialiased=False, shade=False) # only added because of this very limited dataset
	# #m = cm.ScalarMappable(cmap=cm.jet)
	# #m.set_array(C)
	# #cbar = plt.colorbar(m,shrink=0.5,aspect=5)
	# #fig.colorbar(surf,shrink=0.5,aspect=5)
	# #cbar.set_label("The coefficient of threebody term")
	# plt.savefig("3dksandv.png")
	# plt.show()