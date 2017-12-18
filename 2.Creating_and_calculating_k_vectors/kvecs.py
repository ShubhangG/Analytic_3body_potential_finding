#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors

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

if __name__ == '__main__':

  df = pd.read_json('Mykvalues.json')
  axes = 36.*(np.eye(3)-0.5*np.eye(3))
  pos  = df.loc[:,['k1','k2','k3']].values
  val_arr  = df.loc[:,['x']].values.flatten()

  k1,v1 = unfold(pos,val_arr)
  pos = k1
  val_arr = v1

  vmin,vmax= val_arr.min(), val_arr.max()

  fig = plt.figure()
  from mpl_toolkits.mplot3d import Axes3D
  ax = fig.add_subplot(1,1,1,projection="3d")
  ax.set_xlabel("k1")
  ax.set_ylabel("k2")
  ax.set_zlabel("k3")
  from qharv.inspect import crystal
  crystal.draw_cell(ax,axes)

  uValpts = []
  cmap = plt.get_cmap('rainbow')
  for ipos in range(len(pos)):
    val = val_arr[ipos]
    if abs(val) < 2e-5:
      continue
    uval= (val-vmin)/(vmax-vmin)
    rgb = cmap(uval)
    cs = crystal.draw_atoms(ax,np.array([pos[ipos]]),c=rgb)#,ms=5,alpha=0.5)
    uValpts.append(uval)  
  # end for
  fig1 = plt.figure()
  ax1 = fig1.gca()
  s=ax1.scatter(0,0,s=0,facecolors='none')
  cptnorm = colors.Normalize(vmin=min(uValpts),vmax=max(uValpts))
  #scalmap = cm.ScalarMappable(norm=cptnorm,cmap=cmap)
  #cbar = plt.colorbar(s,norm=cptnorm)
  #cbar.set_label("Vk value",fontsize=22)
  #plt.savefig("3d_plotKvsV.png")
  plt.show()

# end __main__
