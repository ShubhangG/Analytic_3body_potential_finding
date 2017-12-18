#!/usr/bin/env python
#########################################
##Code authors
#Shubhang Goswami, sgoswam3@illinois.edu 
#########################################
import pandas as pd 
import numpy as np
import pickle as pk 

if __name__ == '__main__':
	#Read in parsed forces and energy into a pandas df
	allf = pd.read_table("forces.dat", delim_whitespace=True, header=None, skiprows=65, skipfooter=64, engine='python')
	energy = pd.read_table("energy.dat", delim_whitespace=True, header=None,usecols=[4])
	pos = pd.read_table("positions.dat", delim_whitespace=True, header=None, skipfooter=64, engine='python')
	pos.replace(to_replace='H', value=0, inplace=True)	#take only the numbers
	E = energy[4].values/(64*2)							#Convert Rydberg to Hartree/Atom
	split_len=len(allf)/64								#split into number of iterations

	alat = 7.4152										#box length for simulation
	ctr = 0
	n_force = allf.iloc[:,6:9]							
	n_force = n_force.divide(2)							#converting Ry to H 
	net = pd.concat([pos,n_force], axis=1)
	del net[0]
	#convert position range betwenn 0 to L instead of -L/2 to L/2
	net[1] = (net[1]/alat - np.floor(net[1]/alat))*alat 
	net[2] = (net[2]/alat - np.floor(net[2]/alat))*alat
	net[3] = (net[3]/alat - np.floor(net[3]/alat))*alat

	N_eff= pk.load(open("N_eff.pk", "rb"))				#Effective number of uncorrelated points
	N_eff = np.floor(N_eff)
	#print N_eff

	outf = open('Potfit/potfitforce7000K.dat', 'w')			
	pre = "#N 64 1\n#X    7.4152   0.00000000    0.00000000\n#Y    0.00000000    7.4152   0.00000000\n#Z    0.00000000    0.00000000    7.4152\n"

	ctr = 0
	for tstep in np.array_split(net,split_len):
		if(ctr%N_eff !=0):								#Skips correlated entries
			ctr+=1
			continue

		#Writes in the format taken by Potfit
		outf.write(pre + '#E  '+ repr(E[ctr]) + '\n' + '#F\n')
		for item in tstep.values:
			outf.write('0 ')
			outf.write(' '.join(map(str,item)))
			outf.write('\n')
		ctr+=1
		
	outf.close()

