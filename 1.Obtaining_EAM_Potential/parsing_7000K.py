#!/usr/bin/env python
#########################################
##Code authors
#Shubhang Goswami, sgoswam3@illinois.edu 
#########################################
if __name__ == '__main__':
	tempf = "scf_7000K.out" 			#Name of the ab-initio scf calculation
	inf = open(tempf, 'r')
	pos = open("positions.dat", 'w') 	#Name of where we will store positions
	force = open("forces.dat",'w')		#Similarly forces and energy
	ener = open("energy.dat",'w')
	ctr= 64								#Number of atoms in our box
	#Loop through each line looking for energy/positions and forces
	for line in inf:					
		if ctr<64:
			pos.write(line)
			ctr+=1
		if line.lstrip().startswith('atom'):	#Force data starts with 'atom'
			force.write(line)
		if line.startswith("ATOMIC_POSITIONS"):
			ctr=0
		if line.startswith("!"):		#Energy starts with ! in scf output file
			ener.write(line)
	#Close all files
	inf.close()
	pos.close()
	force.close()
	ener.close()
