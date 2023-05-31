import numpy as np
from numpy import polyfit
import matplotlib.pyplot as plt
import re
import time
import glob
import sys
import scipy.constants as scicon
import math


is_sorted = lambda a: np.all(a[:-1] <= a[1:])

pathx = "./eqlbm3.300.*0000.custom"

filex = glob.glob(pathx)

size = np.loadtxt(filex[0], delimiter=' ', comments='#', skiprows=5, max_rows=3)

print('size=\n',size)

coordinates = np.loadtxt(filex[0], delimiter=' ', comments='#', skiprows=9, max_rows=None)

#coordinates = np.delete(coordinates, 1, 1)

print('shape of the mapped array =', np.shape(coordinates))
#print(is_sorted(coordinates[:, 0]))

coordinates = coordinates[coordinates[:, 0].argsort()]

#print(is_sorted(coordinates[:, 0]))

Ntotal=len(coordinates)
print('Total number of atoms=',Ntotal)
local_ratio=0.3				## phonon localization limit
propagate_ratio=0.4
print('Maximum participation ratio of localized phonons=',local_ratio,'\nMinimum participation ratio of propagating phonons=',propagate_ratio)

#################################################################
######################################### Reading the eigenvector file
pathvec= "./eigvec.dat"
filevec = glob.glob(pathvec)
sysdim = 3
Nphonons=sysdim*Ntotal

modes = np.array([])
parti = np.array([])
spatialarr = np.zeros((Ntotal,1))
npropagating_modes = 0

startread = time.time()
print('started reading at time (s)=',startread)

skip = 1
with open(pathvec,"r",encoding='utf-8') as f:
	for i in range(Nphonons):
		for j in range(skip):
			next(f)
		line = f.readline()
#		print(line)
		c0 = re.split(r'\s+|\s|, |,',line)
		c = [ele for ele in c0 if ele.strip()]
		c1 = np.array(c[:])	
		freq = c1[7].astype(np.float32)
		modes=np.append(modes,freq)
		if i%500 == 0:
			print(i,' mode frequency=',modes[-1],' THz')
		next(f)
		use = 8
		maxr = Ntotal	
		atom = 0
		eigvectemp=np.array([])
		for line in f:
#			print('inner loop')
			atom += 1
			c0 = re.split(r'\s+|\s|, |,',line)
			c = [ele for ele in c0 if ele.strip()]
			c1 = np.array(c[:])
			eigvectemp = np.append(eigvectemp, c1[8].astype(np.float32))
			if atom == Ntotal:
				break
#		print('inner loop broke')
		exponent = (scicon.hbar*2*scicon.pi*freq*10**12)/(scicon.k*coordinates[:,5])
	#	print(coordinates[:,5])
	#	print('shape of exponent=', np.shape(exponent))
		noccu = 1/(np.power(math.e,exponent)-1)
	#	print(noccu)
	#	print(math.e,scicon.k,scicon.hbar,scicon.pi,)
		spatialtemp =  (noccu+0.5)*scicon.hbar*2*scicon.pi*freq*10**12*eigvectemp**2
		parti_deno = np.sum(eigvectemp**4)
		partitemp = 1/(Ntotal*parti_deno)
		parti = np.append(parti,partitemp)
	#	eigvecarr = np.add(eigvecarr,spatialtemp.reshape(Ntotal,1))
		if partitemp >= propagate_ratio:
			spatialarr = np.add(spatialarr,spatialtemp.reshape(Ntotal,1))
			npropagating_modes += 1
		
	#	print(tempspatial)
	#	skip += Ntotal+1
	#	print(np.shape(eigvecarr))
	#	print(eigvecarr)
	#	print(np.mean(eigvecarr))
	#	print(skip)

endread = time.time()
print('time endread=',endread)
	
#eigvecarr = np.array(eigvecarr).reshape()
print(np.shape(spatialarr),np.max(spatialarr),np.min(spatialarr),np.mean(spatialarr))
print('number of propagating modes = ',npropagating_modes)
print(len(parti),np.max(parti),np.min(parti),np.mean(parti))

spatialarr = spatialarr*6.2415e18		## convert spatial energy distribution to eV

nlocal_modes = Nphonons - npropagating_modes

spatial = np.append(coordinates, spatialarr, axis=1)
spatialfname = 'spatial_data'+' '+str(propagate_ratio)+' '+str(npropagating_modes)+'.custom'
np.savetxt(spatialfname, spatial, delimiter=' ')

modes = np.append(modes.reshape(-1,1), parti.reshape(-1,1), axis=1)
partifname = 'participation_ratios'+' '+str(local_ratio)+' '+str(nlocal_modes)+'.data'
np.savetxt(partifname, modes, delimiter=' ')









