import numpy as np
from numpy import polyfit
import matplotlib.pyplot as plt
import re
import time

'''
This uses a str.data file which contains just the section on atomic coordinates of a typical lammps structure data file along with the eigvec.data output of phana package and log file from running lammps phonon module 
'''

is_sorted = lambda a: np.all(a[:-1] <= a[1:])


name="eqlbm.log" 		# Change according to simulation output file

xlo=0.029951160202565753	# Change according to model
xhi=5.01404903979739
ylo=0.15562346396823343
yhi=26.053755736031867
zlo=0.0
zhi=100

Ntotal=48			# Change according to model (if not read from the file automatically)

Nlocal=5			
localRatio=Nlocal/Ntotal
print(localRatio)

nulo=0
nuhi=18
ratiolo=0
ratiohi=0.8

xparts=(xhi-xlo)*2//3
xslab=(xhi-xlo)/xparts
yparts=(yhi-ylo)*2//3
yslab=(yhi-ylo)/yparts
zparts=(zhi-zlo)*2//3
zslab=(zhi-zlo)/zparts

#freqlo=7.597
#freqhi=7.598


modes = np.array([])
#eigvec = np.array([])
datapos = np.array([])
parti = np.array([])	
parti2 = np.array([])


######################## Reading Coordinates and Atomic ID

with open("str.data","r",encoding='utf-8') as fb:
	coordinates=np.array([])
	j=0
	for line in fb:
		if j==Ntotal:
			break
#		print(j)
		c0 = re.split(r'\s+|\s|: |, |:|,',line)
		c = [ele for ele in c0 if ele.strip()]
		c1 = np.array(c[:])
		coordinates = np.append(coordinates,np.array([c1[0].astype(int), c1[2].astype(np.float64), c1[3].astype(np.float64), c1[4].astype(np.float64)]))
		j += 1
		
coordinates=coordinates.reshape(-1,4)	

print(coordinates)
print('shape of coordinates array=',np.shape(coordinates))

coordinates=coordinates[coordinates[:, 0].argsort()]
print('New coordinates array=\n',coordinates)

######################## Reading Basis and Atomic ID mapping from phonon log file

with open(name,"r",encoding='utf-8') as fc:
	basis=np.array([])
	j=0
	for i in range(15):
#		print(i)
		next(fc)
	flag3=True
	for line in fc:
		if flag3:
			print('first line of basis=',line)
			flag3=False
		if j==Ntotal:
			break
#		print(j)
		c0 = re.split(r'\s+|\s|: |, |:|,',line)
		c = [ele for ele in c0 if ele.strip()]
		c1 = np.array(c[:])
#		print(c1)
		basis = np.append(basis,np.array([c1[3].astype(int), c1[4].astype(int)]))
		j += 1

basis=basis.reshape(-1,2)
print('shape of basis array=',np.shape(basis))

basis=basis[basis[:, 1].argsort()]
print('New basis array=\n',basis)

######################################################################## Map arrays

idsorted = np.append(coordinates,basis,axis=1)
idsorted = idsorted[:,:-1]
print('shape of combined mapping array=',np.shape(idsorted))
print('combined mapping array=',idsorted)

basissorted = idsorted[idsorted[:, 4].argsort()]
print('Basis sorted combined mapping array=',basissorted)


#####################################################################################
	
"""
using map.in and data.pos files in lammps Fixphonon simulation ensures that the atom basis number k and atomic IDs share the same sequence and there is no ambiguity - hence directly the atomic ID can be used to analyze the phonon output without confusion regarding whether to use the k sequence or the atomic ID sequence for it.
"""
startread = time.time()
print('started reading at time (s)=',startread)
with open("eigvec.dat","r",encoding='utf-8') as fa:
	flag = True
	lmda = 0
	f=0
	for line in fa:
		c0 = re.split(r'\s+|\s|: |, |:|,',line)
		c = [ele for ele in c0 if ele.strip()]
		c1 = np.array(c[:])
#		print(c1)
		if flag:
#			print(line)
#			print(c)
#			print(c1)#[-7])#.split(","))
			Ntotal = c1[-1].astype(int)
			sysdim = c1[-7].astype(int)
#			parti=np.zeros([sysdim*Ntotal])
#			modes=np.zeros([sysdim*Ntotal])
			k = 0
#			print(Ntotal,sysdim)
			flag = False

			continue
		if k==0:
			
			temp = c1[-1].astype(np.float64)
			temp2 = c1
			
			f += 1
			
#			if lmda == 0:
#				print('first frequency line=',temp2)
				
			partinv = 0
			partinum = 0
#			print(temp)
			for _ in range(1):
				next(fa)
			k = 1
			
#			print(k,lmda)
			continue

#		print(c1)
		n1 = c1.astype(np.float64)
		partinv = partinv + n1[-1]**4
#		partinum = partinum + n1[-1]**2

		
#		contriroot = np.append(contriroot,n1[-1])
		
		
			
#		eigvec = np.append(eigvec,n1)#,axis=0)
		k += 1
		if k==Ntotal+1:
			modes = np.append(modes,temp)#,axis=0)
			parti = np.append(parti,1.0/(Ntotal*partinv))
#			parti2 = np.append(parti2,(partinum**2)*1.0/(Ntotal*partinv))
			lmda += 1
#			print(lmda)
			k = 0
			if lmda >= sysdim*Ntotal:
#				print(f'last lamda before stopping')
				break
			for _ in range(1):
				next(fa)
#		print(k)

endread = time.time()
print('time endread=',endread)
print('last frequency line=',temp2)

n1 = np.sum(np.array([1 for i in parti if i <= localRatio]))
#n2 = np.sum(np.array([1 for i in parti2 if i <= 0.05]))
		
print('are the modes sorted?',is_sorted(modes))
is_sorted(modes)
print('number of modes =',len(modes))

#contrisum=np.sum(contriroot**2)
#print('Sum of all atomic contributions to the selected mode=',contrisum)

#########################################  Local modes

Nlow=2
Nfull=4

partilow=parti[modes<4]
modeslow=modes[modes<4]

idlow=np.argpartition(partilow, Nlow)
idlow=idlow[:Nlow]

idfull=np.argpartition(parti, Nfull)
idfull=idfull[:Nfull]

idall=np.append(idlow,idfull)

idall=idall[idall.argsort()]			# Since the modes are already sorted in eigvec.dat

print(f'idall after selection before sorting',idall)

##################################### Plotting 1
nrows = 2
ncolumns = 1
fig, axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(12,8))
#axes = axes.flatten()

j=0
axes[j].scatter(modes, parti, label=f'Number of highly localized modes (P$_\lambda$ <= {localRatio}) = {n1}', color='r', s=10)
axes[j].set_xlabel('Phonon mode Frequency, THz')
axes[j].set_ylabel(f'Participation Ratio (P$_\lambda$)')
axes[j].legend(loc='center right')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

j=1

axes[j].scatter(modes, parti, label=f'Plot limited to frequency range of {nulo} to {nuhi} THz \nand {ratiolo} to {ratiohi} participation ratios', color='r', s=10)
axes[j].scatter(modes[idall], parti[idall], label=f'Selected localized modes', color='b', s=10.5)
axes[j].set_xlabel('Phonon mode Frequency, THz')
axes[j].set_ylabel(f'Participation Ratio (P$_\lambda$)')
axes[j].legend(loc='center right')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
axes[j].set_xlim(nulo,nuhi)
axes[j].set_ylim(ratiolo,ratiohi)

plt.suptitle(f"Participation ratio")

plt.show()
fig.savefig('Participation ratio13.jpg',dpi=300, bbox_inches='tight')
plt.close()


##################################################### Atomic Contributions to selected local modes

startread = time.time()
print('started reading at time (s)=',startread)
with open("eigvec.dat","r",encoding='utf-8') as fa:
	flag = True
	lmda = 0
	f=0
	contriroot=np.array([])
	for line in fa:
		c0 = re.split(r'\s+|\s|: |, |:|,',line)
		c = [ele for ele in c0 if ele.strip()]
		c1 = np.array(c[:])
#		print(c1)
		if flag:
#			print(line)
#			print(c)
#			print(c1)#[-7])#.split(","))
			Ntotal = c1[-1].astype(int)
			sysdim = c1[-7].astype(int)
#			parti=np.zeros([sysdim*Ntotal])
#			modes=np.zeros([sysdim*Ntotal])
			k = 0
#			print(Ntotal,sysdim)
			flag = False

			continue
		if k==0:
			
			temp = c1[-1].astype(np.float64)
			temp2 = c1
			
			f += 1
			
			if lmda == 0:
				print(temp2)
				
#			print(temp)
			for _ in range(1):
				next(fa)
			k = 1
			
#			print(k,lmda)
			continue

#		print(c1)
		n1 = c1.astype(np.float64)

		if lmda in idall:
			contriroot = np.append(contriroot,n1[-1])

		k += 1
		if k==Ntotal+1:

			lmda += 1
#			print(lmda)
			k = 0
			if lmda >= sysdim*Ntotal:
				print(f'last lamda before stopping')
				break
			for _ in range(1):
				next(fa)
#		print(k)

endread = time.time()
print('time endread=',endread)
print('last frequency line=',temp2)

print(f'length of contributions root array = {len(contriroot)}')



#####################################

contri=contriroot**2

contri=contri.reshape(-1,1)

contri=contri.reshape(-1,Ntotal).transpose()

sumcontri=np.sum(contri,axis=0)

print(f'the sum of total atomic contributions to each mode= {np.average(sumcontri)} +/- {np.std(sumcontri)}')

#sumcontri=np.sum(contri,axis=1)

#print(f'the mean of total atomic contributions (transpose) to each mode= {np.average(sumcontri)} +/- {np.std(sumcontri)}')

print('shape of contribution array=',np.shape(contri))

basismap=np.append(basissorted,contri,axis=1)
print('shape of basis mapped array=',np.shape(basismap))
print('basis mapped array=\n',basismap)

###################################################### Selection of modes to map

c=5		#Number of columns in atomic mapping array before selected mode contributions
print(f'idall before adding {c}',idall)

idall=idall+c

print('shape of ID array=',np.shape(idall))
print(f'idall after adding {c}',idall)

###################################################### Grouping along increasing x



basismap=basismap[basismap[:, 1].argsort()]
slab=xslab

poscompare=xlo+slab
temparr=np.array([])

basismaptemp=np.array([])

for i in basismap:
#	print(np.shape(i))
	temp=i[1]
	if temp > poscompare:
		temparr=temparr.reshape(-1,len(idall)+c)
		totalcontri=np.sum(temparr[:,c:],axis=0)
		averagepos=np.average(temparr,axis=0)
		basismaptemp=np.append(basismaptemp,np.array([len(temparr[:,0])]))
		basismaptemp=np.append(basismaptemp,averagepos[1:4])
#		basismaptemp=np.append(basismaptemp,np.array([totalcontri]))
		basismaptemp=np.append(basismaptemp,np.array([totalcontri/len(temparr[:,0])]))
		temparr=np.array([])
		poscompare=poscompare+slab
	temparr=np.append(temparr,i)
basismapx=basismaptemp.reshape(-1,4+len(idall))

print(f'number of x slabs = {len(basismapx[:,0])}')

######################################################## along y


basismap=basismap[basismap[:, 2].argsort()]
slab=yslab


basismaptemp=np.array([])

poscompare=ylo+slab
temparr=np.array([])

for i in basismap:
#	print(np.shape(i))
	temp=i[2]
	if temp > poscompare:
		temparr=temparr.reshape(-1,len(idall)+5)
		totalcontri=np.sum(temparr[:,c:],axis=0)
		averagepos=np.average(temparr,axis=0)
		basismaptemp=np.append(basismaptemp,np.array([len(temparr[:,0])]))
		basismaptemp=np.append(basismaptemp,averagepos[1:4])
#		basismaptemp=np.append(basismaptemp,np.array([totalcontri]))
		basismaptemp=np.append(basismaptemp,np.array([totalcontri/len(temparr[:,0])]))
		temparr=np.array([])
		poscompare=poscompare+slab
	temparr=np.append(temparr,i)
basismapy=basismaptemp.reshape(-1,4+len(idall))

print(f'number of y slabs = {len(basismapy[:,0])}')

######################################################################

basismap=basismap[basismap[:, 3].argsort()]
slab=zslab

basismaptemp=np.array([])

poscompare=zlo+slab
temparr=np.array([])

basismaptemp=np.array([])

for i in basismap:
#	print(np.shape(i))
	temp=i[3]
	if temp > poscompare:
		temparr=temparr.reshape(-1,len(idall)+5)
		totalcontri=np.sum(temparr[:,c:],axis=0)
		averagepos=np.average(temparr,axis=0)
		basismaptemp=np.append(basismaptemp,np.array([len(temparr[:,0])]))
		basismaptemp=np.append(basismaptemp,averagepos[1:4])
#		basismaptemp=np.append(basismaptemp,np.array([totalcontri]))
		basismaptemp=np.append(basismaptemp,np.array([totalcontri/len(temparr[:,0])]))
		temparr=np.array([])
		poscompare=poscompare+slab
	temparr=np.append(temparr,i)
basismapz=basismaptemp.reshape(-1,4+len(idall))

print(f'number of z slabs = {len(basismapz[:,0])}')

############################################################ Plotting 2


nrows = 3
ncolumns = len(idall)

fig, axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(20,25))
axes = axes.flatten()

j=0
for i in range(len(idall)):
	k = i+4
	axes[j].plot(basismapx[:,1], basismapx[:,k],'--bo', label=f'{modes[idall[i]-5]: 0.2f}THz\nvs. x')
	axes[j].set_xlabel('Position, A')
	axes[j].set_ylabel(f'Phonon localization')
	axes[j].legend(loc='upper left')
	axes[j].adjustable='datalim'
	axes[j].set_aspect('auto')


	axes[j+len(idall)].plot(basismapy[:,2], basismapy[:,k], '--bo', label=f'{modes[idall[i]-5]: 0.2f}THz\nvs. y')
	axes[j+len(idall)].set_xlabel('Position, A')
	axes[j+len(idall)].set_ylabel(f'Phonon localization')
	axes[j+len(idall)].legend(loc='upper left')
	axes[j+len(idall)].adjustable='datalim'
	axes[j+len(idall)].set_aspect('auto')



	axes[j+2*len(idall)].plot(basismapz[:,3], basismapz[:,k],'--bo', label=f'{modes[idall[i]-5]: 0.2f}THz\nvs. z')
	axes[j+2*len(idall)].set_xlabel('Position, A')
	axes[j+2*len(idall)].set_ylabel(f'Phonon localization')
	axes[j+2*len(idall)].legend(loc='upper left')
	axes[j+2*len(idall)].adjustable='datalim'
	axes[j+2*len(idall)].set_aspect('auto')
	
	j +=1


plt.suptitle(f"Local phonon contributions")

plt.show()
fig.savefig('Local phonon contributions13.jpg',dpi=300, bbox_inches='tight')
plt.close()


		
