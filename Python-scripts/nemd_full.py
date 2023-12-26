import numpy as np
from numpy import polyfit
import matplotlib.pyplot as plt
import math
import re
import time

## Times are in picoseconds and distances are in angstroms

ts=0.0002
Estep=1000
Tstep=1000
thermosteps=1000
thickness=3.34
######################
initiallength=2400
slabsize=5
Texpected=300
deltaT=157
bathsize=25


meanstarttime=500
meanendtime= 1500
finitemeanflag = True
#####################
Nslabsexpected=int(initiallength/5)
runsteps=1500/ts		## Expected but may be shorter or longer - keep non restrictive
heatsteps=200/ts

avetime=20
avetime1=50
avetime2=100
avetime3=150
avetime4=200
avetime5=250
Ntemperatureplots=6		## Including 0 timestep
runningavTpoints=int(avetime/(ts*Tstep))
runningavEpoints=int(avetime/(ts*Estep))

runningavTpoints1=int(avetime1/(ts*Tstep))
runningavTpoints2=int(avetime2/(ts*Tstep))
runningavTpoints3=int(avetime3/(ts*Tstep))
runningavTpoints4=int(avetime4/(ts*Tstep))
runningavTpoints5=int(avetime5/(ts*Tstep))

runningavTpointswide = runningavTpoints2
avetimewide = avetime2

distancefromsource=60

Tfname = 'tmp.dat'
Efname = 'energy.dat'

##################################	Processing energy data
##################################

Edataori = np.loadtxt(Efname, delimiter=' ', comments='#', skiprows=0, max_rows=None)

#Esource, Esink, Ecouple, Econserve = Edata[:,0], Edata[:,1], Edata[:,2], Edata[:,3] 

Elabels = ['Energy to source','Energy to sink','Energy net thermostat','Energy conserved']
Estdlabels = ['Energy to source std','Energy to sink std','Energy net thermostat std','Energy conserved std']

Eavlen = int((len(Edataori[:,0])//runningavEpoints)*runningavEpoints)

#E100, E100std = np.average(Esource[:Eavlen].reshape(-1,runningavEpoints),axis=1), np.std(Esource[:Eavlen].reshape(-1,runningavEpoints),axis=1) 

#Esource, Esourcestd = np.average(Esource[:Eavlen].reshape(-1,runningavEpoints),axis=1), np.std(Esource[:Eavlen].reshape(-1,runningavEpoints),axis=1) 

#Esink, Esinkstd = np.average(Esink[:Eavlen].reshape(-1,runningavEpoints),axis=1), np.std(Esink[:Eavlen].reshape(-1,runningavEpoints),axis=1)

#Ecouple, Ecouplestd = np.average(Ecouple[:Eavlen].reshape(-1,runningavEpoints),axis=1), np.std(Ecouple[:Eavlen].reshape(-1,runningavEpoints),axis=1)

#Econserve, Econservestd = np.average(Econserve[:Eavlen].reshape(-1,runningavEpoints),axis=1), np.std(Econserve[:Eavlen].reshape(-1,runningavEpoints),axis=1)

#timeE = np.arange(avetime,(len(Esource)+1)*avetime,avetime)

#E100 = [[Esource100,Esource100std][]]

Edatatranspose = Edataori.transpose()

Edata, Edatastd = np.average(Edatatranspose[:,:Eavlen].reshape(7,-1,runningavEpoints),axis=2), np.std(Edatatranspose[:,:Eavlen].reshape(7,-1,runningavEpoints),axis=2)

Edata, Edatastd = Edata.transpose(), Edatastd.transpose()

timeE = np.arange(avetime,(len(Edata)+1)*avetime,avetime)

print(f'shape of Edata = {Edata.shape},\nshape of timeE = {timeE.shape}\n')

#######################################		Finding slope for k calculation

if (timeE[-1] >= meanendtime) and finitemeanflag and (timeE[-1] > meanstarttime):
	meanmask = (timeE < meanendtime) & (timeE >= meanstarttime)
else:
	print('\nWarning: the simulation might not have run for sufficient duration.\n')
	meanmask = (timeE >= 0)
timeEnew = timeE[meanmask]
Edatanew = Edata[meanmask]
xlength = Edatanew[-1,-3]
ylength = Edatanew[-1,-2]
zlength = thickness
	
dEdt , dEdtstd = np.polyfit(timeEnew,(Edatanew[:,1]-Edatanew[:,0])/2,1,cov=1)
dEdtmetric, dEdtstdmetric = dEdt[0]*1.602176634, dEdtstd[0][0]*1.602176634**2		#converting to e-19 J/ps
	
print(f'xlength is {xlength}\nylength is {ylength}\nzlength is \
	{zlength}\ndE/dt is {dEdt[0]} +/- {np.sqrt(dEdtstd[0][0])} eV/ps')
	
deltax = xlength - 2*initiallength*bathsize/xlength - 2*initiallength*slabsize/xlength
Tgradientexpected = 2*deltaT/deltax
flux = dEdtmetric*1e-19/(ylength*zlength)
	
k = (flux/Tgradientexpected)*1e+22 
	
#print(f'{flux}')
	
print(f'Expected temperature gradient over {deltax: 0.2f} A \
	is {Tgradientexpected: .5E} K/A\nHeat flux is {flux: .4E} J/ps.A^2\
	\nk is {k: .6E} W/m.K\n\n\n\n')
	
#	exit()
	

#######################################		Plotting energy convergence
nrows = 2
ncolumns = 4
j = 0
fig,axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(ncolumns*3.5,nrows*3))
axes=axes.flatten()
#print(len(axes))
for j in range(ncolumns):

	axes[j].plot(timeE, Edata[:,j], label=Elabels[j],color='blue',linestyle=':' )
	axes[j].set_xlabel('Time, ps')
	axes[j].set_ylabel(f'Energy, eV')
	axes[j].legend(loc='upper right')
	axes[j].adjustable='datalim'
	axes[j].set_aspect('auto')
	
#	print(j)

	axes[j+ncolumns].plot(timeE, Edatastd[:,j], label=Estdlabels[j],color='blue',linestyle=':'  )
	axes[j+ncolumns].set_xlabel('Time, ps')
	axes[j+ncolumns].set_ylabel(f'Energy, eV')
	axes[j+ncolumns].legend(loc='center right')
	axes[j+ncolumns].adjustable='datalim'
	axes[j+ncolumns].set_aspect('auto')

	j += 1

Edatafullmean = np.average(Edatatranspose, axis=1)
print(f'Mean E net thermostat = {Edatafullmean[2]}')
j = 2	
axes[j].axhline(y=Edatafullmean[2],color='black')

figname= Efname + ' ' + str(int(avetime)) + 'ps' + '.png'
plt.suptitle(f"Energy Convergence")
plt.show()
fig.savefig(figname,dpi=300, bbox_inches='tight')
plt.close()

################################## 	Read slab temperature data file
##################################

readtime=[]
temperatures=[]
coordinates=[]
Ncount=[]
#time=np.append(time,0)

print('time start tmp.dat read=',time.time())

with open(Tfname,"r",encoding='utf-8') as fa:
	for _ in range(3):
		next(fa)
	islab=0
	startflag=1
	totsteps=0
	for line in fa:
		if startflag:
			c0 = re.split(r'\s+|./|\s|: |, |:|,',line)
			c = [ele for ele in c0 if ele.strip()]
#			c = line.split(" ")
			c1 = np.array(c[:])
			readtime.append(c1[0].astype(int)*ts)
			totsteps = c1[0].astype(int)
			
			Nslabs=c1[1].astype(int)
			startflag=0
			continue
#			if slabs==slabsExpected+1:
#				slabsFlag = True
#			elif slabs==slabsExpected:
#				slabsFlag = False
		if islab == Nslabs:
			c0 = re.split(r'\s+|./|\s|: |, |:|,',line)
			c = [ele for ele in c0 if ele.strip()]
#			c = line.split(" ")
			c1 = np.array(c[:])
			readtime.append(c1[0].astype(int)*ts)
			totsteps = c1[0].astype(int)
			islab=0
			continue
		
		c0 = re.split(r'\s+|./|\s|: |, |:|,',line)
		c = [ele for ele in c0 if ele.strip()]
#		c = line.split(" ")
		c1 = np.array(c[:])
#		print(islab)
#		print(c1)
		coordinates.append(c1[1].astype(float))
		temperatures.append(c1[3].astype(float))
		Ncount.append(round(float(c1[2])))	#.astype(int))
		
		islab += 1
		
#endread = time.time()
print('time end tmp.dat read=',time.time())
		
readtime = np.array(readtime)
coordinates = np.array(coordinates).reshape(-1,Nslabs)
temperatures = np.array(temperatures).reshape(-1,Nslabs)
Ncount = np.array(Ncount).reshape(-1,Nslabs)

print(f'total steps = {totsteps}, total time = {readtime[-1]}, \
number of slaps = {Nslabs} \nlast slab index value = {islab}, \
shape of time array = {readtime.shape}, shape of \
temperature/coordinates/Ncount array = {temperatures.shape}')

if Nslabs != Nslabsexpected:
	mask = np.ones(temperatures.shape, dtype=bool)
	mask[:,[Nslabs-1]] = False
#	print(temperatures)
	T = temperatures[mask].reshape(-1,Nslabsexpected)
	C = coordinates[mask].reshape(-1,Nslabsexpected)
	N = Ncount[mask].reshape(-1,Nslabsexpected)
#	print(T)
else:
	T = temperatures
	C = coordinates
	N = Ncount

#print(f'shape of T = {T.shape}')

Cinner = C[:,1:len(C[0,:])-1]
Tinner = T[:,1:len(T[0,:])-1]
Ninner = N[:,1:len(N[0,:])-1]

print(T)
print(Tinner)

print(f'shape of T = {T.shape}\nshape of Tinner = {Tinner.shape}\n')

#########################	Averaging the Temperatures

Tavlen = int((len(Tinner[:,0])//runningavTpoints)*runningavTpoints)

Ttranspose = Tinner.transpose()

print('time start T average=',time.time())

Taverage, Taveragestd = np.average(Ttranspose[:,:Tavlen].reshape(len(Ttranspose[:,0]),-1,runningavTpoints),axis=2), np.std(Ttranspose[:,:Tavlen].reshape(len(Ttranspose[:,0]),-1,runningavTpoints),axis=2)

print('time end T average=',time.time())

Taverage, Taveragestd = Taverage.transpose(), Taveragestd.transpose()

timeT = np.arange(avetime,(len(Taverage)+1)*avetime,avetime)

print(f'shape of Taverage = {Taverage.shape},\nshape of timeT = {timeT.shape}\n')

##########################		Slope dT/dx

slabsfromend = int(bathsize/slabsize + distancefromsource/slabsize)

Clinear = Cinner[:,1:len(Cinner[0,:])-1]
Tlinear = Tinner[:,1:len(Tinner[0,:])-1]
Nlinear = Ninner[:,1:len(Ninner[0,:])-1]

slope, dslope = [], []

print('time start slope calculation=',time.time())

for i in range(len(Clinear)):
	slopei , dslopei = np.polyfit(Clinear[i].transpose(),Tlinear[i].transpose(),1,cov=1)
	slope.append(slopei)
	dslope.append(dslopei)
	
print('time end slope calculation=',time.time())

slope = np.array(slope)
dslope = np.array(dslope)

print(f'shape of slope = {slope.shape}\nshape of dslope = {dslope.shape}\n')

slope = slope[:,0]
dslope = np.sqrt(dslope[:,0,0])
print(f'shape of slope = {slope.shape}\nshape of dslope = {dslope.shape}\n')

##########################		average slope and mean 'T-average'&'T-std'
print('time start slope averaging=',time.time())
slopeaverage, slopeaveragestd = np.average(slope[:Tavlen].reshape(-1,runningavTpoints),axis=1), np.std(slope[:Tavlen].reshape(-1,runningavTpoints),axis=1)
print('time end slope averaging=',time.time())

mean_Taverage, mean_Taveragestd = np.average(Taverage, axis=1), np.average(Taveragestd, axis=1)

print(f'shape of slopeaverage = {slopeaverage.shape}\nshape of slopeaveragestd =\
 {slopeaveragestd.shape}\nshape of mean_Taverage = {mean_Taverage.shape}\nshape of\
  mean_Taveragestd = {mean_Taveragestd.shape}\n')

#######################################		Plotting profile convergence

nrows = 2
ncolumns = 2
j = 0
fig,axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(ncolumns*3.5,nrows*3))
axes=axes.flatten()

axes[j].plot(timeT, mean_Taverage, label='mean_Taverage',color='blue',linestyle=':')
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel(f'Temperature, K')
axes[j].legend(loc='upper right')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
	
#print(j)
axes[j+ncolumns].plot(timeT, mean_Taveragestd, label='mean_Taveragestd',color='blue',linestyle=':')
axes[j+ncolumns].set_xlabel('Time, ps')
axes[j+ncolumns].set_ylabel(f'Temperature, K')
axes[j+ncolumns].legend(loc='center right')
axes[j+ncolumns].adjustable='datalim'
axes[j+ncolumns].set_aspect('auto')

j += 1
axes[j].plot(timeT, slopeaverage, label='slope average',color='blue',linestyle=':')
axes[j].set_xlabel('Time, ps')
axes[j].set_ylabel(f'Temperature gradient, K/A')
axes[j].legend(loc='upper right')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
	
#print(j)

axes[j+ncolumns].plot(timeT, slopeaveragestd, label='slope averagestd',color='blue',linestyle=':')
axes[j+ncolumns].set_xlabel('Time, ps')
axes[j+ncolumns].set_ylabel(f'Temperaturegradient, K/A')
axes[j+ncolumns].legend(loc='center right')
axes[j+ncolumns].adjustable='datalim'
axes[j+ncolumns].set_aspect('auto')


mean_Taveragefullmean = np.average(mean_Taverage)#, axis=1)
print(f'Mean temperature full average = {mean_Taveragefullmean}')
j = 0	
axes[j].axhline(y=mean_Taveragefullmean,color='black')

slopeaveragefullmean = np.average(slopeaverage)#, axis=1)
print(f'slope full average = {slopeaveragefullmean}')
j = 1	
axes[j].axhline(y=slopeaveragefullmean,color='black')


figname= Tfname + ' ' + str(int(avetime)) + 'ps' + '.png'

plt.suptitle(f"Temperature profile convergence")
plt.show()
fig.savefig(figname,dpi=300, bbox_inches='tight')
plt.close()


#######################################		Average profile

dxwithbath = xlength - 2*initiallength*slabsize/xlength

avestartindex=int(meanstarttime/(ts*Tstep)-1)
aveendindex=int(meanendtime/(ts*Tstep)-1)

netTprofile=np.average(Tinner[avestartindex:aveendindex],axis=0)
print(avestartindex,aveendindex,netTprofile)
netCprofile=np.average(Cinner[avestartindex:aveendindex],axis=0)
netTprofilestd=np.std(Tinner[avestartindex:aveendindex],axis=0)
netCprofilestd=np.std(Cinner[avestartindex:aveendindex],axis=0)

meannetTprofilestd=np.average(netTprofilestd)

netTdiff=np.max(netTprofile)-np.min(netTprofile)

Tgradientacrossbath = netTdiff/dxwithbath
krealacrossbath = (flux/Tgradientacrossbath)*1e+22 


netlinearTprofile=np.average(Tlinear[avestartindex:aveendindex],axis=0)
netlinearCprofile=np.average(Clinear[avestartindex:aveendindex],axis=0)
netlinearTprofilestd=np.std(Tlinear[avestartindex:aveendindex],axis=0)
netlinearCprofilestd=np.std(Clinear[avestartindex:aveendindex],axis=0)

netlinearTdiff=np.max(netlinearTprofile)-np.min(netlinearTprofile)

#netCdiff=np.max(netCprofile)-np.min(netCprofile)

Tgradientreal = netlinearTdiff/deltax
kreal = (flux/Tgradientreal)*1e+22 

print(f'\n\n\n\nxlength is {xlength}\nylength is {ylength}\nzlength is \
	{zlength}\ndE/dt is {dEdt[0]} +/- {np.sqrt(dEdtstd[0][0])} eV/ps \
	\nHeat flux is {flux: .4E} J/ps.A^2\
	\nExpected T diff = {2*deltaT}\
	\nActual T diff NOT considering bath = {netlinearTdiff} \
	\nActual T diff considering bath = {netTdiff} \
	\n\nExpected temperature gradient over {deltax: 0.2f} A is {Tgradientexpected: .5E} K/A \
	\nActual temperature gradient BETWEEN baths over {deltax: 0.2f} A is {Tgradientreal: .5E} K/A \
	\nActual temperature gradient across bath over {dxwithbath: 0.2f} A is {Tgradientacrossbath: .5E} K/A \
	\n\nk from delta T expected is {k: .6E} W/m.K\
	\nk from actual delta T BETWEEN baths is {kreal: .6E} W/m.K\
	\nk across/including baths is {krealacrossbath: .6E} W/m.K\n\n\n\n')

fa=open('k_log.txt','w')
fa.write(f'\n\n\n\nxlength is {xlength}\nylength is {ylength}\nzlength is \
	{zlength}\ndE/dt is {dEdt[0]} +/- {np.sqrt(dEdtstd[0][0])} eV/ps \
	\nHeat flux is {flux: .4E} J/ps.A^2\
	\nExpected T diff = {2*deltaT}\
	\nActual T diff NOT considering bath = {netlinearTdiff} \
	\nActual T diff considering bath = {netTdiff} \
	\n\nExpected temperature gradient over {deltax: 0.2f} A is {Tgradientexpected: .5E} K/A \
	\nActual temperature gradient BETWEEN baths over {deltax: 0.2f} A is {Tgradientreal: .5E} K/A \
	\nActual temperature gradient across bath over {dxwithbath: 0.2f} A is {Tgradientacrossbath: .5E} K/A \
	\n\nk from delta T expected is {k: .6E} W/m.K\
	\nk from actual delta T BETWEEN baths is {kreal: .6E} W/m.K\
	\nk across/including baths is {krealacrossbath: .6E} W/m.K\n\n\n\n')
fa.close()

nrows = 1
ncolumns = 2
j = 0
fig,axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(ncolumns*7,nrows*6))
#axes=axes.flatten()

axes[j].errorbar(netCprofile, netTprofile, yerr=netTprofilestd, label=str(int(meanstarttime))+'ps to '+str(int(meanendtime))+'ps averaged T profile',color='blue')
axes[j].set_xlabel('Position, A')
axes[j].set_ylabel(f'Temperature, K')
axes[j].legend(loc='center right')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

j += 1
axes[j].plot(netCprofile, netTprofile, label=str(int(meanstarttime))+'ps to '+str(int(meanendtime))+'ps averaged T profile',color='blue')
axes[j].set_xlabel('Position, A')
axes[j].set_ylabel(f'Temperature, K')
axes[j].legend(loc='center right')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')

figname= Tfname + ' ' + str(int(meanstarttime))+'ps to '+str(int(meanendtime)) + 'ps with and without error' + '.png'

plt.suptitle(f"Temperature profiles through k calculation-averaging time")
plt.show()
fig.savefig(figname,dpi=300, bbox_inches='tight')
plt.close()

#######################################		 Wider Temperature average

Tavlenwide = int((len(Tinner[:,0])//runningavTpointswide)*runningavTpointswide)

Ttransposewide = Tinner.transpose()

print('time start T average=',time.time())

Taveragewide, Taveragestdwide = np.average(Ttransposewide[:,:Tavlenwide].reshape(len(Ttransposewide[:,0]),-1,runningavTpointswide),axis=2), np.std(Ttransposewide[:,:Tavlenwide].reshape(len(Ttransposewide[:,0]),-1,runningavTpointswide),axis=2)

print('time end T average=',time.time())

Taveragewide, Taveragestdwide = Taveragewide.transpose(), Taveragestdwide.transpose()

timeTwide = np.arange(avetimewide,(len(Taveragewide)+1)*avetimewide,avetimewide)


Cavlenwide = int((len(Cinner[:,0])//runningavTpointswide)*runningavTpointswide)

Ctransposewide = Cinner.transpose()

print('time start C average=',time.time())

Caveragewide, Caveragestdwide = np.average(Ctransposewide[:,:Cavlenwide].reshape(len(Ctransposewide[:,0]),-1,runningavTpointswide),axis=2), np.std(Ctransposewide[:,:Cavlenwide].reshape(len(Ctransposewide[:,0]),-1,runningavTpointswide),axis=2)

print('time end C average=',time.time())

Caveragewide, Caveragestdwide = Caveragewide.transpose(), Caveragestdwide.transpose()

print(f'shape of Caveragewide = {Caveragewide.shape}\n\n')


#######################################		Plotting profiles

di = int(len(timeTwide)//Ntemperatureplots)

print('di=',di)

nrows = min(di, 3)
ncolumns = Ntemperatureplots
j = 0
i = 0
fig,axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(ncolumns*3.5,nrows*3))
axes=axes.flatten()

for i in range(nrows):
	axes[j+i*ncolumns].errorbar(Caveragewide[i], Taveragewide[i], yerr=Taveragestdwide[i], label='at '+str(timeTwide[i])+'ps: '+str(avetimewide)+'ps averaged T\'s',color='blue')
	axes[j+i*ncolumns].set_xlabel('Position, A')
	axes[j+i*ncolumns].set_ylabel(f'Temperature, K')
	axes[j+i*ncolumns].legend(loc='center right')
	axes[j+i*ncolumns].adjustable='datalim'
	axes[j+i*ncolumns].set_aspect('auto')

i = max(di-3,di)
j += 1
for k in range(Ntemperatureplots-1):
	for l in range(nrows):
		axes[j+l*ncolumns].errorbar(Caveragewide[i+l], Taveragewide[i+l], yerr=Taveragestdwide[i+l], label='at '+str(timeTwide[i+l])+'ps: '+str(avetimewide)+'ps averaged T\'s',color='blue')
		axes[j+l*ncolumns].set_xlabel('Position, A')
		axes[j+l*ncolumns].set_ylabel(f'Temperature, K')
		axes[j+l*ncolumns].legend(loc='center right')
		axes[j+l*ncolumns].adjustable='datalim'
		axes[j+l*ncolumns].set_aspect('auto')

	i += di-nrows+1 
	j += 1

figname= Tfname + ' ' + str(int(avetimewide)) + 'ps with error' + '.png'

plt.suptitle(f"Temperature profile convergence")
plt.show()
fig.savefig(figname,dpi=300, bbox_inches='tight')
plt.close()
######################### No error bar
j = 0
i = 0
fig,axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(ncolumns*3.5,nrows*3))
axes=axes.flatten()

for i in range(nrows):
	axes[j+i*ncolumns].plot(Caveragewide[i], Taveragewide[i], label='at '+str(timeTwide[i])+'ps: '+str(avetimewide)+'ps averaged T\'s',color='blue')
	axes[j+i*ncolumns].set_xlabel('Position, A')
	axes[j+i*ncolumns].set_ylabel(f'Temperature, K')
	axes[j+i*ncolumns].legend(loc='center right')
	axes[j+i*ncolumns].adjustable='datalim'
	axes[j+i*ncolumns].set_aspect('auto')

i = max(di-3,di)
j += 1
for k in range(Ntemperatureplots-1):
	for l in range(nrows):
		axes[j+l*ncolumns].plot(Caveragewide[i+l], Taveragewide[i+l], label='at '+str(timeTwide[i+l])+'ps: '+str(avetimewide)+'ps averaged T\'s',color='blue')
		axes[j+l*ncolumns].set_xlabel('Position, A')
		axes[j+l*ncolumns].set_ylabel(f'Temperature, K')
		axes[j+l*ncolumns].legend(loc='center right')
		axes[j+l*ncolumns].adjustable='datalim'
		axes[j+l*ncolumns].set_aspect('auto')

	i += di-nrows+1 
	j += 1

figname= Tfname + ' ' + str(int(avetimewide)) + 'ps' + '.png'

plt.suptitle(f"Temperature profile convergence")
plt.show()
fig.savefig(figname,dpi=300, bbox_inches='tight')
plt.close()






exit()



