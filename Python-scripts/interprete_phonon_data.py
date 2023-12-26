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

def read_phonon_data(path):
	dataset = []
	names = []
	ratios = []
	n_modesarr = []
	for fname in glob.glob(path):
		c0 = re.split(r'\s+|.txt|.custom|.data|./|\s|: |, |:|,',fname)
		c = [ele for ele in c0 if ele.strip()]
		ratio = float(c[1])
		n_modes = float(c[2])
	
		data = np.loadtxt(fname, delimiter=' ', comments='#', skiprows=0, max_rows=None)
		names.append(c[0])
		dataset.append(data)
		ratios.append(ratio)
		n_modesarr.append(n_modes)
	return names, np.array(dataset), np.array(ratios), np.array(n_modesarr)

###################### Mapping limits for local participation plot
nulo=0
nuhi=20
ratiolo=0
ratiohi=0.8
#################################################################
######################################### Reading the eigenvector file

pathspatial = './spatial*.custom'
spatialfiles, spatialset, propagate_ratios, npropagating_modes = read_phonon_data(pathspatial)
	
pathparti = './participation*.data'
partifiles, modesarr, local_ratios, nlocal_modes = read_phonon_data(pathparti)

#print(np.shape(modesarr),np.shape(modesarr[0,:,1]))

for i in range(len(partifiles)):
	print(i)
	
	freq=modesarr[i,:,0]
	parti=modesarr[i,:,1]
#	meanfreq = np.mean(modesarr[i,:,0])
	medianfreq = np.median(freq)
	freq25 = np.percentile(freq, 25)
	freq75 = np.percentile(freq, 75)
	
	meanparti = np.mean(parti)
	medianparti = np.median(parti)
	parti25 = np.percentile(parti, 25)
#	print(freq<=freq25)
	meanparti25 = np.mean(parti[freq<=freq25])
	meanparti25_75 = np.mean(parti[(freq>=freq25)&(freq<=freq75)])
	meanparti75 = np.mean(parti[freq>=freq75])
	minparti = np.min(parti)
	meanfreq20THz = np.mean(freq[freq<25])
	meanparti20THz = np.mean(parti[freq<25])
		
	f=open('participation_stats.dat','w')
	f.write(f'median_frequency {medianfreq}\nfreq25 {freq25}\nfreq75 {freq75}\nmeanparti {meanparti}\nmedianparti {medianparti}\nparti25 {parti25}\nmeanparti25 {meanparti25}\nmeanparti25_75 {meanparti25_75}\nmeanparti75 {meanparti75}\nminparti {minparti}\nmeanparti20THz {meanparti20THz}\nmeanfreq20THz {meanfreq20THz}')
	f.close()
	
	
	fname= partifiles[i]+' '+str(local_ratios[i])+' '+str(nlocal_modes[i])+'.dat'
	
	nrows = 2
	ncolumns = 1
	fig,axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(9,6))

	j=0
	axes[j].scatter(modesarr[i,:,0], modesarr[i,:,1], label=f'Number of localized modes (P$_\lambda$ <= {local_ratios[i]}) = {nlocal_modes[i]}', color='r', s=8)
	axes[j].set_xlabel('Phonon mode Frequency, THz')
	axes[j].set_ylabel(f'Participation Ratio (P$_\lambda$)')
	axes[j].legend(loc='center right')
	axes[j].adjustable='datalim'
	axes[j].set_aspect('auto')

	j=1
	axes[j].scatter(modesarr[i,:,0], modesarr[i,:,1], label=f'Plot limited to frequency range of {nulo} to {nuhi} THz \nand {ratiolo} to {ratiohi} participation ratios', color='r', s=8)
	axes[j].set_xlabel('Phonon mode Frequency, THz')
	axes[j].set_ylabel(f'Participation Ratio (P$_\lambda$)')
	axes[j].legend(loc='center right')
	axes[j].adjustable='datalim'
	axes[j].set_aspect('auto')
	axes[j].set_xlim(nulo,nuhi)
	axes[j].set_ylim(ratiolo,ratiohi)
	
	figname= partifiles[i]+' '+str(local_ratios[i])+' '+str(nlocal_modes[i])+'.png'
	
	plt.suptitle(f"Participation ratio")
	plt.show()
	fig.savefig(figname,dpi=300, bbox_inches='tight')
	plt.close()
	
	nrows = 2
	ncolumns = 1
	
	fig,axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(9,6))
	j=0
	axes[j].boxplot(modesarr[i,:,0], meanline=True, autorange=False, vert=False, showmeans=True) 
#	axes[j].legend(loc='upper right')
	axes[j].set_xlabel('Phonon frequencies (THz)')
	
	
	j=1
	axes[j].boxplot(modesarr[i,:,1], meanline=True, autorange=False, vert=False, showmeans=True) 
#	axes[j].legend(loc='upper right')
	axes[j].set_xlabel('Participation Ratios')
	
	figname= partifiles[i]+'_boxplot'+' '+str(local_ratios[i])+' '+str(nlocal_modes[i])+'.png'
	
	plt.suptitle(f"Statistical distribution of phonon energies and frequencies")
	plt.show()
	fig.savefig(figname,dpi=300, bbox_inches='tight')
	plt.close()


#for i in range(len(spatialfiles)):
#	Xs, Ys, energy = spatialset[i,:,2], spatialset[i,:,3], spatialset[i,:,6]


	







