import numpy as np
from numpy import polyfit
import matplotlib.pyplot as plt
import math
import re
import time

N = 10
max = 15
min = 1
xs = np.linspace(min, max, N, endpoint = True)
ys = np.linspace(min, max, N, endpoint = True)

#print('xs=\n',xs)
#	print('xcoordinates=',xs,'\nycoordinates=',ys)

Xs, Ys = np.meshgrid(xs, ys)  #shape: rn x cn

Xsflat = Xs.reshape(-1)
Ysflat = Ys.reshape(-1)

#print('Xsflat=\n',Xsflat)

pairs = np.array([a for a in zip(Xsflat,Ysflat)])

#print('pairs=\n',pairs)

mask = (pairs[:,0]-2*pairs[:,1])>-max

#print('mask=',mask)

Xsflatnew = Xsflat[mask]
Ysflatnew = Ysflat[mask]
sizes1 = 55 - 3.5*Xsflatnew + Ysflatnew**1.9
#sizes2 = Xsflatnew**1.3 + 15*Ysflatnew
sizes2 = 10+12*Ysflatnew
sizes3 = 2*(sizes2-sizes1)
sizes3opp = -sizes3
maskmech1 = sizes3opp>0
maskmech2 = sizes3>0
Xsflat1 = Xsflatnew[maskmech1]
Ysflat1 = Ysflatnew[maskmech1]
Xsflat2 = Xsflatnew[maskmech2]
Ysflat2 = Ysflatnew[maskmech2]
sizes3mech1 = sizes3opp[maskmech1]
sizes3mech2 = sizes3[maskmech2]


#exit()

#######################################		Plotting 
nrows = 1
ncolumns = 3
fig,axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(ncolumns*3.5,nrows*3))
axes=axes.flatten()

j = 0
axes[j].scatter(Xsflatnew, Ysflatnew, label='large to small grains',c='red',s=sizes1)
axes[j].set_xlabel('Length (arbitrary u.)')
axes[j].set_ylabel('Grain size asymmetry (arbitrary u.)')
axes[j].set_ylim([min-(max-min)/N,max])
axes[j].set_xlim([min-(max-min)/N,max+(max-min)/N])
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
axes[j].set_xticklabels([])
axes[j].set_yticklabels([])
axes[j].set_xticks([])
axes[j].set_yticks([])
axes[j].title.set_text('Mechanism 1 Strength')



j += 1
axes[j].scatter(Xsflatnew, Ysflatnew, label='small to large grains',c='blue',s=sizes2)
axes[j].set_xlabel('Length (arbitrary u.)')
#axes[j].set_ylabel('TSD')
axes[j].set_ylim([min-(max-min)/N,max])
axes[j].set_xlim([min-(max-min)/N,max+(max-min)/N])
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
axes[j].set_xticklabels([])
axes[j].set_yticklabels([])
axes[j].set_xticks([])
axes[j].set_yticks([])
axes[j].title.set_text('Mechanism 2 Strength')


j += 1
axes[j].scatter(Xsflat2, Ysflat2, label='small to large grains',c='blue',s=sizes3mech2)
axes[j].scatter(Xsflat1, Ysflat1, label='large to small grains',c='red',s=sizes3mech1)
axes[j].set_xlabel('Length (arbitrary u.)')
#axes[j].set_ylabel('TSD')
axes[j].set_ylim([min-(max-min)/N,max])
axes[j].set_xlim([min-(max-min)/N,max+(max-min)/N])
axes[j].legend(loc='upper left')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')
axes[j].set_xticklabels([])
axes[j].set_yticklabels([])
axes[j].set_xticks([])
axes[j].set_yticks([])
axes[j].title.set_text('Resultant Thermal Rectification')


figname= 'mechanisms3' + '.png'
#plt.suptitle(f"Energy Convergence")
plt.show()
fig.savefig(figname,dpi=300, bbox_inches='tight')
plt.close()



