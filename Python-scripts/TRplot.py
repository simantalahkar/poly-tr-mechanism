import numpy as np
from numpy import polyfit
import matplotlib.pyplot as plt
import math
import re


nrows = 1
ncolumns = 2
j = 0
fig,axes = plt.subplots(nrows,ncolumns,squeeze=True,constrained_layout=True,figsize=(ncolumns*7,nrows*6))
#axes=axes.flatten()

axes[j].scatter(netCprofile, netTprofile, yerr=netTprofilestd, label=str(int(meanstarttime))+'ps to '+str(int(meanendtime))+'ps averaged T profile',color='blue', elinewidth=0.5)
axes[j].set_xlabel('Position, A')
axes[j].set_ylabel(f'Temperature, K')
axes[j].legend(loc='lower right')
axes[j].adjustable='datalim'
axes[j].set_aspect('auto')



figname= 'average plus error with break' + '.png'

plt.suptitle(f"Rectification in gradient grain density polycrystalline graphene")
plt.show()
fig.savefig(figname,dpi=300, bbox_inches='tight')
plt.close()