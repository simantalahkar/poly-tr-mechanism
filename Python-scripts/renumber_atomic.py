import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys

filename = 'pos.data'
atoms = 518238
#atoms_original = 1021937
#atoms_deleted = atoms_original - atoms

coord2d = np.loadtxt(filename, delimiter=' ', comments='#', skiprows=9, max_rows=atoms)
#coord2d = coord2d[coord2d[:,0].argsort()]
coord2d[:,0] = np.arange(1,atoms+1).astype(int)
np.savetxt('coord.txt', coord2d, fmt='%i %i %.18e %.18e %.18e %i %i %i', delimiter=' ')
del coord2d





