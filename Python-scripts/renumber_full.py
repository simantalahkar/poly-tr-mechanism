import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys

filename = 'trimmed.data'
atoms = 749757
#atoms_original = 1021937
#atoms_deleted = atoms_original - atoms

coord2d = np.loadtxt(filename, delimiter=' ', comments='#', skiprows=9, max_rows=atoms)
coord2d = coord2d[coord2d[:,0].argsort()]
coord2d[:,0] = np.arange(1,atoms+1).astype(int)
np.savetxt('coord.txt', coord2d, fmt='%i %i %i %.3e %.18e %.18e %.18e %i %i %i', delimiter=' ')
del coord2d


vel2d = np.loadtxt(filename, delimiter=' ', comments='#', skiprows=9+atoms+3, max_rows=atoms)
vel2d = vel2d[vel2d[:,0].argsort()]
vel2d[:,0] = np.arange(1,atoms+1).astype(int)
np.savetxt('vel.txt', vel2d, fmt='%i %.18e %.18e %.18e', delimiter=' ')
del vel2d



