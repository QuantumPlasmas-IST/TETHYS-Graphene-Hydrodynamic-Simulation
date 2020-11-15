# import required libraries
import glob
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from pydmd import DMD

for filename in glob.glob("*.h5"):
	print('Reading file:')
	print(filename)
	f = h5.File(filename, "r")
	data_group =  f['Data']
	density_group = f['Data/Density']
	label = filename.replace('hdf5_2D_','').replace('.h5','')
	print('Saving snapshots:')
	snapshots = [np.array(density_group[key]) for key in density_group.keys()]
	dmd = DMD(exact=True, opt=False)
	print('Fitting snapshots:')
	dmd.fit(snapshots)
	dmd.plot_eigs()
f.close()

