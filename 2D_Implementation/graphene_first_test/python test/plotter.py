
# import required libraries
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt


# Read H5 file
f = h5.File("OPEN_BC_RichtmyerCLASSTEST_2D_Graphene_DS.h5", "r")

density_group = f['Data/Density']
#velocity_x_group = f['Data/VelocityX']
#velocity_y_group = f['Data/VelocityX']

#plt.figure()
#plt.clim(0.2, 1.8);
#plt.colorbar()


for key in density_group.keys():
	data_density=density_group[key]
	fig_name = key + '.jpeg'
	plt.pcolormesh(data_density,cmap='Spectral')
#	plt.clim(0.2, 2.0);
	plt.clim(0.3, 1.9);
	cbar=plt.colorbar()
	cbar.set_label('n/n0')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.savefig(fig_name,quality=60)
	plt.close("all")


 



f.close()
