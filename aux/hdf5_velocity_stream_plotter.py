# import required libraries
import glob
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

# Read H5 file
for filename in glob.glob("*.h5"):
	print('Reading file:')
	print(filename)
	f = h5.File(filename, "r")
	data_group =  f['Data']
	density_group = f['Data/Density']
	velocity_x_group = f['Data/VelocityX']
	velocity_y_group = f['Data/VelocityY']
	label = filename.replace('hdf5_2D_','').replace('.h5','')

	X = np.linspace(0, 1.0, num=data_group.attrs.get('Number of spatial points x')[0])
	Y = np.linspace(0, 1.0*data_group.attrs.get('Aspect ratio')[0], num=data_group.attrs.get('Number of spatial points x')[0])
	print('Plotting graphs:')
	for key in velocity_x_group.keys():
		print(key)
		u = velocity_x_group[key]
		v = velocity_y_group[key]
		U = np.array(u)
		V = np.array(v)
		norma = np.sqrt(U**2, V**2)
		fig_name = 'velocity_stream_' + label + key + '.png'
		plt.pcolormesh(X,Y,norma,cmap='coolwarm')
		cbar=plt.colorbar()
		strm = plt.streamplot(X, Y, U, V, density=(.3,.5),color='k',linewidth=0.5)
		plt.title('Velocity streamplot')
		plt.xlabel('x')
		plt.ylabel('y')
		plt.savefig(fig_name,quality=100)
		plt.close("all")
f.close()
