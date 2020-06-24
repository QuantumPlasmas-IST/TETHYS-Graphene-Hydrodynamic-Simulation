
# import required libraries
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt


# Read H5 file
f = h5.File("OPEN_BC_RichtmyerCLASSTEST_2D_Graphene_DS.h5", "r")

density_group = f['Data/Density']
velocity_x_group = f['Data/VelocityX']
velocity_y_group = f['Data/VelocityY']

# for key in density_group.keys():
	# data_density=density_group[key]
	# fig_name = key + '.jpeg'
	# plt.pcolormesh(data_density,cmap='Spectral')
	# plt.clim(0.3, 1.9);
	# cbar=plt.colorbar()
	# cbar.set_label('n/n0')
	# plt.xlabel('x')
	# plt.ylabel('y')
	# plt.savefig(fig_name,quality=60)
	# plt.close("all")

for key in velocity_x_group.keys():
	u = velocity_x_group[key]
	v = velocity_y_group[key]
	U = np.array(u)
	V = np.array(v)
	Y = np.linspace(0, 1, num=201)
	X = np.linspace(0, 1, num=201)
	norma = np.sqrt(U**2, V**2)
	#strm = plt.streamplot(X, Y, U, V,color=norma,cmap='Blues', density = 1)
	strm = plt.streamplot(X, Y, U, V, density=1)
	#cbar = plt.colorbar(strm.lines)
	#cbar.set_label('|v|')
	fig_name = 'vel_' + key + '.jpeg'
	plt.savefig(fig_name, quality=60)
	plt.close("all")

#total = np.loadtxt('OpenBC_slice_test.dat')
#t=total[:,0]
#den=total[:,2]
#plt.plot(t,den)
#plt.show()
#from scipy import signal
#fs=1.0/t[1]
#sinal = den-np.average(den)
#Fspec, Tspec, Sxx = signal.spectrogram(sinal, fs)
#plt.pcolormesh(Tspec, Fspec, Sxx)

f.close()
