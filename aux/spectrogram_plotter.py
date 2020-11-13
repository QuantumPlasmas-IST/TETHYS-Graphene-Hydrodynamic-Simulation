# import required libraries
import glob
import numpy as np
import matplotlib.pyplot as plt

# Read file
for filename in glob.glob("preview_*.dat"):
	data = np.loadtxt(filename)
	label = filename.replace('preview_','').replace('.dat','')
	time = data[:, 0]
	denL = data[:, 1]
	velL = data[:, 2]
	den0 = data[:, 3]
	vel0 = data[:, 4]
	dt = time[1]-time[0]
	fs = int(1.0/dt)
	nfft = 512
	tmax = time[-1]
	fig, (ax1,ax2) = plt.subplots(nrows=2)
	ax1.plot(time, denL)
	ax1.set_xlim(0.0, tmax)
	ax1.set_ylabel(r'$n(L)/n_0$')
	Pxx, freqs, bins, im = ax2.specgram(denL-np.mean(denL), Fs=fs, NFFT=nfft, noverlap=nfft/2, mode='psd', scale='dB', detrend='linear',cmap='coolwarm')
	ax2.set_xlim(0.0, tmax)
	ax2.set_xlabel(r'$t [L/v_0]$')
	ax2.set_ylabel(r'$f$')
	fig_name = 'spectrogram_' + label + '.png'
	plt.savefig(fig_name, quality=100)
	plt.close("all")
	#spectrum, freqs, im_spectrum = plt.magnitude_spectrum(denL-np.mean(denL),Fs=fs,scale='linear')
	#Pxx_psd, freqs_psd = plt.psd(denL-np.mean(denL),NFFT=nfft,noverlap=nfft/2,Fs=fs,scale_by_freq=False)
	#fig_name = 'psd_' + label + '.png'
	#plt.savefig(fig_name, quality=100)
	#plt.close("all")
