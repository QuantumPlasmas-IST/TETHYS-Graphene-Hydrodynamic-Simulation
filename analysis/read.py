import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmath
import scipy as scp
import scipy.io as sio
import csv
import scipy.signal as signal
import scipy.fft as fourier
import pybaselines


sum=0
with open('../build/preview_1D_S=20.00vF=15.00vis=0.00l=0.00.dat') as data:
    data_reader = csv.reader(data, delimiter='\t')
    for row in data_reader:
        sum+=1
        
time=np.zeros(sum)
n0 = np.zeros(sum)
v0 = np.zeros(sum)
nL = np.zeros(sum)
vL = np.zeros(sum)
sum=0

with open('../build/preview_1D_S=20.00vF=15.00vis=0.00l=0.00.dat') as data:
    data_reader = csv.reader(data, delimiter='\t')
    for row in data_reader:
        time[sum]=row[0]
        nL[sum]=row[1]
        vL[sum]=row[2]
        n0[sum]=row[3]
        v0[sum]=row[4]
        sum+=1

#Draw density and velocity

ax = plt.gca()
plt.plot(time,n0,linewidth=.5)
ax.set_xlabel("Time [L/v0]")
ax.set_ylabel("Density [n0]")
plt.savefig("../img/n0.png",dpi=250)
plt.clf()

ax = plt.gca()
plt.plot(time,nL,linewidth=.5)
ax.set_xlabel("Time [L/v0]")
ax.set_ylabel("Density [n0]")
plt.savefig("../img/nL.png",dpi=250)
plt.clf()

ax = plt.gca()
plt.plot(time,v0,linewidth=.5)
ax.set_xlabel("Time [L/v0]")
ax.set_ylabel("Velocity [v0]")
plt.savefig("../img/v0.png",dpi=250)
plt.clf()

ax = plt.gca()
plt.plot(time,vL,linewidth=.5)
ax.set_xlabel("Time [L/v0]")
ax.set_ylabel("Velocity [v0]")
plt.savefig("../img/vL.png",dpi=250)
plt.clf()


#frequency analysis

dt = time[1]-time[0]
fs = int(1.0/dt)
nfft = 2048*8
tmax = time[-1]

print(fs)

fig, (ax1,ax2) = plt.subplots(nrows=2)
ax1.plot(time, nL,linewidth=.5)
ax1.set_xlim(0.0, tmax)
ax1.set_ylabel(r'$n(L)/n_0$')
ax1.set_xlabel(r'$t [L/v_0]$')
#Pxx, freqs, bins, im = ax2.specgram(nL-np.mean(nL), Fs=fs, NFFT=nfft, noverlap=nfft/2, mode='psd', scale='dB', detrend='linear',cmap='viridis')
#ax2.set_xlim(0.0, tmax)
#ax2.set_ylim(0.0, 100)
#ax2.set_ylabel(r'$f$')
spectrum, freqs, im_spectrum = ax2.magnitude_spectrum(nL-np.mean(nL),Fs=fs,scale='linear')
plt.xlim(0.0, 100)
fig_name = '../img/spectrogram_nL.png'
plt.savefig(fig_name, dpi=250)
plt.close("all")

print("FREQ =",freqs[np.argmax(spectrum)])

print("AMP =",np.max(nL))

print("return 0")

delta=0.26

f= open("delay_data_noDS.txt","a")
f.write("%f;%f;%f\n" % (delta,np.max(nL),freqs[np.argmax(spectrum)]))
#f.write(delta,";",np.max(nL),";",freqs[np.argmax(spectrum)])
f.close()

