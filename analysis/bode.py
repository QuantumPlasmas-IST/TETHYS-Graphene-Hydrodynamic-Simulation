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
with open('delay_data.txt') as data:
    data_reader = csv.reader(data, delimiter=';')
    for row in data_reader:
        sum+=1
        
Vdelta=np.zeros(sum)
Vamp=np.zeros(sum)
Vfreq=np.zeros(sum)
sum=0

with open('delay_data.txt') as data:
    data_reader = csv.reader(data, delimiter=';')
    for row in data_reader:
        Vdelta[sum]=row[0]
        Vamp[sum]=row[1]
        Vfreq[sum]=row[2]
        sum+=1



#Draw freq/amp vs delay

ax = plt.gca()
plt.plot(Vdelta,Vamp/2)
ax.set_xlabel("Delay [L/v0]")
ax.set_ylabel("Gain")
plt.savefig("../img/delay_new.png",dpi=250)
plt.clf()

ax = plt.gca()
plt.plot(Vdelta,Vfreq)
ax.set_xlabel("Delay [L/v0]")
ax.set_ylabel("Freq")
plt.ylim(0.0, 10)
plt.savefig("../img/freq_new.png",dpi=250)
plt.clf()

########################################################
sum=0
with open('delay_data_noDS.txt') as data:
    data_reader = csv.reader(data, delimiter=';')
    for row in data_reader:
        sum+=1
        
Vdelta_noDS=np.zeros(sum)
Vamp_noDS=np.zeros(sum)
Vfreq_noDS=np.zeros(sum)
sum=0

with open('delay_data_noDS.txt') as data:
    data_reader = csv.reader(data, delimiter=';')
    for row in data_reader:
        Vdelta_noDS[sum]=row[0]
        Vamp_noDS[sum]=row[1]
        Vfreq_noDS[sum]=row[2]
        sum+=1



#Draw freq/amp vs delay

ax = plt.gca()
plt.plot(Vdelta_noDS,Vamp_noDS/2)
ax.set_xlabel("Delay [L/v0]")
ax.set_ylabel("Gain")
plt.savefig("../img/delay_new_noDS.png",dpi=250)
plt.clf()

ax = plt.gca()
plt.plot(Vdelta_noDS,Vfreq_noDS)
ax.set_xlabel("Delay [L/v0]")
ax.set_ylabel("Freq")
plt.ylim(0.0, 10)
plt.savefig("../img/freq_new_noDS.png",dpi=250)
plt.clf()
#######################################################

## ANALYTICAL

S=20
vf=15
L=1
W=1
C=0.01
I= 0.+1.j
eps=0.2

theta=np.sqrt(1+16*S*S+8*vf*vf)/(2*(1-2*S*S-vf*vf))*L
phi=3*L/(2*(1-2*S*S-vf*vf))

T=np.zeros([2,2],dtype=complex)

G=np.zeros([2,2],dtype=complex)

H=np.zeros([2,2],dtype=complex)


s1=np.zeros(1000)
s2=np.zeros(1000)

#BODE PLOT

freq=np.linspace(0,999,1000)

for w in range(0,999):
    T[0,0]=np.exp(I*phi*w)/theta*(theta*np.cos(theta*w) + I*phi*np.sin(theta*w))
    T[0,1]=np.exp(I*phi*w)/theta*I*(theta*theta - phi*phi)/(W*C)*np.sin(theta*w)
    T[1,0]=np.exp(I*phi*w)/theta*I*W*C*np.sin(theta*w)
    T[1,1]=np.exp(I*phi*w)/theta*(theta*np.cos(theta*w) - I*phi*np.sin(theta*w))


    H[0,0]=eps
    H[0,1]=0
    H[1,0]=0
    H[1,1]=0

    G[:,:]=np.matmul(np.linalg.inv(np.identity(2)-np.matmul(T,H)),T)

    s1[w]=(np.linalg.svd(G)[1])[0]
    s2[w]=(np.linalg.svd(G)[1])[1]

ax = plt.gca()
plt.plot(freq*2*cmath.pi,s1,linewidth=1)
ax.set_xlabel("Frequency [v0/L]")
ax.set_ylabel("Gain")
plt.xlim(0.0, 100)
#plt.savefig("../img/s1.png",dpi=250)
#plt.clf()

#ax = plt.gca()
plt.plot(freq*2*cmath.pi,s2,linewidth=1)
#ax.set_xlabel("Frequency [v0/L]")
#ax.set_ylabel("Gain")
plt.savefig("../img/svd.png",dpi=250)
plt.clf()

#NYQUIST PLOT

l1=np.zeros(2000+100, dtype=complex)
l2=np.zeros(2000+100, dtype=complex)

for w in range(0,1999+100):
    T[0,0]=np.exp(I*phi*w)/theta*(theta*np.cos(theta*w) + I*phi*np.sin(theta*w))
    T[0,1]=np.exp(I*phi*w)/theta*I*(theta*theta - phi*phi)/(W*C)*np.sin(theta*w)
    T[1,0]=np.exp(I*phi*w)/theta*I*W*C*np.sin(theta*w)
    T[1,1]=np.exp(I*phi*w)/theta*(theta*np.cos(theta*w) - I*phi*np.sin(theta*w))

    H[0,0]=eps
    H[0,1]=0
    H[1,0]=0
    H[1,1]=0

    G[:,:]=np.matmul(T,H)

    l1[w]=(np.linalg.eigvals(G))[0]
    l2[w]=(np.linalg.eigvals(G))[1]

ax = plt.gca()
plt.plot(np.real(l1),np.imag(l1),linewidth=1)
ax.set_xlabel("Real")
ax.set_ylabel("Imag")
plt.savefig("../img/l1.png",dpi=250)
plt.clf()

ax = plt.gca()
plt.plot(np.real(l2),np.imag(l2),linewidth=1)
ax.set_xlabel("Real")
ax.set_ylabel("Imag")
plt.savefig("../img/l2.png",dpi=250)
plt.clf()

#delay bode plot

w=5*2*cmath.pi

count=0

for d in np.linspace(0,1,1000):
    T[0,0]=np.exp(I*phi*w)/theta*(theta*np.cos(theta*w) + I*phi*np.sin(theta*w))
    T[0,1]=np.exp(I*phi*w)/theta*I*(theta*theta - phi*phi)/(W*C)*np.sin(theta*w)
    T[1,0]=np.exp(I*phi*w)/theta*I*W*C*np.sin(theta*w)
    T[1,1]=np.exp(I*phi*w)/theta*(theta*np.cos(theta*w) - I*phi*np.sin(theta*w))

    H[0,0]=eps*np.exp(-I*w*d)
    H[0,1]=0
    H[1,0]=0
    H[1,1]=0

    G[:,:]=np.matmul(np.linalg.inv(np.identity(2)-np.matmul(T,H)),T)

    s1[count]=(np.linalg.svd(G)[1])[0]
    s2[count]=(np.linalg.svd(G)[1])[1]
    count+=1

ax = plt.gca()
plt.plot(np.linspace(0,1,1000),s1,linewidth=1)
ax.set_xlabel("Delay [L/v0]")
ax.set_ylabel("Gain")
plt.xlim(0.0, 1)
plt.plot(np.linspace(0,1,1000),s2,linewidth=1)
plt.savefig("../img/delay.png",dpi=250)
plt.clf()

#delay nyquist plot

#overlap

ax = plt.gca()
plt.plot(np.linspace(0,1,1000),(s1-min(s1))/(max(s1)-min(s1)),linewidth=1)
plt.plot(Vdelta,(Vamp-min(Vamp))/(max(Vamp)-min(Vamp)))
ax.set_xlabel("Delay [L/v0]")
ax.set_ylabel("Gain")
plt.xlim(0.0, 1.1)
plt.savefig("../img/overlap.png",dpi=250)
plt.clf()

ax = plt.gca()
plt.plot(np.linspace(0,1,1000),(s1-min(s1))/(max(s1)-min(s1)),linewidth=1)
plt.plot(Vdelta_noDS,(Vamp_noDS-min(Vamp_noDS))/(max(Vamp_noDS)-min(Vamp_noDS)))
ax.set_xlabel("Delay [L/v0]")
ax.set_ylabel("Gain")
plt.xlim(0.0, 1.1)
plt.savefig("../img/overlap_noDS.png",dpi=250)
plt.clf()