#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 13:29:53 2020

@author: showard
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.optimize
#mpl.use('pgf')
mpl.rcParams['font.family'] = 'serif'

default_width = 4 #5.78853 # in inches
default_ratio = (np.sqrt(5.0) - 1.0) / 2.0 # golden mean

#mpl.rcParams.update({
#    "text.usetex": True,
#    "pgf.texsystem": "xelatex",
#    "pgf.rcfonts": False,
#    "font.family": "serif",
#    "font.serif": [],
#    "font.sans-serif": [],
#    "font.monospace": [],
#    "figure.figsize": [default_width, default_width * default_ratio],
#    "pgf.preamble": [
#        # put LaTeX preamble declarations here
#        r"\usepackage[utf8x]{inputenc}",
#        r"\usepackage[T1]{fontenc}",
#        # macros defined here will be available in plots, e.g.:
#        r"\newcommand{\vect}[1]{#1}",
#        # You can use dummy implementations, since your LaTeX document
#        # will render these properly, anyway.
#    ],
#})





plt.close("all")

freq0 = 7.1e6*2*np.pi
omega = np.logspace(-1,1,1000)*freq0
#omega = np.linspace(.25,3,1000)*freq0

# LC
Q= [10,20]
R=50

#
#fig1 = texfig.figure()
#ax11= plt.subplot(2,1,1)
#ax12 = plt.subplot(2,1,2)
fig1,(ax11,ax12) = plt.subplots(2,1)
for qi in Q:
    L=(qi)*R/freq0
    C=1/qi/R/freq0
    Z = 1/(1j*omega*C)+1j*L*omega+R
    H = R/(Z)
    ax11.loglog(omega/np.pi/2,abs(Z))
    valueatfreq0= H[omega.searchsorted(freq0, 'left')]
    ax12.semilogx(omega/np.pi/2,20*np.log10(np.abs(H)/np.abs(valueatfreq0)))

ax11.grid(which="both",axis="both")
ax11.axvline(freq0/2/np.pi, color='green') # cutoff frequency
ax11.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
ax11.set_ylabel(r"$|Z|\, (\Omega)$")

ax12.grid(which="both",axis="both")
ax12.axvline(freq0/2/np.pi, color='g') # cutoff frequency
ax12.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
ax12.axhline(-43, color='r') # cutoff frequency
ax12.axhline(-43+5.7, ls=':', color='r') # cutoff frequency
ax12.set_ylabel(r"$|H|$")
ax12.set_xlabel("Frequency (Hz)")

plt.savefig("LCbp.pdf")

# L-Impedance Match
Ri=10 #5x boost
Qboost=np.sqrt(R/Ri-1)
Cp=Qboost/freq0/R
#Lextra=1/freq0**2/Cp*(Qboost**2/(1+Qboost**2))
Lextra=Qboost*Ri/freq0

Q=(10, 20) #4 looks like 20, 2 looks like 10

fig2, (ax21,ax22) = plt.subplots(2,1)

for qi in Q:
    L=(qi)*Ri/freq0
    C=1/qi/Ri/freq0
    Zlast = 1/(1/R+1j*omega*Cp)
    Ztot = Zlast+1j*omega*Lextra+1/(1j*omega*C)+1j*L*omega
    H = Zlast/(Ztot)
    ax21.loglog(omega/np.pi/2,abs(Ztot))
    valueatfreq0= H[omega.searchsorted(freq0, 'left')]
    ax22.semilogx(omega/np.pi/2,20*np.log10(np.abs(H)/np.abs(valueatfreq0)))

ax21.grid(which="both",axis="both")
ax21.axvline(freq0/2/np.pi, color='green') # cutoff frequency
ax21.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
ax21.set_ylabel(r"$|Z|\, (\Omega)$")

ax22.grid(which="both",axis="both")
ax22.axvline(freq0/2/np.pi, color='g') # cutoff frequency
ax22.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
ax22.axhline(-43, color='r') # cutoff frequency
ax22.axhline(-43+5.7, ls=':', color='r') # cutoff frequency
ax22.set_ylabel(r"$|H|$")
ax22.set_xlabel("Frequency (Hz)")


plt.savefig("LCbpimp1.pdf")



#Tee network
Qtee = 5
Ri=25
R=50

fig3, (ax31,ax32) = plt.subplots(2,1)
fig4, (ax41,ax42) = plt.subplots(2,1)


Q=Qtee
L=(Q)*Ri/freq0
C=1/Q/Ri/freq0

Rv = scipy.optimize.fsolve(lambda x: np.sqrt(x/Ri-1)+np.sqrt(x/R-1)-Qtee, 100)
Q1 = np.sqrt(Rv/Ri-1)
Q2 = np.sqrt(Rv/R-1)
L2tee = Q1*Ri/freq0
L3tee = Q2*R/freq0
C2tee = Q1/Rv/freq0+Q2/Rv/freq0
print(L2tee,L3tee,C2tee)

Z1 = R+L3tee*1j*omega
Z12 = 1/(1j*C2tee*omega+1/Z1)
Z123 = Z12 + 1j*L2tee*omega+ 1j*L*omega+1/(1j*omega*C)

H1= R/(Z1)
H2 = Z12/(Z123)
H=H1*H2

ax31.loglog(omega/np.pi/2,abs(Z123))
valueatfreq0= H[omega.searchsorted(freq0, 'left')]
ax32.semilogx(omega/np.pi/2,20*np.log10(np.abs(H)/np.abs(valueatfreq0)))

L3teenotch=3/4*L3tee
C3tee = 1/(2*freq0)**2/L3teenotch
Z2ndfilt = 1/(1/(L3teenotch*1j*omega)+1j*omega*C3tee)
Z2nd1 = R+Z2ndfilt
Z2nd12= 1/(1j*C2tee*omega+1/Z2nd1)
Z2nd123 = Z2nd12 + 1j*L2tee*omega+ 1j*L*omega+1/(1j*omega*C)

H2nd1 = R/(Z2nd1)
H2nd2 = Z2nd12/Z2nd123
H2nd = H2nd1*H2nd2


ax41.loglog(omega/np.pi/2,abs(Z2nd123))
valueatfreq0= H2nd[omega.searchsorted(freq0, 'left')]
ax42.semilogx(omega/np.pi/2,20*np.log10(np.abs(H2nd)/np.abs(valueatfreq0)))

#ax31.grid(which="both",axis="both")
#ax31.axvline(freq0/2/np.pi, color='green') # cutoff frequency
#ax31.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
#ax31.set_ylabel(r"$|Z|\, (\Omega)$")
#
#
#ax32.grid(which="both",axis="both")
#ax32.axvline(freq0/2/np.pi, color='g') # cutoff frequency
#ax32.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
#ax32.axhline(-43, color='r') # cutoff frequency
#ax32.axhline(-43+5.7, ls=':', color='r') # cutoff frequency
#ax32.set_ylabel(r"$|H|$")
#ax32.set_xlabel("Frequency (Hz)")

#Pi
Qpi = 5
Ri=25
R=50

Q=Qpi
L=(Q)*Ri/freq0
C=1/Q/Ri/freq0

Rv = scipy.optimize.fsolve(lambda x: np.sqrt(Ri/x-1)+np.sqrt(R/x-1)-Qpi, 1)
Q1 = np.sqrt(Ri/Rv-1)
Q2 = np.sqrt(R/Rv-1)
C2pi = Q1/Ri/freq0
C3pi = Q2/R/freq0
L2pi = Rv/freq0*(Q1+Q2)


Z1 = 1/(1/R+C3pi*1j*omega)
Z12 = Z1+1j*L2pi*omega
Z123 = 1/(1/Z12+1j*C2pi*omega)
Z1234 = Z123+ 1j*L*omega+1/(1j*omega*C)


H1= Z1/(Z12)
H2 = Z123/(Z1234)
H=H1*H2

Z1234noLf=Z123+1/(1j*omega*100e-9)
Hnolof=H1*Z123/Z1234noLf
plt.figure("pitwithnoLf")
plt.subplot(211)
plt.semilogx(omega/np.pi/2,np.abs(Z1234noLf))
plt.subplot(212)
plt.semilogx(omega/np.pi/2,np.angle(Z1234noLf)*180/np.pi)


ax31.loglog(omega/np.pi/2,abs(Z1234))
valueatfreq0= H[omega.searchsorted(freq0, 'left')]
ax32.semilogx(omega/np.pi/2,20*np.log10(np.abs(H)/np.abs(valueatfreq0)))



ax32.grid(which="both",axis="both")
ax32.axvline(freq0/2/np.pi, color='g') # cutoff frequency
ax32.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
ax32.axhline(-43, color='r') # cutoff frequency
ax32.axhline(-43+5.7, ls=':', color='r') # cutoff frequency
ax32.set_ylabel(r"$|H|$")
ax32.set_xlabel("Frequency (Hz)")

ax31.grid(which="both",axis="both")
ax31.axvline(freq0/2/np.pi, color='green') # cutoff frequency
ax31.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
ax31.set_ylabel(r"$|Z|\, (\Omega)$")

plt.figure(fig3.number)
plt.savefig("LCbpimp2.pdf")
plt.figure(fig4.number)

L2pinotch=L2pi*3/4
C4pi = 1/(2*freq0)**2/L2pinotch
Z2nd1 = Z1
Z2ndfilt = 1/(1/(L2pinotch*1j*omega)+1j*omega*C4pi)
Z2nd12 = Z2nd1+Z2ndfilt
Z2nd123 = 1/(1j*C2pi*omega+1/Z2nd12)
Z2nd1234= Z2nd123 + 1j*L*omega+1/(1j*omega*C)

H2nd1 = Z2nd1/(Z2nd12)
H2nd2 = Z2nd123/Z2nd1234
H2nd = H2nd1*H2nd2






ax41.loglog(omega/np.pi/2,abs(Z2nd1234))
valueatfreq0= H2nd[omega.searchsorted(freq0, 'left')]
ax42.semilogx(omega/np.pi/2,20*np.log10(np.abs(H2nd)/np.abs(valueatfreq0)))




ax41.grid(which="both",axis="both")
ax41.axvline(freq0/2/np.pi, color='green') # cutoff frequency
ax41.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
ax41.set_ylabel(r"$|Z|\, (\Omega)$")


ax42.grid(which="both",axis="both")
ax42.axvline(freq0/2/np.pi, color='g') # cutoff frequency
ax42.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
ax42.axhline(-43, color='r') # cutoff frequency
ax42.axhline(-43+5.7, ls=':', color='r') # cutoff frequency
ax42.set_ylabel(r"$|H|$")
ax42.set_xlabel("Frequency (Hz)")

plt.savefig("LCbpimp3.pdf")


#QCX Filter
fig5, (ax51,ax52) = plt.subplots(2,1)
R=50

C29val=100e-9
C2728val=270e-12
C2526val=680e-12
L13val = 1.4e-6
L2val=1.7e-6

Z1 = 1/(1/R+C2728val*1j*omega)
Z2 = Z1+1j*L13val*omega
Z3 = 1/(1/Z2+1j*C2526val*omega)
Z4 = Z3+ 1j*L2val*omega
Z5 = 1/(1/Z4+1j*C2526val*omega)
Z6 = Z5+1j*L13val*omega
Z7 = 1/(1/Z6+C2728val*1j*omega)
Ztot = Z7+1/(1j*C29val*omega)

H1= Z1/(Z2)
H2 = Z3/(Z4)
H3 = Z5/Z6
H4=Z7/Ztot
H=H1*H2*H3*H4


ax51.loglog(omega/np.pi/2,abs(Ztot))
valueatfreq0= H[omega.searchsorted(freq0, 'left')]
ax52.semilogx(omega/np.pi/2,20*np.log10(np.abs(H)/np.abs(valueatfreq0)))



ax52.grid(which="both",axis="both")
ax52.axvline(freq0/2/np.pi, color='g') # cutoff frequency
ax52.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
ax52.axhline(-43, color='r') # cutoff frequency
ax52.axhline(-43+5.7, ls=':', color='r') # cutoff frequency
ax52.set_ylabel(r"$|H|$")
ax52.set_xlabel("Frequency (Hz)")

ax51.grid(which="both",axis="both")
ax51.axvline(freq0/2/np.pi, color='green') # cutoff frequency
ax51.axvline(freq0*2/2/np.pi, ls='--', color='g') # cutoff frequency
ax51.set_ylabel(r"$|Z|\, (\Omega)$")


plt.savefig("QCXfilt.pdf")

##Pi with no Lf