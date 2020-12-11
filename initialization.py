# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 22:26:59 2020

@author: az
"""
import numpy as np 
import matplotlib.pyplot as plt
from rusanov2D import rusanov2D
from rusanov1D import rusanov1D
from maccormack2D import maccormack2D
from maccormack1D import maccormack1D
from Riemann import Riemann

R=Riemann()
rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
case= 1 #case numbers are summarized above 
gam=1.4
xexact,rhoexact,uexact,pexact,eexact = R.exactsol(gam,max_time[case],case,"rusanov"+str(case))  

L = 1
Rc = L/10
xc = 2*L/3
yc = L/2
T0 = 298
p0 = 101300
gamma = 1.4
Mac = 0.3
dx = L/100
dy = L/100
N = int(L/dx+1)
xs = int(L/4)
Ns = int(N*xs)
vortex = True
shock = True
# SOD PROBLEM
rhoL = 1
uL = 0
pL = 1
aL = np.sqrt(gam*pL/rhoL)
rhoR = 0.125*rhoL
uR = 0
pR = 0.1*pL
aR = np.sqrt(gam*pR/rhoR)

rhoexact = np.roll(rhoexact,-250)
rhoexact[750:1001] = rhoR
uexact = np.roll(uexact,-250)
uexact[750:1001] = uR
pexact = np.roll(pexact,-250)
pexact[750:1001] = pR
eexact = np.roll(eexact,-250)
eR = pR/((gamma-1)*rhoR)
eexact[750:1001] = eR

CFL = 0.9
[pstar,ustar,rhostarL,astarL,SL,SHL,STL,rhostarR,astarR,SR,SHR,STR,rhoL,pL,aL,uL,rhoR,pR,aR,uR,gam,gamm1,gamp1] = R.find_star_state(gam, rhoL, pL, aL, uL, rhoR, pR, aR, uR)  
R = 287.058
t_final = (L/4)/SR
# First, establish some storage vectors 
x = np.linspace(0,L,N) # vectors for storing x and y coordinates
y = np.linspace(0,L,N) 
u = np.zeros((N,N))
v = np.zeros((N,N))
T = np.zeros((N,N))
rho = np.zeros((N,N))
p = np.zeros((N,N))
omegaex = np.zeros((N,N))

uinf = np.zeros((N,N))
vinf = np.zeros((N,N))
rhoinf = np.zeros((N,N))
pinf = np.zeros((N,N))

if shock == True:
    rhoinf[:,0:Ns] = rhoL
    rhoinf[:,Ns:N] = rhoR
    uinf[:,0:xs] = uL
    uinf[:,xs:N] = uR
    pinf = pR
    Tinf = pR/(rhoR*R)
    
if vortex == False:
    u = uinf
    v = vinf
    p = pinf
    rho = rhoinf
    e = p/((gamma-1)*rho)
    et = e+(u**2+v**2)/2
    H = e+p/rho+(u**2+v**2)/2

# This nested loop initializes the vortex
if vortex == True:
    Tc = Tinf/(1+((gamma-1)/2)*Mac**2) # Tc and ac are calculated since we need ac to use equation 4a and 4b
    ac = np.sqrt(gamma*R*Tc)
    for j in range(0,N):
        for i in range(0,N):
            rstar = np.sqrt((x[i]-xc)**2+(y[j]-yc)**2)/Rc #Calculate rstar, ystar, and xstar here since they are used repeatedly
            ystar = (y[j]-yc)/Rc
            xstar = (x[i]-xc)/Rc
            rstarexp = np.exp((1-(rstar)**2)/2) #This is the exponential function and everything inside it
            u[j,i] = uinf[j,i]-Mac*ac*ystar*rstarexp #Calculate u and v using equation 4 at every node - note that the indexing is j,i since rows in the array correspond to the y-axis
            v[j,i] = vinf[j,i]+Mac*ac*xstar*rstarexp 
            T[j,i] = Tinf*(1-(((gamma-1)/2)*Mac**2/(1+((gamma-1)/2)*Mac**2))*rstarexp*rstar**2) # Calculate temperature using eq. 6
            omegaex[j,i] = (Mac*ac/Rc)*rstarexp*(2-xstar**2-ystar**2) 
            
    p = pinf*(T/Tinf)**((gamma)/(gamma-1)) # Can now calculate p, rho, e based on the other values we have
    rho = p/(R*T)
    e = p/((gamma-1)*rho)
    et = e+(u**2+v**2)/2
    H = e+p/rho+(u**2+v**2)/2
# Initialization complete

origin = 'lower'
fig1, ax1 = plt.subplots(2,2,figsize=(8,8),constrained_layout=True,sharex=True)
fig1.suptitle('test',fontsize=10, weight='bold')
CS1 = ax1[0,0].contourf(x, y, u, 10, origin=origin)
CS2 = ax1[1,0].contourf(x, y, T, 10 ,origin=origin)
CS3 = ax1[0,1].contourf(x, y, rho, 10, origin=origin)
CS4 = ax1[1,1].contourf(x, y, p, 10 ,origin=origin)
ax1[0,0].set_title('velocity',fontsize=10)
ax1[0,1].set_title('density',fontsize=10)
ax1[1,0].set_title('temperature',fontsize=10)
ax1[1,1].set_title('pressure',fontsize=10)
cbar1 = fig1.colorbar(CS1,ax=ax1[0,0],location='bottom',pad=0)
cbar2 = fig1.colorbar(CS2,ax=ax1[1,0],location='bottom',pad=0)
cbar3 = fig1.colorbar(CS3,ax=ax1[0,1],location='bottom',pad=0)
cbar4 = fig1.colorbar(CS4,ax=ax1[1,1],location='bottom',pad=0)

u1, v1, rho1, p1, T1, e1, omega1 = rusanov2D(N,L,t_final,CFL,u,v,rho,p,gamma,R) 
rho2, u2, p2, e2, x2 = rusanov1D("Case 1",t_final,CFL,dx,gamma)
u3, v3, rho3, p3, T3, e3, omega3 = maccormack2D(N,L,t_final,CFL,u,v,rho,p,gamma,R) 
rho4, u4, p4, e4, x4 = maccormack1D("Case 1",t_final,CFL,dx,gamma)

fig2, ax1 = plt.subplots(2,2,figsize=(8,8),constrained_layout=True,sharex=True)
fig2.suptitle('test',fontsize=10, weight='bold')
CS1 = ax1[0,0].contourf(x, y, u1, 10, origin=origin)
CS2 = ax1[1,0].contourf(x, y, e1, 10 ,origin=origin)
CS3 = ax1[0,1].contourf(x, y, rho1, 10, origin=origin)
CS4 = ax1[1,1].contourf(x, y, p1, 10 ,origin=origin)
ax1[0,0].set_title('velocity',fontsize=10)
ax1[0,1].set_title('density',fontsize=10)
ax1[1,0].set_title('energy',fontsize=10)
ax1[1,1].set_title('pressure',fontsize=10)
cbar1 = fig2.colorbar(CS1,ax=ax1[0,0],location='bottom',pad=0)
cbar2 = fig2.colorbar(CS2,ax=ax1[1,0],location='bottom',pad=0)
cbar3 = fig2.colorbar(CS3,ax=ax1[0,1],location='bottom',pad=0)
cbar4 = fig2.colorbar(CS4,ax=ax1[1,1],location='bottom',pad=0)

fig3, ax1 = plt.subplots(2,2,figsize=(8,4),constrained_layout=True,sharex=True)
fig3.suptitle('rusanov test',fontsize=10, weight='bold')
ax1[0,0].plot(xexact,uexact,'k--',label = 'Exact')
ax1[0,0].plot(x, u1[int(np.floor(N/2)),:],label='2D')
ax1[0,0].plot(x2, u2,label='1D')
ax1[1,0].plot(xexact,eexact,'k--',label = 'Exact')
ax1[1,0].plot(x, e1[int(np.floor(N/2)),:],label='2D')
ax1[1,0].plot(x2, e2,label='1D')
ax1[0,1].plot(xexact,rhoexact,'k--',label = 'Exact')
ax1[0,1].plot(x, rho1[int(np.floor(N/2)),:],label='2D')
ax1[0,1].plot(x2, rho2,label='1D')
ax1[1,1].plot(xexact,pexact,'k--',label = 'Exact')
ax1[1,1].plot(x, p1[int(np.floor(N/2)),:],label='2D')
ax1[1,1].plot(x2, p2,label='1D')
ax1[0,0].set_title('velocity',fontsize=10)
ax1[0,1].set_title('density',fontsize=10)
ax1[1,0].set_title('energy',fontsize=10)
ax1[1,1].set_title('pressure',fontsize=10)
ax1[0,0].legend(loc='best', bbox_to_anchor=(0.5,-0.3))

fig5, ax1 = plt.subplots(2,2,figsize=(8,4),constrained_layout=True,sharex=True)
fig5.suptitle('maccormack test',fontsize=10, weight='bold')
ax1[0,0].plot(xexact,uexact,'k--',label = 'Exact')
ax1[0,0].plot(x, u3[int(np.floor(N/2)),:],label='2D')
ax1[0,0].plot(x4, u4,label='1D')
ax1[1,0].plot(xexact,eexact,'k--',label = 'Exact')
ax1[1,0].plot(x, e3[int(np.floor(N/2)),:],label='2D')
ax1[1,0].plot(x4, e4,label='1D')
ax1[0,1].plot(xexact,rhoexact,'k--',label = 'Exact')
ax1[0,1].plot(x, rho3[int(np.floor(N/2)),:],label='2D')
ax1[0,1].plot(x4, rho4,label='1D')
ax1[1,1].plot(xexact,pexact,'k--',label = 'Exact')
ax1[1,1].plot(x, p3[int(np.floor(N/2)),:],label='2D')
ax1[1,1].plot(x4, p4,label='1D')
ax1[0,0].set_title('velocity',fontsize=10)
ax1[0,1].set_title('density',fontsize=10)
ax1[1,0].set_title('energy',fontsize=10)
ax1[1,1].set_title('pressure',fontsize=10)
ax1[0,0].legend(loc='best', bbox_to_anchor=(0.5,-0.3))

fig4, ax1 = plt.subplots(2,2,figsize=(8,4),constrained_layout=True,sharex=True)
fig4.suptitle('test',fontsize=10, weight='bold')
ax1[0,0].plot(x, omega1[int(np.floor(N/2)),:],label='rusanov')
ax1[0,0].plot(x, omega3[int(np.floor(N/2)),:],label='maccormack')
ax1[0,0].plot(x, omegaex[int(np.floor(N/2)),:],label='exact')
ax1[0,0].set_title('velocity',fontsize=10)
ax1[0,1].set_title('density',fontsize=10)
ax1[1,0].set_title('energy',fontsize=10)
ax1[1,1].set_title('pressure',fontsize=10)
ax1[0,0].legend(loc='best', bbox_to_anchor=(0.5,-0.3))

origin = 'lower'
fig6, ax1 = plt.subplots(2,2,figsize=(8,8),constrained_layout=True,sharex=True)
fig6.suptitle('test',fontsize=10, weight='bold')
CS1 = ax1[0,0].contourf(x, y, u1, 10, origin=origin)
CS2 = ax1[1,0].contourf(x, y, u3, 10 ,origin=origin)
CS3 = ax1[0,1].contourf(x, y, omega1, 10, origin=origin)
CS4 = ax1[1,1].contourf(x, y, omega3, 10 ,origin=origin)
ax1[0,0].set_title('velocity',fontsize=10)
ax1[0,1].set_title('density',fontsize=10)
ax1[1,0].set_title('temperature',fontsize=10)
ax1[1,1].set_title('pressure',fontsize=10)
cbar1 = fig6.colorbar(CS1,ax=ax1[0,0],location='bottom',pad=0)
cbar2 = fig6.colorbar(CS2,ax=ax1[1,0],location='bottom',pad=0)
cbar3 = fig6.colorbar(CS3,ax=ax1[0,1],location='bottom',pad=0)
cbar4 = fig6.colorbar(CS4,ax=ax1[1,1],location='bottom',pad=0)