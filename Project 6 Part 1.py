# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 22:26:59 2020

@author: az
"""
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from rusanov2D import rusanov2D
from rusanov1D import rusanov1D
from maccormack2D import maccormack2D
from maccormack1D import maccormack1D
from Riemann import Riemann

# Computational domain parameters and flow parameters
L = 1
gamma = 1.4
dx = L/25
dy = L/25
N = int(L/dx+1)

xs = L/4
Ns = int(N*xs)

R=Riemann()
rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
case= 1 #case numbers are summarized above 
gamma=1.4
xexact,rhoexact,uexact,pexact,eexact = R.exactsol(gamma,xs,max_time[case],case,"rusanov"+str(case)) 

# SOD PROBLEM
rhoL = 1
uL = 0
pL = 1
aL = np.sqrt(gamma*pL/rhoL)
rhoR = 0.125*rhoL
uR = 0
pR = 0.1*pL
aR = np.sqrt(gamma*pR/rhoR)
aL = np.sqrt(gamma*pL/rhoL)
aR = np.sqrt(gamma*pR/rhoR)

[pstar,ustar,rhostarL,astarL,SL,SHL,STL,rhostarR,astarR,SR,SHR,STR,rhoL,pL,aL,uL,rhoR,pR,aR,uR,gamma,gamm1,gamp1] = R.find_star_state(gamma, rhoL, pL, aL, uL, rhoR, pR, aR, uR)  

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

rho[:,0:Ns] = rhoL
rho[:,Ns:N] = rhoR
u[:,0:Ns] = uL
u[:,Ns:N] = uR
p[:,0:Ns] = pL
p[:,Ns:N] = pR

T = p/(R*rho)
e = p/((gamma-1)*rho)
et = e+(u**2+v**2)/2
H = e+p/rho+(u**2+v**2)/2

CFL = 1.0
u4, v4, rho4, p4, T4, e4, omega4,timevect4,ens4,Sdil4,Stor4,dilatation_fin4,baroclinic_fin4 = rusanov2D(N,L,t_final,CFL,u,v,rho,p,gamma,R) 
rho2, u2, p2, e2, x2, = rusanov1D("Case 1",t_final,CFL,dx,gamma)
u3, v3, rho3, p3, T3, e3, omega3,timevect3,ens3,Sdil3,Stor3,dilatation_fin3,baroclinic_fin3 = maccormack2D(N,L,t_final,CFL,u,v,rho,p,gamma,R) 
rho1, u1, p1, e1, x1 = maccormack1D("Case 1",t_final,CFL,dx,gamma)

yL2 = int(N/2)

fig1, axs1 = plt.subplots(4,figsize=(6.5,8))
fig1.suptitle("Figure 1: Sod problem, Δx = Δy = L/25, $x_s$ = L/4, t = (L/4)/S, $c_{max}$ = 1.0",fontsize=10)
axs1[0].plot(xexact, rhoexact,"k--")
axs1[1].plot(xexact, uexact,"k--")
axs1[2].plot(xexact, pexact,"k--")
axs1[3].plot(xexact, eexact,"k--",label="Analytical")

axs1[0].plot(x1, rho1)
axs1[0].grid()
axs1[1].plot(x1, u1)
axs1[1].grid()
axs1[2].plot(x1, p1)
axs1[2].grid()
axs1[3].plot(x1, e1,label="1D MacCormack")
axs1[3].grid()

axs1[0].plot(x, rho3[yL2,:],'-.')
axs1[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs1[0].set_xlim([0, 1])
axs1[1].plot(x, u3[yL2,:],'-.')
axs1[1].set(ylabel="$u$ $(m/s)$")
axs1[1].set_xlim([0, 1])
axs1[1].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs1[2].plot(x, p3[yL2,:],'-.')
axs1[2].set(ylabel="$p$ $(Pa)$")
axs1[2].set_xlim([0, 1])
axs1[2].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs1[3].plot(x, e3[yL2,:],'-.',label="2D MacCormack")
axs1[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs1[3].set_xlim([0, 1])
axs1[3].yaxis.set_major_locator(ticker.MultipleLocator(0.2))

axs1[0].plot(x, rho2)
axs1[0].grid()
axs1[1].plot(x, u2)
axs1[1].grid()
axs1[2].plot(x, p2)
axs1[2].grid()
axs1[3].plot(x, e2,label="1D Rusanov")
axs1[3].grid()

axs1[0].plot(x, rho4[yL2,:],'-.')
axs1[0].grid()
axs1[1].plot(x, u4[yL2,:],'-.')
axs1[1].grid()
axs1[2].plot(x, p4[yL2,:],'-.')
axs1[2].grid()
axs1[3].plot(x, e4[yL2,:],'-.',label="2D Rusanov")
axs1[3].grid()
axs1[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=3)

fig1.subplots_adjust(top=0.95)
fig1.savefig('HW6 Figure 3.png', dpi=300)