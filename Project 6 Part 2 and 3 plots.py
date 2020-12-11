# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 22:26:59 2020

@author: az
"""
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import string
from project6part2 import project6part2

L = 1
dx1 = L/25
dx2 = L/50
dx3 = L/100
time = 2
N1 = int(L/dx1+1)
N2 = int(L/dx2+1)
N3 = int(L/dx3+1)

x1,y1,omegaex,omegam1,omegar1,timevectm1,timevectr1,ensm1,ensr1,Sdilm1,Sdilr1,Storm1,Storr1,dilatation_finm1,dilatation_finr1,baroclinic_finm1,baroclinic_finr1 = project6part2(dx1,time)
x2,y2,omegaex,omegam2,omegar2,timevectm2,timevectr2,ensm2,ensr2,Sdilm2,Sdilr2,Storm2,Storr2,dilatation_finm2,dilatation_finr2,baroclinic_finm2,baroclinic_finr2 = project6part2(dx2,time)
x3,y3,omegaex,omegam3,omegar3,timevectm3,timevectr3,ensm3,ensr3,Sdilm3,Sdilr3,Storm3,Storr3,dilatation_finm3,dilatation_finr3,baroclinic_finm3,baroclinic_finr3 = project6part2(dx3,time)

yh1 = int(N1/2)
yh2 = int(N2/2)
yh3 = int(N3/2)

#fig1, axs1 = plt.subplots(4,figsize=(6,8))
#fig1.suptitle("Figure 3: Sod problem, Δx = Δy = L/100, $x_s$ = L/4, t = (L/4)/S, $c_{max}$ = 1.0",fontsize=10)

origin = 'lower'
fig1, ax1 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig1.suptitle('Figure 4: Comparison of $\omega$ contours, Δx = Δy = L/100, $x_s$ = L/4, t = (L/4)/S',fontsize=10)
ax1[0].set_title('MacCormack method',fontsize=10)
CS1 = ax1[0].contour(x3, y3, omegaex, 10, origin=origin ,linestyles='dashed')
CS2 = ax1[0].contour(x3, y3, omegam3, 10 ,origin=origin)
ax1[0].set_ylabel('y (m)')
CS3 = ax1[1].contour(x3, y3, omegaex, 10, origin=origin ,linestyles='dashed')
CS4 = ax1[1].contour(x3, y3, omegar3, 10 ,origin=origin)
ax1[1].set_title('Rusanov method',fontsize=10)
ax1[0].clabel(CS1, CS1.levels, inline=True, fontsize=10)
ax1[0].clabel(CS2, CS2.levels, inline=True, fontsize=10)
ax1[1].clabel(CS3, CS3.levels, inline=True, fontsize=10)
ax1[1].clabel(CS4, CS4.levels, inline=True, fontsize=10)
ax1[1].set_xlabel('x (m)')
ax1[1].set_ylabel('y (m)')
h1,_ = CS1.legend_elements()
h2,_ = CS2.legend_elements()
h3,_ = CS3.legend_elements()
h4,_ = CS4.legend_elements()
ax1[0].legend([h1[0], h2[0]], ['Initial', 'MacCormack'],loc = 'best')
ax1[1].legend([h3[0], h4[0]], ['Initial', 'Rusanov'],loc = 'best')
for i, ax in enumerate(ax1):
    ax.text(1.025, 0.5, string.ascii_uppercase[i], transform=ax.transAxes,size=10)
    ax.text(1.0175, 0.5, '(  )', transform=ax.transAxes,size=10)
  
fig2, ax2 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig2.suptitle('Figure 7: Comparison of $\omega$ profiles at $y/L = 1/2$, $x_s$ = L/4, t = (L/2)/S',fontsize=10)
ax2[0].plot(x3, omegaex[yh3,:],"--k",label='Initial')
ax2[0].plot(x1, omegam1[yh1,:],label='Δx = Δy = L/25')
ax2[0].plot(x2, omegam2[yh2,:],label='Δx = Δy = L/50')
ax2[0].plot(x3, omegam3[yh3,:],label='Δx = Δy = L/100')
ax2[1].plot(x3, omegaex[yh3,:],"--k",label='Initial')
ax2[1].plot(x1, omegar1[yh1,:],label='Δx = Δy = L/25')
ax2[1].plot(x2, omegar2[yh2,:],label='Δx = Δy = L/50')
ax2[1].plot(x3, omegar3[yh3,:],label='Δx = Δy = L/100')
ax2[0].set_title('MacCormack method',fontsize=10)
ax2[0].set_ylabel('$\omega$ (1/s)')
ax2[1].set_title('Rusanov method',fontsize=10)
ax2[1].set_xlabel('x (m)')
ax2[1].set_ylabel('$\omega$ (1/s)')
ax2[0].grid(True)
ax2[1].grid(True)
ax2[0].set_xlim([0, 1])
ax2[1].set_xlim([0, 1])
ax2[0].legend(loc='best')
ax2[1].legend(loc='best')
for i, ax in enumerate(ax2):
    ax.text(1.025, 0.5, string.ascii_uppercase[i], transform=ax.transAxes,size=10)
    ax.text(1.0175, 0.5, '(  )', transform=ax.transAxes,size=10)
    
fig3, ax2 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig3.suptitle('Figure 8: Time history of $\epsilon$, $x_s$ = L/4, t = (L/2)/S',fontsize=10)
ax2[0].plot(timevectm1, ensm1,label='Δx = Δy = L/25')
ax2[0].plot(timevectm2, ensm2,label='Δx = Δy = L/50')
ax2[0].plot(timevectm3, ensm3,label='Δx = Δy = L/100')
ax2[1].plot(timevectr1, ensr1,label='Δx = Δy = L/25')
ax2[1].plot(timevectr2, ensr2,label='Δx = Δy = L/50')
ax2[1].plot(timevectr3, ensr3,label='Δx = Δy = L/100')
ax2[0].set_title('MacCormack method',fontsize=10)
ax2[0].set_ylabel('$\epsilon$ ($m^{2}/s^{2}$)')
ax2[1].set_title('Rusanov method',fontsize=10)
ax2[1].set_xlabel('t (s)')
ax2[1].set_ylabel('$\epsilon$ ($m^{2}/s^{2}$)')
ax2[0].grid(True)
ax2[1].grid(True)
ax2[1].legend(loc='best')
ax2[0].set_xlim([0, 0.29])
ax2[1].set_xlim([0, 0.29])
for i, ax in enumerate(ax2):
    ax.text(1.025, 0.5, string.ascii_uppercase[i], transform=ax.transAxes,size=10)
    ax.text(1.0175, 0.5, '(  )', transform=ax.transAxes,size=10)
ax2[0].legend(loc='best')

fig4, ax2 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig4.suptitle('Figure 9: Time history of $\dot S_{dil}$, $x_s$ = L/4, t = (L/2)/S',fontsize=10)
ax2[0].plot(timevectm1, Sdilm1,label='Δx = Δy = L/25')
ax2[0].plot(timevectm2, Sdilm2,label='Δx = Δy = L/50')
ax2[0].plot(timevectm3, Sdilm3,label='Δx = Δy = L/100')
ax2[1].plot(timevectr1, Sdilr1,label='Δx = Δy = L/25')
ax2[1].plot(timevectr2, Sdilr2,label='Δx = Δy = L/50')
ax2[1].plot(timevectr3, Sdilr3,label='Δx = Δy = L/100')
ax2[0].set_title('MacCormack method',fontsize=10)
ax2[0].set_ylabel('$\dot S_{dil}$ ($m^{2}/s^{3}$)')
ax2[1].set_title('Rusanov method',fontsize=10)
ax2[1].set_xlabel('t (s)')
ax2[1].set_ylabel('$\dot S_{dil}$ ($m^{2}/s^{3}$)')
ax2[0].grid(True)
ax2[1].grid(True)
ax2[1].legend(loc='best')
ax2[0].set_xlim([0, 0.29])
ax2[1].set_xlim([0, 0.29])
for i, ax in enumerate(ax2):
    ax.text(1.025, 0.5, string.ascii_uppercase[i], transform=ax.transAxes,size=10)
    ax.text(1.0175, 0.5, '(  )', transform=ax.transAxes,size=10)
ax2[0].legend(loc='best')

fig5, ax2 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig5.suptitle('Figure 10: Time history of $\dot S_{tor}$, $x_s$ = L/4, t = (L/2)/S',fontsize=10)
ax2[0].plot(timevectm1, Storm1,label='Δx = Δy = L/25')
ax2[0].plot(timevectm2, Storm2,label='Δx = Δy = L/50')
ax2[0].plot(timevectm3, Storm3,label='Δx = Δy = L/100')
ax2[1].plot(timevectr1, Storr1,label='Δx = Δy = L/25')
ax2[1].plot(timevectr2, Storr2,label='Δx = Δy = L/50')
ax2[1].plot(timevectr3, Storr3,label='Δx = Δy = L/100')
ax2[0].set_title('MacCormack method',fontsize=10)
ax2[0].set_ylabel('$\dot S_{tor}$ ($m^{2}/s^{3}$)')
ax2[1].set_title('Rusanov method',fontsize=10)
ax2[1].set_xlabel('t (s)')
ax2[1].set_ylabel('$\dot S_{tor}$ ($m^{2}/s^{3}$)')
ax2[0].grid(True)
ax2[1].grid(True)
ax2[1].legend(loc='best')
ax2[0].set_xlim([0, 0.29])
ax2[1].set_xlim([0, 0.29])
for i, ax in enumerate(ax2):
    ax.text(1.025, 0.5, string.ascii_uppercase[i], transform=ax.transAxes,size=10)
    ax.text(1.0175, 0.5, '(  )', transform=ax.transAxes,size=10)
ax2[0].legend(loc='best')

origin = 'lower'
fig6, ax1 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig6.suptitle('Figure 13: Comparison of $\dot S_{dil}$ contours, Δx = Δy = L/100, $x_s$ = L/4, t = (L/2)/S',fontsize=10)
ax1[0].set_title('MacCormack method',fontsize=10)
CS1 = ax1[0].contourf(x3, y3, dilatation_finm3, 10, origin=origin)
ax1[0].set_ylabel('y (m)')
CS2 = ax1[1].contourf(x3, y3, dilatation_finr3, 10 ,origin=origin)
ax1[1].set_title('Rusanov method',fontsize=10)
cbar1 = fig6.colorbar(CS1,ax=ax1[0],location='bottom',pad=0)
cbar2 = fig6.colorbar(CS2,ax=ax1[1],location='bottom',pad=0)
ax1[1].set_xlabel('x (m)')
ax1[1].set_ylabel('y (m)')
for i, ax in enumerate(ax1):
    ax.text(1.025, 0.5, string.ascii_uppercase[i], transform=ax.transAxes,size=10)
    ax.text(1.0175, 0.5, '(  )', transform=ax.transAxes,size=10)
    
origin = 'lower'
fig7, ax1 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig7.suptitle('Figure 14: Comparison of $\dot S_{tor}$ contours, Δx = Δy = L/100, $x_s$ = L/4, t = (L/2)/S',fontsize=10)
ax1[0].set_title('MacCormack method',fontsize=10)
CS1 = ax1[0].contourf(x3, y3, baroclinic_finm3, 10, origin=origin)
ax1[0].set_ylabel('y (m)')
CS2 = ax1[1].contourf(x3, y3, baroclinic_finr3, 10 ,origin=origin)
ax1[1].set_title('Rusanov method',fontsize=10)
cbar1 = fig7.colorbar(CS1,ax=ax1[0],location='bottom',pad=0)
cbar2 = fig7.colorbar(CS2,ax=ax1[1],location='bottom',pad=0)
ax1[1].set_xlabel('x (m)')
ax1[1].set_ylabel('y (m)')
for i, ax in enumerate(ax1):
    ax.text(1.025, 0.5, string.ascii_uppercase[i], transform=ax.transAxes,size=10)
    ax.text(1.0175, 0.5, '(  )', transform=ax.transAxes,size=10)

maxes = [np.amax(omegam1),np.amax(omegam2),np.amax(omegam3),np.amax(omegar1),np.amax(omegar2),np.amax(omegar3)]