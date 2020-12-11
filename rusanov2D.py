# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 22:26:59 2020

@author: az
"""
def rusanov2D(N,L,t_final,CFL,u,v,rho,p,gamma,R):
    import numpy as np
    # First, establish some storage vectors 
    omega = np.zeros((N,N)) # vorticity using centered differencing
    dilatation = np.zeros((N,N))
    baroclinic = np.zeros((N,N))
    dilatation_fin = np.zeros((N,N))
    baroclinic_fin = np.zeros((N,N))
    dx = L/(N-1)
    dy = L/(N-1)   
    # Establish some more storage vectors for Q, F, G, Fhalf, Ghalf
    Q = np.zeros((4,N,N))
    #Q_update = np.zeros((4,N,N))
    F = np.zeros((4,N,N))
    F_half = np.zeros((4,N,N))
    G = np.zeros((4,N,N))
    G_half = np.zeros((4,N,N))
    x = np.linspace(0,L,N) # vectors for storing x and y coordinates
    y = np.linspace(0,L,N) 
    ens1D = np.zeros(N)
    baroclinic1D = np.zeros(N)
    dilatation1D = np.zeros(N)
    
    # Calculate a, e, et, H at the initial condition for use in flux vectors
    e = p/((gamma-1)*rho)
    et = e+(u**2+v**2)/2
    H = e+p/rho+(u**2+v**2)/2
    a = np.sqrt(gamma*p/rho)
    
    # Set up initial time
    t = 0
    timestep = 0
    timevect = []
    ens = []
    Sdil = []
    Stor = []
    # Calculate dtx and dty similarly to the 1D method. I suspect this is where my problem is, since the solution appears to be valid if CFL number is made low enough but blows up and becomes nonsense even with CFL = 0.5 or so
    dtx = CFL*dx/np.max(np.abs(u)+a)
    dty = CFL*dy/np.max(np.abs(v)+a)
    
    # Take the smallest value of the three and use that as the timestep
    dt = min(dtx,dty)
            
    #calculate initial Q and F vectors
    Q[0,:,:] = rho
    Q[1,:,:] = rho*u
    Q[2,:,:] = rho*v
    Q[3,:,:] = rho*et
    
    F[0,:,:] = rho*u
    F[1,:,:] = rho*u**2+p
    F[2,:,:] = rho*u*v
    F[3,:,:] = rho*u*H
    
    G[0,:,:] = rho*v
    G[1,:,:] = rho*u*v
    G[2,:,:] = rho*v**2+p
    G[3,:,:] = rho*v*H
    
    # This loop calculates the solution
    while t <= t_final:
         # This loop calculates vorticity dilatation and baroclinic 
        for j in range(0,N):
            for i in range(0,N):
                
                # Same implementation of periodic BCs
                im=i-1
                ip=i+1
                jm=j-1
                jp=j+1
                
                if i == 0:
                    im = 1
                    #im = N-2
                if i == N-1:
                    #ip = 1
                    ip = N-2
                if j == 0:
                    jm = N-2
                if j == N-1:
                    jp = 1
                
                # Calculate vorticity numerically using central differencing
                omega[j,i] = (v[j,ip]-v[j,im])/(2*dx)-(u[jp,i]-u[jm,i])/(2*dy)
                dilatation[j,i] = -(u[j,ip]-u[j,im])/(2*dx)+(v[jp,i]-v[jm,i])/(2*dy)
                baroclinic[j,i] = ((rho[j,ip]-rho[j,im])/(2*dx))*((p[jp,i]-p[jm,i])/(2*dy))-((rho[jp,i]-rho[jm,i])/(2*dy))*((p[j,ip]-p[j,im])/(2*dx))
        
        for j in range(0,N):
            ens1D[j] = np.trapz(omega[j,:]**2,x)
            dilatation1D[j] = np.trapz(dilatation[j,:],x)
            baroclinic1D[j] = np.trapz(baroclinic[j,:],x)
            
        dilatation2D = np.trapz(dilatation1D,y)
        baroclinic2D = np.trapz(baroclinic1D,y)
        ens2D = np.trapz(ens1D,y)
        
        ens.append(ens2D) 
        Sdil.append(dilatation2D)
        Stor.append(baroclinic2D)
        timevect.append(t)
        
        for j in range(0,N):
            for i in range(0,N): # PREDICTOR
                
                # This is just using a dummy index so I can implement the boundary conditions a bit more elegantly - I can just make ip=1 (second node) if i = N-1 (last node)
                im=i-1 
                ip=i+1
                jm=j-1
                jp=j+1
                
                # Periodic BCs implemented here using a series of if statements
                if i == 0:
                    im = 1
                    #im = N-2
                if i == N-1:
                    #ip = 1
                    ip = N-2
                if j == 0:
                    jm = N-2
                if j == N-1:
                    jp = 1
                            
                sx = np.maximum(np.abs(u[j,i])+a[j,i],np.abs(u[j,ip])+a[j,ip])
                sy = np.maximum(np.abs(v[j,i])+a[j,i],np.abs(v[jp,i])+a[jp,i])
                
                F_half[:,j,i] = 0.5*(F[:,j,i]+F[:,j,ip]-sx*(Q[:,j,ip]-Q[:,j,i])) # calculate F_half
                G_half[:,j,i] = 0.5*(G[:,j,i]+G[:,jp,i]-sy*(Q[:,jp,i]-Q[:,j,i])) # calculate G_half           
        
        for j in range(0,N):
            for i in range(0,N):
                
                im=i-1
                ip=i+1
                jm=j-1
                jp=j+1
            
                if i == 0:
                    im = 1
                    #im = N-2
                if i == N-1:
                    #ip = 1
                    ip = N-2
                if j == 0:
                    jm = N-2
                if j == N-1:
                    jp = 1
                
                Q[:,j,i] = Q[:,j,i]-dt*((F_half[:,j,i]-F_half[:,j,im])/dx+(G_half[:,j,i]-G_half[:,jm,i])/dy) # update value of Q
        
        # Calculate all the update values of properties we care about
        rho = Q[0,:,:]
        u = Q[1,:,:]/rho
        v = Q[2,:,:]/rho
        et = Q[3,:,:]/rho
        e = et-(u**2+v**2)/2
        p = e*(gamma-1)*rho # calculate values of rho, u, e, and p based on predictor Q
        T = p/(rho*R)
        H = e+p/rho +(u**2+v**2)/2
        a = np.sqrt(gamma*p/rho)
    
        F[0,:,:] = rho*u
        F[1,:,:] = rho*u**2+p
        F[2,:,:] = rho*u*v
        F[3,:,:] = rho*u*H
    
        G[0,:,:] = rho*v
        G[1,:,:] = rho*u*v
        G[2,:,:] = rho*v**2+p
        G[3,:,:] = rho*v*H
        
        dtx = CFL*dx/np.max(np.abs(u)+a)
        dty = CFL*dy/np.max(np.abs(v)+a)
        dt = min(dtx,dty)
        t = t+dt
        timestep = timestep+1
        print(t)
        
    dilatation_fin = ((2*omega**2)*dilatation)
    baroclinic_fin = ((omega/(2*rho**2))*(baroclinic))
    
    return u,v,rho,p,T,e,omega,timevect,ens,Sdil,Stor,dilatation_fin,baroclinic_fin
        
        