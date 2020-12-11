# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 22:26:59 2020

@author: az
"""
def init2Dvortex(N,L,Rc,xc,yc,Tinf,pinf,gamma,R,Mac,uinf,vinf):
    import numpy as np 

    # First, establish some storage vectors 
    u = np.zeros((N,N))
    v = np.zeros((N,N))
    T = np.zeros((N,N))
    
    omega = np.zeros((N,N)) # vorticity using centered differencing
    omegaex = np.zeros((N,N)) # vorticity using analytical expression
    
    dx = L/(N-1)
    dy = L/(N-1)
    
    x = np.linspace(0,L,N) # vectors for storing x and y coordinates
    y = np.linspace(0,L,N) 
    
    Tc = Tinf/(1+((gamma-1)/2)*Mac**2) # Tc and ac are calculated since we need ac to use equation 4a and 4b
    ac = np.sqrt(gamma*R*Tc)

    # This nested loop initializes the flow field
    for j in range(0,N):
        for i in range(0,N):
            rstar = np.sqrt((x[i]-xc)**2+(y[j]-yc)**2)/Rc #Calculate rstar, ystar, and xstar here since they are used repeatedly
            ystar = (y[j]-yc)/Rc
            xstar = (x[i]-xc)/Rc
            rstarexp = np.exp((1-(rstar)**2)/2) #This is the exponential function and everything inside it
            u[j,i] = uinf-Mac*ac*ystar*rstarexp #Calculate u and v using equation 4 at every node - note that the indexing is j,i since rows in the array correspond to the y-axis
            v[j,i] = vinf+Mac*ac*xstar*rstarexp 
            T[j,i] = Tinf*(1-(((gamma-1)/2)*Mac**2/(1+((gamma-1)/2)*Mac**2))*rstarexp*rstar**2) # Calculate temperature using eq. 6
            omegaex[j,i] = (Mac*ac/Rc)*rstarexp*(2-xstar**2-ystar**2) # Calculate vorticity at each node using analytical expression
            
    p = pinf*(T/Tinf)**((gamma)/(gamma-1)) # Can now calculate p, rho, e based on the other values we have
    rho = p/(R*T)
    e = p/((gamma-1)*rho)
    # Initialization complete
    
    # This loop calculates vorticity
    for j in range(0,N):
        for i in range(0,N):
            
            # Same implementation of periodic BCs
            im=i-1
            ip=i+1
            jm=j-1
            jp=j+1
            
            if i == 0:
                im = N-2
            if i == N-1:
                ip = 1
            if j == 0:
                jm = N-2
            if j == N-1:
                jp = 1
            
            # Calculate vorticity numerically using central differencing
            omega[j,i] = (v[j,ip]-v[j,im])/(2*dx)-(u[jp,i]-u[jm,i])/(2*dy)
    
    # Calculate L2 vorticity error and circulation error here        
    L2omega = np.sqrt(np.sum((omega-omegaex)**2))/(N**2)
    circ1D = np.zeros(N)
    eps1D = np.zeros(N)
    
    for j in range(0,N):
        circ1D[j] = np.trapz(omega[j,:],x)
        eps1D[j] = np.trapz(omega[j,:]**2,x)
        
    circ2D = np.trapz(circ1D,y)
    eps2D = np.trapz(eps1D,y)     
    circerr = circ2D/np.sqrt(eps2D)
    
    return u, v, rho, p, T, e, omega, omegaex, L2omega, circerr
        
        