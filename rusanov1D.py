def rusanov1D(case,t_final,c_max,dx,gamma):
    import numpy as np

    if case == "Case 1":
        # case 1 - Sod problem
        rhoL=1.0 
        uL=0.0 
        pL=1.0 
        rhoR=0.125 
        uR=0.0 
        pR=0.1 
    if case == "Case 2":
        # case 2 - 123 problem - expansion left and expansion right
        rhoL=1.0 
        uL=-2.0
        pL=0.4
        rhoR=1.0
        uR=2.0
        pR=0.4
    if case == "Case 3":
        # case 3 - blast problem - shock right, expansion left
        rhoL=1.0 
        uL=0.0
        pL=1000.
        rhoR=1.0
        uR=0.
        pR=0.01
    if case == "Case 4":
        # case 4 - blast problem - shock left, expansion right
        rhoL=1.0 
        uL=0.0
        pL=0.01
        rhoR=1.0
        uR=0.
        pR=100.
    if case == "Case 5":
        # case 5 - shock collision - shock left and shock right
        rhoL=5.99924
        uL=19.5975
        pL=460.894
        rhoR=5.99242
        uR=-6.19633
        pR=46.0950
        
    L = 1.0
    ix = int(L/dx+1)
    half = int(np.floor(ix/4))
    x = np.linspace(0,L,ix)

    # storage vectors
    rho = np.zeros(ix)
    u = np.zeros(ix)
    p = np.zeros(ix)
    e = np.zeros(ix)
    Q = np.zeros((3,ix))
    Q_update = np.zeros((3,ix))
    F = np.zeros((3,ix))
    F_half = np.zeros((3,ix))
    
    #establish initial condition vectors
    rho[0:half] = rhoL
    rho[half:ix] = rhoR
    u[0:half] = uL
    u[half:ix] = uR
    p[0:half] = pL
    p[half:ix] = pR
    e = p/((gamma-1)*rho)
    a = np.sqrt(gamma*p/rho)
    
    t = 0
    dt = c_max*dx/np.max(np.abs(u)+a)
    
    #calculate initial Q and F vectors
    Q[0,:] = rho
    Q[1,:] = rho*u
    Q[2,:] = rho*(e+u**2/2)
    
    F[0,:] = rho*u
    F[1,:] = rho*u**2+p
    F[2,:] = rho*u*(e+p/rho+u**2/2)
        
    while t <= t_final:
        for i in range(0,ix): # F_half  
                
            im=i-1 
            ip=i+1
                
            # Neumann BCs implemented here using a series of if statements
            if i == 0:
                im = 1
            if i == ix-1:
                ip = ix-2
                    
            s = np.maximum(np.abs(u[i])+a[i],np.abs(u[ip])+a[ip])
            F_half[:,i] = 0.5*(F[:,i]+F[:,ip]-s*(Q[:,ip]-Q[:,i])) # calculate F_half
        
        for i in range(0,ix): # Q_update
            im=i-1 
            ip=i+1
                
            # Neumann BCs implemented here using a series of if statements
            if i == 0:
                im = 1
            if i == ix-1:
                ip = ix-2
                    
            Q_update[:,i] = Q[:,i]-(dt/dx)*(F_half[:,i]-F_half[:,im])
        
        Q = Q_update[:,:]
        rho = Q[0,:]
        u = Q[1,:]/rho
        e = Q[2,:]/rho-u**2/2
        p = e*(gamma-1)*rho
        a = np.sqrt(gamma*p/rho)
    
        F[0,:] = rho*u
        F[1,:] = rho*u**2+p
        F[2,:] = rho*u*(e+p/rho+u**2/2)
        
        dt = c_max*dx/np.max(np.abs(u)+a)
        t = t+dt
        s = dx/dt
        
    return rho, u, p, e, x