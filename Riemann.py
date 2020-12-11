#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 17:24:27 2019

@author: ped3

Calculation of exact Riemann problem for an arbitrary left and right states.
Much of the algorithm is from Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics"
(see pgs. 119-128 in 2nd edition)

Using the variable notation from Toro.... 

gam = gamma 
gamp1 = gamma + 1
gamm1 = gamma - 1
uR = velocity in right (R) state
uL = velocity in left (L) state
pLR = pressure in L or R states
rhoLR = density in L or R state
pstar = pressure across contact surface
ustart = velocity across contact surface
astarL = speed of sound left of contact, a*L
astarR = speed of sound right of contact, a*R
SL = shock running left
SR = shock running right
SHL = head of left running rarefaction fan, uL-aL
STL = tail of left running rarefaction fan, u* - a*L
SHR = head of right running rarefaction fan, uR-aR
STR = tail of right running rarefaction fan, u*-a*R

"""

class Riemann():
    
    EPS=1.e-6
    MAXIT=100
    
    def f_and_fprm_shock(self, pstar,pLR,rhoLR,gam,gamm1,gamp1):
          # compute value of pressure function for shock
     import math
     A = 2./gamp1/rhoLR
     B = gamm1*pLR/gamp1
     sqrtterm = math.sqrt(A/(pstar+B))
     return[(pstar-pLR)*sqrtterm,sqrtterm*(1.-0.5*(pstar-pLR)/(B+pstar))]
        
    def f_and_fprm_rarefaction(self, pstar, pLR, aLR, gam, gamm1, gamp1):
     # compute value of pressure function for rarefaction
     return[((2.*aLR)/gamm1)*(pow(pstar/pLR,0.5*gamm1/gam)-1.),(aLR/pLR/gam)*pow(pstar/pLR,-0.5*gamp1/gam)]
        
    def find_star_state(self, gam, rhoL, pL, aL, uL, rhoR, pR, aR, uR):
      
      import math
      import numpy as np
      
      # first guess pstar based on two-rarefacation approximation
      pstar = aL+aR - 0.5*(gam-1.)*(uR-uL)
      pstar = pstar / (aL/pow(pL,0.5*(gam-1.)/gam) + aR/pow(pR,0.5*(gam-1.)/gam) )    
      pstar = pow(pstar,2.*gam/(gam-1.))
      gamm1 = gam-1.
      gamp1 = gam+1.
        
      if pstar<=pL: 
        f_L = self.f_and_fprm_rarefaction(pstar,pL,aL,gam,gamm1,gamp1)
      else:
        f_L = self.f_and_fprm_shock(pstar,pL,rhoL,gam,gamm1,gamp1)
      if (pstar<=pR):
        f_R = self.f_and_fprm_rarefaction(pstar,pR,aR,gam,gamm1,gamp1)
      else:
        f_R = self.f_and_fprm_shock(pstar,pR,rhoR,gam,gamm1,gamp1)
      delu = uR-uL
        
      if (f_L[0]+f_R[0]+delu)>self.EPS:
            # iterate using Newton-Rapson
        for it in range(0,self.MAXIT): 
          pold = pstar 
          pstar = pold - (f_L[0]+f_R[0]+delu)/(f_L[1]+f_R[1])
          if pstar<0.: 
             pstar=self.EPS
#          err = 2.*abs(pstar-pold)/(pstar+pold);
#          print ("it, err="+str(it)+" "+str(err))
          if (2.*abs(pstar-pold)/(pstar+pold)<self.EPS):
              break
          else:
             if pstar<=pL: 
                f_L = self.f_and_fprm_rarefaction(pstar,pL,aL,gam,gamm1,gamp1)
             else:
                f_L = self.f_and_fprm_shock(pstar,pL,rhoL,gam,gamm1,gamp1)
             if (pstar<=pR):
                f_R = self.f_and_fprm_rarefaction(pstar,pR,aR,gam,gamm1,gamp1)
             else:
                f_R = self.f_and_fprm_shock(pstar,pR,rhoR,gam,gamm1,gamp1)
        if it>self.MAXIT:
          print("error in Riemann.find_pstar")
          print("did not converage for pstar")
        
        # determine rest of star state
      ustar = 0.5*(uL+uR+f_R[0]-f_L[0])
      
      # intialize variables to something crazy so code flags if not set right.
      rhostarL = np.nan
      rhostarR = np.nan
      SL=np.nan
      SR=np.nan
      SHL=np.nan
      STL=np.nan
      SHR=np.nan
      STR=np.nan
      pratio=np.nan
      astarL=np.nan
      astarR=np.nan
        
        # left star state
      pratio = pstar/pL
      if pstar<=pL: # rarefaction
            rhostarL = rhoL*pow(pratio,1./gam)
            astarL = aL*pow(pratio,0.5*gamm1/gam)
            SHL = uL-aL
            STL = ustar - astarL
      else: #shock
            rhostarL = rhoL*(pratio+gamm1/gamp1)/(gamm1*pratio/gamp1+1.)
            SL = uL - aL*math.sqrt(0.5*gamp1/gam*pratio+0.5*gamm1/gam)
        
        # right star state
      pratio = pstar/pR
      if pstar<=pR: # rarefaction
            rhostarR = rhoR*math.pow(pratio,1./gam)
            astarR = aR*math.pow(pratio,0.5*gamm1/gam)
            SHR = uR+aR
            STR = ustar + astarR
      else: #shock
            rhostarR = rhoR*(pratio+gamm1/gamp1)/(gamm1*pratio/gamp1+1.)
            SR = uR + aR*math.sqrt(0.5*gamp1/gam*pratio+0.5*gamm1/gam)    
        
      return [pstar,ustar,rhostarL,astarL,SL,SHL,STL,rhostarR,astarR,SR,SHR,STR,
                            rhoL,pL,aL,uL,rhoR,pR,aR,uR,gam,gamm1,gamp1]
    
    def sample(self, state, xDt):
        
      # sample the Riemann solution state
        pstar = state[0]
        ustar = state[1]
        rhostarL = state[2]
        astarL = state[3]
        SL = state[4]
        SHL = state[5]
        STL = state[6]
        rhostarR = state[7]
        astarR = state[8]
        SR = state[9]
        SHR = state[10]
        STR = state[11]
        rhoL = state[12]
        pL = state[13]
        aL = state[14]
        uL = state[15]
        rhoR = state[16]
        pR = state[17]
        aR = state[18]
        uR = state[19]
        gam = state[20]
        gamm1 = state[21]
        gamp1 = state[22]
        
        if (xDt <= ustar): # left of contact surface
          if (pstar<=pL): # rarefaction
            if (xDt<= SHL):
              rho = rhoL
              p = pL
              u = uL
            elif (xDt <=STL): # SHL < x/t < STL
              tmp = 2./gamp1 + (gamm1/gamp1/aL)*(uL-xDt)
              rho = rhoL*pow(tmp,2./gamm1)
              u = (2./gamp1)*(aL + 0.5*gamm1*uL+xDt)
              p = pL*pow(tmp,2.*gam/gamm1)
            else: # STL < x/t < u*
              rho = rhostarL
              p = pstar
              u = ustar
          else: # shock
             if xDt<= SL: # xDt < SL
                 rho = rhoL
                 p = pL
                 u = uL
             else: # SL < xDt < ustar
                 rho = rhostarL
                 p = pstar
                 u = ustar
        else: # right of contact surface 
          if pstar<=pR: # rarefaction
                if xDt>= SHR:
                  rho = rhoR
                  p = pR
                  u = uR
                elif (xDt >= STR): # SHR < x/t < SHR
                  tmp = 2./gamp1 - (gamm1/gamp1/aR)*(uR-xDt)
                  rho = rhoR*pow(tmp,2./gamm1)
                  u = (2./gamp1)*(-aR + 0.5*gamm1*uR+xDt)
                  p = pR*pow(tmp,2.*gam/gamm1)
                else: # u* < x/t < STR
                  rho = rhostarR
                  p = pstar
                  u = ustar               
          else: # shock
                 if (xDt>= SR): # xDt > SR
                     rho = rhoR
                     p = pR
                     u = uR
                 else: # ustar < xDt < SR
                     rho = rhostarR
                     p = pstar
                     u = ustar
        e=p/gamm1/rho;
        return [rho,p,u,e]
 
    def exactsol(self,gam,xs,time,case,filename):
 
     import numpy as np
     import math
     import matplotlib.pyplot as plt
     
     #analytical Riemann result
     L = 1.
     NX=1000
     dx= L/NX
     xexact = np.arange(0,L+dx,dx)
     nx = np.size(xexact)
     rhoexact=np.zeros(nx)
     pexact=np.zeros(nx)
     uexact=np.zeros(nx)
     eexact=np.zeros(nx)
     rhoL,uL,pL,rhoR,uR,pR,max_time = self.get_cases()
     aL = math.sqrt(gam*pL[case]/rhoL[case]) 
     aR = math.sqrt(gam*pR[case]/rhoR[case])
     star_state = self.find_star_state(gam,rhoL[case],pL[case],aL,uL[case],rhoR[case],pR[case],aR,uR[case])
#     print("pstar ="+str(star_state[0])+"    ustar ="+str(star_state[1]))
     for i in range(0,nx):
       xDt = (xexact[i]-xs)/time
       out = self.sample(star_state,xDt)
       rhoexact[i]=out[0]
       pexact[i]=out[1] 
       uexact[i]=out[2]
       eexact[i]=out[3]
       
     return xexact,rhoexact,uexact,pexact,eexact
     
#     plt.savefig(filename+'.png',bbox_inches='tight')
     
    def plot_compare(self,x,rho,p,u,e,gam,time,case,filename):
 
     import numpy as np
     import math
     import matplotlib.pyplot as plt
     
     #analytical Riemann result
     L = 1.
     NX=1000
     dx= L/NX
     xexact = np.arange(0,L+dx,dx)
     nx = np.size(xexact)
     rhoexact=np.zeros(nx)
     pexact=np.zeros(nx)
     uexact=np.zeros(nx)
     eexact=np.zeros(nx)
     rhoL,uL,pL,rhoR,uR,pR,max_time = self.get_cases()
     aL = math.sqrt(gam*pL[case]/rhoL[case]) 
     aR = math.sqrt(gam*pR[case]/rhoR[case])
     star_state = self.find_star_state(gam,rhoL[case],pL[case],aL,uL[case],rhoR[case],pR[case],aR,uR[case])
     print("pstar ="+str(star_state[0])+"    ustar ="+str(star_state[1]))
     for i in range(0,nx):
       xDt = (xexact[i]-L/2)/max_time[case]
       out = self.sample(star_state,xDt)
       rhoexact[i]=out[0]
       pexact[i]=out[1] 
       uexact[i]=out[2]
       eexact[i]=out[3]
       
     #plot rho,u,p,e comparisons
     f, ax = plt.subplots(2,2,figsize=(12 ,5))
     
      # rho
     ax[0][0].plot(xexact,rhoexact,label='t='+str(max_time[case]),linestyle='--',color='black')
     ax[0][0].plot(x,rho,label='t='+str(time),linestyle='-',color='black')
     ax[0][0].set_xlabel('$x$',size=20)
     ax[0][0].set_ylabel(r'$\rho (kg/m^3)$',size=20)
     ax[0][0].grid()
     ax[0][0].legend(fontsize=16)
    
    # u
     ax[0][1].plot(xexact,uexact,label='t='+str(max_time[case]),linestyle='--',color='black')
     ax[0][1].plot(x,u,label='t='+str(time),linestyle='-',color='black')
     ax[0][1].set_xlabel('$x$',size=20)
     ax[0][1].set_ylabel('$u (m/s)$',size=20)
     ax[0][1].grid()
     ax[0][1].legend(fontsize=16)
     
     # p
     ax[1][0].plot(xexact,pexact,label='t='+str(max_time[case]),linestyle='--',color='black')
     ax[1][0].plot(x,p,label='t='+str(time),linestyle='-',color='black')
     ax[1][0].set_xlabel('$x(m)$',size=20)
     ax[1][0].set_ylabel('$p(Pa)$',size=20)
     ax[1][0].grid()
     ax[1][0].legend(fontsize=16)
    
    # ener
     ax[1][1].plot(xexact,eexact,label='t='+str(max_time[case]),linestyle='--',color='black')
     ax[1][1].plot(x,e,label='t='+str(time),linestyle='-',color='black')
     ax[1][1].set_xlabel('$x(m)$',size=20)
     ax[1][1].set_ylabel('$e (J/kg-K)$',size=20)
     ax[1][1].grid()
     ax[1][1].legend(fontsize=16)
     
     plt.show()
    
#     plt.savefig(filename+'.png',bbox_inches='tight')
      
     
    def get_cases(self):
         
        rhoL = dict()
        uL = dict()
        pL = dict()
        rhoR = dict()
        uR = dict()
        pR = dict()
        max_time = dict()
       
        # case 1 - Sod problem
        rhoL[1]=1.0 
        uL[1]=0.0 
        pL[1]=1.0 
        rhoR[1]=0.125 
        uR[1]=0.0 
        pR[1]=0.1 
        max_time[1] = 0.14268138124362462
        
        # case 2 - 123 problem - expansion left and expansion right
        rhoL[2]=1.0 
        uL[2]=-2.
        pL[2]=0.4
        rhoR[2]=1.0
        uR[2]=2.
        pR[2]=0.4
        max_time[2] = 0.15
        
        # case 3 - blast problem - shock right, expansion left
        rhoL[3]=1.0 
        uL[3]=0.0
        pL[3]=1000.
        rhoR[3]=1.0
        uR[3]=0.
        pR[3]=0.01
        max_time[3] = 0.012
        
        # case 4 - blast problem - shock left, expansion right
        rhoL[4]=1.0 
        uL[4]=0.0
        pL[4]=0.01
        rhoR[4]=1.0
        uR[4]=0.
        pR[4]=100.
        max_time[4] = 0.035
        
        # case 5 - shock collision - shock left and shock right
        rhoL[5]=5.99924
        uL[5]=19.5975
        pL[5]=460.894
        rhoR[5]=5.99242
        uR[5]=-6.19633
        pR[5]=46.0950
        max_time[5] = 0.035
        
        return(rhoL,uL,pL,rhoR,uR,pR,max_time)
        
#########################
# example of usage.... 
#########################
if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case= 5 #case numbers are summarized above 
 gam=1.4
 R.plot(gam,max_time[case],case,"rusanov"+str(case))        