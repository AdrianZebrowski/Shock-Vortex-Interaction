import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import string
from project6part2 import project6part2

L = 1
dx1 = L/150
time = 2
N1 = int(L/dx1+1)

x1,y1,omegaex,omegam1,omegar1,timevectm1,timevectr1,ensm1,ensr1,Sdilm1,Sdilr1,Storm1,Storr1,dilatation_finm1,dilatation_finr1,baroclinic_finm1,baroclinic_finr1 = project6part2(dx1,time)

maxes = [np.amax(omegam1),np.amax(omegar1)]