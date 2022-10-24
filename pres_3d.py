import numpy as np
import pandas as pd
import scipy.special as special
from numba import jit
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)

# import scipy.special as special

# Import math Library
import math 

from sympy import *

# PARA USAR A LA HORA DE GUARDAR DOS COLUMNAS EN UN ARCHIVO
import csv

tvsd= pd.read_csv("temp_vs_DenPart-1.dat",header=0,sep = '\s+', names=["den","ts","err"])

# tvsd.columns = tvsd.iloc[0]

# T = tvsd['Ts1']
print(tvsd['den'])


class Panel:
    def __init__(self,sigma=1.0,w=0.001,m=1.0,h = 29):
        self.sigma = sigma
        self.w=  w
        self.m = m 
        self.h = h
    def  betas(self,alpha):   
        return  10*(1-alpha)/(1+3*alpha)
    
    def htilde(self,dens):
        return self.h*dens*self.sigma**2

    def v1s(self,dens,alpha):
        
        return 4 * self.w * np.sqrt( 1 + self.betas(alpha) ) / (np.sqrt(np.pi) *
                    ( np.sqrt((1+self.betas(alpha))**2 + (4*np.sqrt(2)/15)*(1+alpha)*np.sqrt(1+self.betas(alpha))*(10+17*self.betas(alpha)-alpha*(10+9*self.betas(alpha)))*self.htilde(dens) )  -(1+self.betas(alpha) ) )  )
    def v1s_wdes(self,dens,alpha):
            
        return 6 * self.w * ( 1 + self.betas(alpha) ) / (np.sqrt(2*np.pi) * (1+alpha) * 
                     (2*(1-alpha)+(17-9*alpha)*self.betas(alpha)/5)*self.htilde(dens)       )

    def p(self,dens,alpha):
        
        return dens*self.m*0.5* (self.v1s(dens,alpha))**2
    def pp(self,dens,alpha):
            
        return dens*self.m*0.5* (self.v1s_wdes(dens,alpha))**2
fig23=plt.figure()


# dens = np.linspace(1.5, 4.5, 1000)

dens = np.linspace(0.005 , max(tvsd['den']) , 1000)
# dens  = 0.02 
w = 0.001
sigma = 1.0
h = 29*sigma
m = 1.0
alpha = 0.9
# alpha = np.linspace(0.80,0.975,10000)



# plt.plot(alpha,Panel(1.0,0.001,1.0).p(dens,alpha),linewidth=1.5,linestyle=":",color="C0",label=" $p$  "  )
# plt.plot(tvsd['den'],tvsd['Ts1'],marker= 'o',linestyle = 'None',color="C2",label=" $v_{1}^s/\omega$ (MD) "  )

plt.errorbar(tvsd['den'], tvsd['ts'], yerr=tvsd['err'], color='C2',marker="o",linestyle="",label="$n_z$ (MD)") 
plt.plot(dens,Panel(1.0,2*w,1.0,29).v1s(dens,alpha)/(2*w),linewidth=1.5,color="C0",label=" $v_{1}^s/\omega$  "  )
plt.plot(dens,Panel(1.0,2*w,1.0,29).v1s_wdes(dens,alpha)/(2*w),linewidth=1.5,linestyle=":",color="C1",label=" $v_{1}^s/\omega$ $O(w^2)$  "  )

# plt.yscale("log")
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

plt.ylabel ( r'  ',rotation=0.0,fontsize=30)

plt.xlabel( r' $n$ ', fontsize=30)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.title ( r' \textbf {Presión como función de la densidad}' ,fontsize=40)
plt.legend(loc=0,fontsize=30)

plt.show()