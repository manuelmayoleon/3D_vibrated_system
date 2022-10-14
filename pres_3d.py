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


class Panel:
    def __init__(self,sigma=1.0,w=0.001,m=1.0):
        self.sigma = sigma
        self.w=  w
        self.m = m 
        
    def  betas(self,alpha):   
        return  10*(1-alpha)/(1+3*alpha)
    
    def htilde(self,densa):
        return densa*self.sigma**2

    def v1s(self,densa,alpha):
        
        return 4 * self.w * np.sqrt( 1 + self.betas(alpha) ) / (np.sqrt(np.pi) *
                    ( np.sqrt((1+self.betas(alpha))**2 + (4*np.sqrt(2)/15)*(1+alpha)*np.sqrt(1+self.betas(alpha))*(10+17*self.betas(alpha)-alpha*(10+9*self.betas(alpha)))*self.htilde(densa) )  -(1+self.betas(alpha) ) )  )
    def v1s_wdes(self,densa,alpha):
            
        return 6 * self.w * ( 1 + self.betas(alpha) ) / (np.sqrt(2*np.pi) * (1+alpha) * 
                     (2*(1-alpha)+(17-9*alpha)*self.betas(alpha)/5)*self.htilde(densa)       )

    def p(self,densa,alpha):
        
        return densa*self.m*0.5* (self.v1s(densa,alpha))**2
    def pp(self,densa,alpha):
            
        return densa*self.m*0.5* (self.v1s_wdes(densa,alpha))**2
fig23=plt.figure()


densa = np.linspace(0.0001, 0.02, 10000)
# densa  = 0.02 
alpha = 0.8
# alpha = np.linspace(0.01,0.999,10000)

w = 0.001

# plt.plot(alpha,Panel(1.0,0.001,1.0).p(densa,alpha),linewidth=1.5,linestyle=":",color="C0",label=" $p$  "  )
plt.plot(densa,Panel(1.0,w,1.0).p(densa,alpha),linewidth=1.5,color="C0",label=" $p$  "  )
plt.plot(densa,Panel(1.0,w,1.0).pp(densa,alpha),linewidth=1.5,linestyle=":",color="C1",label=" $p$ approx  "  )

# plt.yscale("log")
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

plt.ylabel ( r' $p$ ',rotation=0.0,fontsize=30)

plt.xlabel( r' $N/L^2$ ', fontsize=30)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.title ( r' \textbf {Presión como función de la densidad}' ,fontsize=40)
plt.legend(loc=0,fontsize=30)

plt.show()