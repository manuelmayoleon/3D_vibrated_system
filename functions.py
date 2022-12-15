#functions

import numpy as np
import pandas as pd
import scipy.special as special
from numba import jit
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import random
# import scipy.special as special

# Import math Library
import math 

from sympy import *

import csv

from scipy.optimize import root, fsolve




class Panel1:
    def __init__(self,a=0.950,sigma=1.0):
        self.a = a
        self.sigma = sigma
    def htilde(self,dens,h):
            return h*dens*self.sigma**2
    def polinomio(self,b):
        return np.sqrt(b+1)*(4 * ( 1 + self.a ) * b ** 2 + ( 16 * self.a - 8 ) * b + ( 9 - 3 * self.a ) )/( ( 30 - 18 * self.a) * b + ( 9 -3 * self.a ) )
    def  bets_1(alpha):   
            return  10*(1-alpha)/(1+3*alpha)
            # return 32*(1-a)/(5+9*a)
    def  bets_2(alpha):   
            # return  10*(1-a)/(1+3*a)
            return 32*(1-alpha)/(5+9*alpha)
    def transcendental_eq(self,b):
        return  np.sqrt(b+1)*(4 * ( 1 + self.a ) * b ** 2 + ( 16 * self.a - 8 ) * b + ( 9 - 3 * self.a ) )/( ( 30 - 18 * self.a) * b + ( 9 -3 * self.a ) )- np.arcsinh(np.sqrt(b))/np.sqrt(b)
    def b1(alfa,b):
        return  - np.arcsinh( np.sqrt(b) ) * ( 12 * b + ( 9 - 3 * alfa ) ) / np.sqrt(b)  
    def b2(alfa,b):
            return np.sqrt(1+b )*(( 16 - 8 * alfa ) * b**2 + ( 22 - 14 * alfa ) * b + ( 9 - 3 * alfa ) )
    def btilde(alfa,b): 
          return (1+alfa)*(Panel1.b1(alfa,b) + Panel1.b2(alfa,b))/( 96 * b ) 
    def btilde(alfa,b): 
          return (1+alfa)*(Panel1.b1(alfa,b) + Panel1.b2(alfa,b))/( 96 * b )
    def b_coef(self,alfa,b,dens,h):
        return (1+b)/(np.sqrt(2*np.pi)*Panel1.btilde(alfa,b)*Panel1().htilde(dens,h))
    def v1s_w(self,alfa,b,dens,h):
        return 0.5 * Panel1().b_coef(alfa,b,dens,h)*( 1 + np.sqrt(1 + 8 / ( Panel1().b_coef(alfa,b,dens,h)*np.sqrt((1+b)*np.pi) ) ) )



    
    
class Panel2:
    def __init__(self,sigma=1.0,w=0.001,m=1.0,h = 29):
        self.sigma = sigma
        self.w=  w
        self.m = m 
        self.h = h
    def  betas(self,alpha):   
        return  10*(1-alpha)/(1+3*alpha)
    
    def htilde(self,dens):
        return self.h*dens*self.sigma**2

    def a1(self,dens,alpha):
        
        return (24* (np.sqrt(1+self.betas(alpha)))/(np.sqrt(2*np.pi)*(1+alpha)*self.htilde(dens)*(2*(1-alpha)+self.betas(alpha)*(17-9*alpha)/5)))
    def b1(self,dens,alpha):
        
        return (12* (np.sqrt(1+self.betas(alpha)))/(np.sqrt(2*np.pi)*(1+alpha)*self.htilde(dens)*(2*(1-alpha)+self.betas(alpha)*(17-9*alpha)/5)))

    def a2(self,dens,alpha):
            
        return  self.a1(dens,alpha) * np.sqrt(self.betas(alpha)+1)
    def vsw(self,dens,alpha):
        
        return  0.25*self.a1(dens,alpha)*self.a2(dens,alpha)+ 6.32455532033676e-8*np.sqrt(self.a1(dens,alpha)*(15625000000000.0*self.a1(dens,alpha)*self.a2(dens,alpha)**2 + 141047395886939.0))
    
    def vsw_mano(self,dens,alpha):
            
        return self.b1(dens,alpha)*np.sqrt(self.betas(alpha)+1)/2 *( 1 + np.sqrt( 1 + ( 8 ) / ( ( 1 + self.betas( alpha ) ) * self.b1( dens , alpha ) * np.sqrt( np.pi ) )  ) )
    