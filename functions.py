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
    
class Panel3:
    def __init__(self,sigma=1.0,w=0.001,m=1.0,h = 29):
        self.sigma = sigma
        self.w=  w
        self.m = m 
        self.h = h
        
    def gama(self,h,sigma):
        
        return (h -2*sigma)/(h-sigma)
    
    def eps(self,h,sigma):
        
        return (h -sigma)/(sigma)
    
    def gg(self,d):
    
        return -( 3 + d  )/ ( 3 + 4 * d + d ** 2 ) + ( ( 2 + d ) * np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 4 * special.gamma(2 + d / 2 ) )
    
    def hh(self,d):
        
            return - 2  / ( 3 + 4 * d + d ** 2 ) + ( np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 4 * special.gamma(2 + d / 2 ) )
            
    def gp(self,d):
    
            return - ( 6 + 2 * d  )  / ( ( 3 + d) * (- 1 + d ** 2 ) ) - (   np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 2 * (1 - d )  * special.gamma(1 + d / 2 ) )
    
    def hp(self,d):
        
            return - 8 / ( ( 3 + d) * (- 1 + d ** 2 ) ) - ( 3 * np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 4 * (1 - d )  * special.gamma(2 + d / 2 ) )
    
    def factor1(self,d):
        
        return np.sqrt( np.pi ) * special.gamma( ( 1 + d ) / 2 ) / ( 4 * special.gamma( 2 + d / 2  )   )
    
    def factor2(self,d):
            
        return np.sqrt( np.pi ) * special.gamma( ( d - 1 ) / 2 ) / ( 8 * special.gamma( 2 + d / 2  )   )
        
    def  betas(self,alpha):   
    
        return (10*(1-alpha)/(1+3*alpha))*(5/(self.eps(self.h,self.sigma))+8*self.gama(self.h,self.sigma))/(3/self.eps(self.h,self.sigma)+self.gama(self.h,self.sigma)*8)
   
    def  betas_h2(self,alpha):   
    
        return (5 / 3 )*(10*(1-alpha)/(1+3*alpha))
    
    def betas_d(self, alpha,d ):
        
        return ( 2 * ( 1 - alpha ) / (1 + 3 * alpha )  ) * ( ( 2 + d ) * self.gama(self.h,self.sigma) + self.gg(d)/( self.factor1( d ) * self.eps(self.h,self.sigma) ))/( self.gama(self.h,self.sigma) + self.hh(d)/( self.factor1( d ) * self.eps(self.h,self.sigma) ) ) 
    
    def solid_angle_surface(self,d):
        
        return 2 * ( np.sqrt( np.pi  ) ) ** (d-1) / special.gamma( (d-1) / 2 ) 
    
    def  ts_w_uds_d(self,alpha,dens,d):  
         
        return self.m * ( ( 16 * ( 1 + self.betas_d(alpha,d) ) *  special.gamma( 2 + d / 2   ) )  / (   special.gamma(  ( d - 1 ) / 2   ) * self.solid_angle_surface( d ) * ( 1 + alpha ) * dens * ( self.h - self.sigma ) * self.sigma **( d-1 )  *( (1 - alpha )* ( (2 + d ) * self.gama(self.h,self.sigma) + self.gp(d)/( self.eps(self.h,self.sigma) * self.factor2( d )  ) )  + 0.5 *  self.betas_d(alpha,d) * ( ( 5 + 4 * d - 9 * alpha ) * self.gama(self.h,self.sigma) - ( ( 1 + 3 * alpha )* self.hp(d) - 4 * self.gp(d) )/ ( self.eps(self.h,self.sigma) * self.factor2( d ) ) )  )  )   ) ** 2

        
    def  ts_w_uds(self,alpha,dens):   
    
        return self.m * ( 6 * ( 1 + self.betas(alpha) )  / ( np.sqrt(np.pi) * ( 1 + alpha ) * dens * ( self.h - self.sigma ) * self.sigma ** 2  *( (1 - alpha )* ( 2 * self.gama(self.h,self.sigma) + 0.5/self.eps(self.h,self.sigma) )  +  self.betas(alpha) * ( ( 17 - 9 * alpha ) * self.gama(self.h,self.sigma) + 3 * 0.5 *( 3 - alpha )/ self.eps(self.h,self.sigma) ) / 5 )  )   ) ** 2
    
    def  ts_w_uds_2sigma(self,alpha,dens):   
    
        return self.m * ( 12 * ( 1 + self.betas_h2(alpha) )  / ( np.sqrt(np.pi) * ( 1 + alpha ) * dens *  self.sigma ** 3  *( (1 - alpha ) +  (3/5) *  self.betas_h2(alpha) * ( ( 3 -  alpha ) )  )  )   ) ** 2
    

class Panel4:
    def __init__(self,sigma=1.0,a=0.950,h = 29):
        self.sigma = sigma
        self.a = a 
        self.h = h
    def gama(self,h,sigma):
            
        return (h -2*sigma)/(h-sigma)
    
    def eps(self,h,sigma):
        
        return (h -sigma)/(sigma)
    def coef0(self,b):
            return (-(1/(self.eps(self.h,self.sigma)*np.sqrt(1+b)*b))  *( ( 1 - np.sqrt(1+b) ) * (352 - 96*self.a ) + (784 -336*self.a ) * b ) + (315 * self.gama(self.h,self.sigma) +923 *(1/self.eps(self.h,self.sigma))  ) ) 
    def trans(self,b):
        return ( (self.coef0(b) + 315*self.gama(self.h,self.sigma)+923*(1/self.eps(self.h,self.sigma)) - self.a* (105*self.gama(self.h,self.sigma) +393/self.eps(self.h,self.sigma)) - 8 * b * (35*self.gama(self.h,self.sigma) +15/self.eps(self.h,self.sigma) - self.a*( 70 *self.gama(self.h,self.sigma) - 34 / self.eps(self.h,self.sigma)  ) ) - b**2 * ( 1 + self.a ) * ( 144*self.gama(self.h,self.sigma) + 44 / self.eps(self.h,self.sigma)) ) / (105 * ( 3-self.a + b* ( 10 - 6 * self.a))) - np.arcsinh(np.sqrt(b))/(np.sqrt(b*(1+b))) ) * (np.sqrt(1+b))/(b)
    def trans1(self,alpha,b):
        return  (self.coef0(b) + 315*self.gama(self.h,self.sigma)+923*(1/self.eps(self.h,self.sigma)) - alpha* (105*self.gama(self.h,self.sigma) +393/self.eps(self.h,self.sigma)) - 8 * b * (35*self.gama(self.h,self.sigma) +15/self.eps(self.h,self.sigma) - alpha*( 70 *self.gama(self.h,self.sigma) - 34 / self.eps(self.h,self.sigma)  ) ) - b**2 * ( 1 +alpha ) * ( 144*self.gama(self.h,self.sigma) + 44 / self.eps(self.h,self.sigma)) ) / (105 * ( 3-alpha + b* ( 10 - 6 * alpha))) 
    def trans3(self,b):

        return     (2/(3360 *  b**2 *(3 - self.a + 10  * b - 6 * self.a * b)))*(-784 *  b + 352 * (-1 + np.sqrt(1 + b)) + 
                        b * np.sqrt(1 + b) *(923 + 4 *  b *  (-30 + 11 * b)) + 
                         self.a *  (-96 *  (-1 + np.sqrt(1 + b)) + 
                          b * (336 - 393 * np.sqrt(1+b) + 4 *  b *  np.sqrt(1+b) * (68 + 11 * b))) )  -  np.arcsinh(np.sqrt(b))/(48*b*np.sqrt(b))
        
        #  self.gama(self.h,self.sigma) * (np.sqrt(b+1)/(b *48))*( (4 * ( 1 + self.a ) * b ** 2 + ( 16 * self.a - 8 ) * b + ( 9 - 3 * self.a ) )/( ( 30 - 18 * self.a) * b + ( 9 -3 * self.a ) )- np.arcsinh(np.sqrt(b))/np.sqrt(b*(b+1)) )+
    def trans4(self,b):
    
        return    ( np.sqrt((b+1))*(4 * ( 1 + self.a ) * b ** 2 + ( 16 * self.a - 8 ) * b + ( 9 - 3 * self.a ) )/((48 *b )*( ( 30 - 18 * self.a) * b + ( 9 -3 * self.a ) ) ) )* self.gama(self.h,self.sigma)+(1/(1680 *  b**2 *( ( 30 - 18 * self.a) * b + ( 9 -3 * self.a ) ) * self.eps(self.h,self.sigma) ))*(-784 *  b + 352 * (-1 + np.sqrt(1 + b)) + 
                        b * np.sqrt(1 + b) *(923 + 4 *  b *  (-30 + 11 * b)) + 
                         self.a *  (-96 *  (-1 + np.sqrt(1 + b)) + 
                          b * (336 - 393 * np.sqrt(1+b) + 4 *  b *  np.sqrt(1+b) * (68 + 11 * b))) ) -  np.arcsinh(np.sqrt(b))/(48*b*np.sqrt(b))
        
class Stat_Temp:
    def __init__(self,sigma=1.0,a=0.950,h = 29):
        self.sigma = sigma
        self.a = a 
        self.h = h
    
    def gama(self,h,sigma):
            
        return (h -2*sigma)/(h-sigma)
    
    def eps(self,h,sigma):
        
        return (h -sigma)/(sigma)
    
    def integral_z1(self,b):
        
        return  ( 1/(48 * np.sqrt(b**4* (1 + b))) ) * (b *  (1 + b)  * (-9 + self.a * (3 + 2 *  b) *  (1 + 4 *  b) - 
      2 * b *  (11 + 8 *  b)) + (9 *  np.sqrt(b *  (1 + b)) - 3 *  self.a *  np.sqrt(b * (1 + b)) + 
      12 *  np.sqrt(b**3 *  (1 + b))) *  np.arcsinh(np.sqrt(b)))
      
    def integral_z2(self,b):
        
        return  (1/(1680 *  b** 2)) * (448 *  b - 352 *  (-1 + np.sqrt(1 + b)) - 
                b  * np.sqrt(1 + b) *  (587 + 2 *  b *  (129 + 64 *  b)) + 
                self.a *  (b  *  np.sqrt(1 + b) *  (3 + 4  * b) *  (19 + 10 *  b) + 96 *  (-1 + np.sqrt(1 + b))) + 
                105 * np.sqrt(b)*  (3 - self.a + 4 *  b) *  np.arcsinh(np.sqrt(b)))
        
    def  ts_w_uds(self,b,dens):   
    
        return  ( ( 1 + b )  / (   np.sqrt( np.pi) * ( 1 + self.a ) * dens * ( self.h - self.sigma ) * self.sigma ** 2  * (  self.integral_z1(b) * self.gama(self.h,self.sigma) +self.integral_z2(b)/self.eps(self.h,self.sigma)    )  )   ) ** 2

class Evol_eq:
    def __init__(self,sigma=1.0,rho = 0.02,w=0.001,alfa=0.950,altura = 29,d=3):
        self.sigma = sigma
        self.w = w
        self.alfa = alfa
        self.rho = rho 
        self.altura = altura
        self.d = d
    def gama(self,altura,sigma):
            
        return (altura -2*sigma)/(altura-sigma)
    
    def eps(self,altura,sigma):
        
        return (altura -sigma)/(sigma)
    def gg(self,d):
        
        return -( 3 + d  )/ ( 3 + 4 * d + d ** 2 ) + ( ( 2 + d ) * np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 4 * special.gamma(2 + d / 2 ) )
    
    def hh(self,d):
        
            return - 2  / ( 3 + 4 * d + d ** 2 ) + ( np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 4 * special.gamma(2 + d / 2 ) )
            
    def ggp(self,d):
    
            return - ( 6 + 2 * d  )  / ( ( 3 + d) * (- 1 + d ** 2 ) ) - (   np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 2 * (1 - d )  * special.gamma(1 + d / 2 ) )
    
    def hhp(self,d):
        
            return - 8 / ( ( 3 + d) * (- 1 + d ** 2 ) ) - ( 3 * np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 4 * (1 - d )  * special.gamma(2 + d / 2 ) )
    
    def factor1(self,d):
        
        return np.sqrt( np.pi ) * special.gamma( ( 1 + d ) / 2 ) / ( 4 * special.gamma( 2 + d / 2  )   )
    
    def factor2(self,d):
            
        return np.sqrt( np.pi ) * special.gamma( ( d - 1 ) / 2 ) / ( 8 * special.gamma( 2 + d / 2  )   )
        
   
    def solid_angle_surface(self,d):
        
        return 2 * ( np.sqrt( np.pi  ) ) ** (d-1) / special.gamma( (d-1) / 2 ) 
    
    def fpd_dim(self,t1,t2,t):
            
        return    (2/((self.d-1)*np.sqrt(np.pi))) * (1+self.alfa)*self.rho*np.sqrt(t1)*self.solid_angle_surface(self.d)*self.factor1(self.d)*( (  -(1-self.alfa) * ((self.d +2) *   self.gama(self.altura,self.sigma) + self.gg(self.d)/(self.eps(self.altura,self.sigma)*self.factor1(self.d)) )*t1 +0.5 * (1 + 3 *self.alfa) *  ( self.gama(self.altura,self.sigma) + self.hh(self.d)  / (self.factor1(self.d)* self.eps(self.altura,self.sigma)) )*(t2 -t1 )) )
 

    def gpd_dim(self,t1,t2,t):
    
        return  (2/(np.sqrt(np.pi)))  * (1+self.alfa)*self.rho*np.sqrt(t1)*self.solid_angle_surface(self.d)*self.factor2(self.d)*(   -(1-self.alfa) *   ( (self.d +2) * self.gama(self.altura,self.sigma) + self.ggp(self.d) /(self.factor2(self.d) * self.eps(self.altura,self.sigma))  ) * t1 +0.5 *(  -((5+4*self.d-9*self.alfa)) +    ((1+3*self.alfa)*self.hhp(self.d) - 4 *self.ggp(self.d))/(self.factor2(self.d)*self.eps(self.altura,self.sigma) )   )* ( t2 - t1 )  )  +( 4.0*self.w*t2/(self.altura-self.sigma) ) 
    
    
    
    def fp(self,t1,t2,t):
            
        return (1+self.alfa)*self.rho*np.sqrt(t1*np.pi)*( (  -(1-self.alfa) * (4 *   self.gama(self.altura,self.sigma) + 5/(2*self.eps(self.altura,self.sigma)) )*t1 +(1 + 3 *self.alfa) *  ( 2 *self.gama(self.altura,self.sigma)/5 + 3  / (20 * self.eps(self.altura,self.sigma)) )*(t2 -t1 )) )/3
 

    def gp(self,t1,t2,t):
    
        return   (1+self.alfa)*self.rho*np.sqrt(t1*np.pi)*( - (4 * ( 1 - self.alfa ) * self.gama(self.altura,self.sigma) + 3*(3- self.alfa) /(5 * self.eps(self.altura,self.sigma))  ) * t1 -(2*( 17 - 9 * self.alfa )* self.gama(self.altura,self.sigma)/5 + (1-self.alfa) / self.eps(self.altura,self.sigma) )* ( t2 - t1 )  ) / 3 +( 4.0*self.w*t2/(self.altura-self.sigma) ) 


class hom_eigen:
    def __init__(self,sigma=1.0,alfa=0.950, h  = 29,d = 3):
        self.sigma = sigma    
        self.h = h
        self.alfa = alfa
        self.d = d
    def gama(self,h,sigma):
            
        return (h -2*sigma)/(h-sigma)
    
    def eps_inv(self,h,sigma):
        
        return (sigma)/(h -sigma)
    
  
    def factor1(self,d):
            
        return np.sqrt( np.pi ) * special.gamma( ( 1 + d ) / 2 ) / ( 4 * special.gamma( 2 + d / 2  )   )
    
    def factor2(self,d):
            
        return np.sqrt( np.pi ) * special.gamma( ( d - 1 ) / 2 ) / ( 8 * special.gamma( 2 + d / 2  )   )
    
    
    def factor(self,d):
            
        return np.sqrt( np.pi ) * special.gamma( ( d - 1 ) / 2 ) / ( 2 * special.gamma( 2 + d / 2  )   )
    
    
    
    def gg(self,d):
        
        return (-( 3 + d  )/ ( 3 + 4 * d + d ** 2 ) + ( ( 2 + d ) * np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 4 * special.gamma(2 + d / 2 )  ) )/self.factor1(d)
    
    def hh(self,d):
        
            return (- 2  / ( 3 + 4 * d + d ** 2 ) + ( np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 4 * special.gamma(2 + d / 2 ) )) /self.factor1(d)
            
    def ggp(self,d):
    
            return (- ( 6 + 2 * d  )  / ( ( 3 + d) * (- 1 + d ** 2 ) ) - (   np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 2 * (1 - d )  * special.gamma(1 + d / 2 ) ) )  /self.factor2(d)
    
    def hhp(self,d):
        
            return (- 8 / ( ( 3 + d) * (- 1 + d ** 2 ) ) - ( 3 * np.sqrt(np.pi) *special.gamma( ( d + 1 )/2 )  )/ ( 4 * (1 - d )  * special.gamma(2 + d / 2 ) ) ) /self.factor2(d)
    
    def betas_d(self, alfa,d ):
            
        return ( 2 * ( 1 - alfa ) / (1 + 3 * alfa )  ) * ( ( 2 + d ) * self.gama(self.h,self.sigma) + self.gg(d) * self.eps_inv(self.h,self.sigma) )/( self.gama(self.h,self.sigma) + self.hh(d) * self.eps_inv(self.h,self.sigma)  ) 
    
    def  betas(self,alfa):   
    
        return (10*(1-alfa)/(1+3*alfa))*(5*(self.eps_inv(self.h,self.sigma))+8*self.gama(self.h,self.sigma))/(3*self.eps_inv(self.h,self.sigma)+self.gama(self.h,self.sigma)*8)
    
    def A(self,alfa): 
        
        return  (( 2 + self.d ) * self.gama(self.h,self.sigma)  +  self.gg(self.d) * self.eps_inv(self.h,self.sigma))
    
    def B(self,alfa): 
        
        return (1+3 * alfa  )* (self.gama(self.h,self.sigma) +self.hh(self.d) * self.eps_inv(self.h,self.sigma))
      
    def Ap(self,alfa): 
        
        return   (( 2 + self.d ) * self.gama(self.h,self.sigma)  +  self.ggp(self.d) * self.eps_inv(self.h,self.sigma))
    
    def Bp(self,alfa): 
        
        return (5 + 4 *self.d -9* alfa ) * self.gama(self.h,self.sigma) -((1+3*alfa)*self.hhp(self.d)- 4 * self.ggp(self.d)  ) * self.eps_inv(self.h,self.sigma)
      
    def a11 (self,alfa): 
        
        return (1-alfa) * self.A( alfa ) + 0.5 * self.B(alfa )
    
    def a12 (self,alfa): 
        
        return  -0.5 * self.B(alfa )
    def a21 (self,alfa): 
        
        # return    3*(1-alfa) * self.Ap( alfa ) / 2 +0.25*(self.betas_d(alfa,self.d) -2 ) * self.Bp(alfa )
        return 3*(1-alfa) * self.Ap( alfa ) / 2 +0.25*(self.betas(alfa) -2 ) * self.Bp(alfa )
    def a22 (self,alfa): 
        
        return   (0.5* self.Bp(alfa ) - (1-alfa) * self.Ap( alfa )) /(1 + self.betas_d( alfa,self.d ) )
class hom_eigen3d:
    def __init__(self,sigma=1.0,alfa=0.950, h  = 29):
        self.sigma = sigma    
        self.h = h
        self.alfa = alfa
   
    def gama(self,h,sigma):
            
        return (h -2*sigma)/(h-sigma)
    
    def eps_inv(self,h,sigma):
        
        return (sigma)/(h -sigma)

    def  betas(self,alfa):   
    
        return (10*(1-alfa)/(1+3*alfa))*(5*(self.eps_inv(self.h,self.sigma))+8*self.gama(self.h,self.sigma))/(3*self.eps_inv(self.h,self.sigma)+self.gama(self.h,self.sigma)*8)
    
  
    def A_3d(self,alfa):
        
        return 5  * self.gama(self.h,self.sigma) + 25 * self.eps_inv(self.h,self.sigma) / 8  
    
    def B_3d(self,alfa):
        
        return (1+3* alfa)*(   self.gama(self.h,self.sigma) + 3 * self.eps_inv(self.h,self.sigma) / 8   )
    
    def Ap_3d(self,alfa):
            
        return 5  * self.gama(self.h,self.sigma) + 5 * self.eps_inv(self.h,self.sigma) / 4  
    def Bp_3d(self,alfa):
        
        return (  ( 17 - 9 * alfa ) * self.gama(self.h,self.sigma) + 3 * 0.5 *( 3 - alfa )* self.eps_inv(self.h,self.sigma)     )
    
    def a11_3d (self,alfa): 
        
        return -(4 / 15) *  ( (1 - alfa ) * self.A_3d(alfa) +0.5 * self.B_3d(alfa)  ) 
        # return  - ((1-alfa) * (4 * self.gama(self.h,self.sigma) + 5 *0.5 * self.eps_inv(self.h,self.sigma) )  + (1 + 3 * alfa ) * ( 4 * self.gama(self.h,self.sigma) + 3 *0.5 * self.eps_inv(self.h,self.sigma) )/10  )/ 3
    
    def a12_3d (self,alfa): 
        
          return -(4 / 15) *  ( - 0.5 * self.B_3d(alfa)  ) 
  
    def a21_3d (self,alfa): 
        
        return   -(4 / 15) *  (  3* 0.5 *  (1 - alfa ) * self.Ap_3d(alfa) + 0.25 *(self.betas(alfa)-2) *  self.Bp_3d(alfa) ) 
    
    def a22_3d (self,alfa): 
        
         return    -(4 / 15) *  (  - (1 - alfa ) * self.Ap_3d(alfa) +  self.Bp_3d(alfa) * 0.5  )  / ( 1 + self.betas( alfa) )