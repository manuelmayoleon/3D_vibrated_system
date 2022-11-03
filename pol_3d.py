from numpy.core.function_base import linspace
import pandas as pd
import numpy as np
import scipy . stats as ss
import math
import matplotlib . mlab as mlab
import matplotlib . pyplot as plt
from scipy import optimize
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
from sympy import *
import random
# import scipy.linalg as la




# x = Symbol('x')
# A1 =Symbol('A1')
# A2 =Symbol('A2')

# # alfa = Symbol('alfa')
# # h = Symbol('h')
# bet = Symbol('bet')
# # bet =  10*(1-alfa)/(1+3*alfa) 
# #! A2 = A1 * np.sqrt(bet+1)
# A1 = (24* (np.sqrt(1+bet))/(np.sqrt(2*np.pi)*(1+alfa)*h*(2*(1-alfa)+bet*(9-17*alfa)/5)))

# polinomio = x**2 - A1*A2/2 * x - A1/np.sqrt(np.pi)


# 0.25*A1*A2 + 6.32455532033676e-8*sqrt(A1*(15625000000000.0*A1*A2**2 + 141047395886939.0))


# quadratic_equation = Eq(polinomio, 0)

# solucion=solve(quadratic_equation, x)

# print(solucion[1])


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



fig23=plt.figure()


# dens = np.linspace(1.5, 4.5, 1000)

# dens = np.linspace(0.005 , max(tvsd['den']) , 1000)
dens  = 0.02 
w = 0.001
sigma = 1.0
# h =np.linspace( 10*sigma,29*sigma,10)
h = [10.0* sigma , 12.0* sigma , 14.0* sigma , 29.0* sigma]

m = 1.0
# alpha = 0.9
alpha = np.linspace(0.60,0.98,10000)

for x in h:
    
    r = random.random()

    b =      random.random()

    g = random.random()

    color = (r, g, b)
    # plt.plot(kk,functions.eigenvalue1(functions.Panel(1,sigma,epsilon,vp).mu(x),functions.Panel(1,sigma,epsilon,vp).kapa(x),q,kk),color='C2',linestyle="--",label="$\lambda_1$ ")
    # plt.plot(kk,np.real(functions.eigenvalue1(functions.Panel(1,sigma,epsilon,vp).mu(x),functions.Panel(1,sigma,epsilon,vp).kapa(x),abs(functions.lamda1(x,epsilon)),kk)),linewidth=1.5,color=color,label=" $\lambda_1$ %1.3f" % x)
    # plt.plot(kk,np.real(functions.eigenvalue2(functions.Panel(1,sigma,epsilon,vp).mu(x),functions.Panel(1,sigma,epsilon,vp).kapa(x),abs(functions.lamda1(x,epsilon)),kk)),linewidth=1.5,linestyle="--",color=color,label=" $\lambda_2$ %1.3f" % x)
    
    # plt.plot(alpha,(Panel(1.0,w,1.0,x).vsw(dens,alpha)*w)**2,linewidth=1.5,color=color,label= " $H = $ %2.1f  $\sigma $" % x)
    
    plt.plot(alpha,m*(Panel(1.0,w,1.0,x).vsw_mano(dens,alpha))**2/2,linestyle = ":",linewidth=1.5,color=color,label= " $H = $ %2.1f  $\sigma $" % x)

    #! calcular si es cero 
    cero = Panel(1.0,w,1.0,x).vsw_mano(dens,alpha)**2 - Panel(1.0,w,1.0,x).b1(dens,alpha)*Panel(1.0,w,1.0,x).vsw_mano(dens,alpha)*np.sqrt(1+Panel(1.0,w,1.0,x).betas(alpha)) - 2*Panel(1.0,w,1.0,x).b1(dens,alpha)/np.sqrt(np.pi)
    for i in cero:
        if i >10e-12:
            print( i )
    # plt.plot(kk,functions.eigenvalue3(functions.Panel(1,sigma,epsilon,vp).mu(x),functions.Panel(1,sigma,epsilon,vp).kapa(x),functions.lamda1(x,epsilon),kk),color='C3',linestyle="--",label="$\lambda_3$ ")


# plt.plot(alpha,(Panel(1.0,w,1.0,29).vsw(dens,alpha))**2,linewidth=1.5,color='C0',label=" $\H %1.3f")

# plt.plot(alpha,(Panel(1.0,w,1.0,29).a1(dens,alpha)),linewidth=1.5,color='C0',label=" $\H %1.3f")

# plt.plot(alpha,Panel(1.0,0.001,1.0).p(dens,alpha),linewidth=1.5,linestyle=":",color="C0",label=" $p$  "  )
# plt.plot(tvsd['den'],tvsd['Ts1'],marker= 'o',linestyle = 'None',color="C2",label=" $v_{1}^s/\omega$ (MD) "  )

# plt.errorbar(tvsd['den'], tvsd['ts'], yerr=tvsd['err'], color='C0',marker="o",linestyle="",label="$n_z$ (MD)") 
# plt.plot(dens,Panel(1.0,w,1.0,29).v1s(dens,alpha)/(2*w),linewidth=1.5,color="C0",label=" $v_{1}^s/\omega$  "  )
# plt.plot(dens,Panel(1.0,w,1.0,29).vsw(dens,alpha)**2,linewidth=1.5,linestyle=":",color="C1",label=" $v_{1}^s/\omega$ $O(w^2)$  "  )

# plt.yscale("log")
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

plt.ylabel ( r'  $T_s$ ',rotation=0.0,fontsize=30)

plt.xlabel( r' $\alpha$ ', fontsize=30)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.title ( r' \textbf {Temperatura como función de $\alpha$}' ,fontsize=40)
plt.legend(loc=0,fontsize=30)

plt.show()