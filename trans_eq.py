# trans_eq
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


ba= pd.read_csv("b_a.txt",header=None,sep='\s+',names=["alpha","beta"])
ba_sim8 = pd.read_csv("dataBeta/beta_H8.dat",header=0,sep='\s+')
ba_sim10 = pd.read_csv("dataBeta/beta_H10.dat",header=0,sep='\s+')
ba_sim12 = pd.read_csv("dataBeta/beta_H12.dat",header=0,sep='\s+')
ba_sim14 = pd.read_csv("dataBeta/beta_H14.dat",header=0,sep='\s+')
ba_sim29 = pd.read_csv("dataBeta/beta_H29.dat",header=0,sep='\s+')


betas = np.linspace(0.01,10.0,1000)
alpha = np.linspace(0.010,1.000,1000)

sigma = 1.0
rho = 0.02*sigma

class Panel:
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
          return (1+alfa)*(Panel.b1(alfa,b) + Panel.b2(alfa,b))/( 96 * b ) 
    def btilde(alfa,b): 
          return (1+alfa)*(Panel.b1(alfa,b) + Panel.b2(alfa,b))/( 96 * b )
    def b_coef(self,alfa,b,dens,h):
        return (1+b)/(np.sqrt(2*np.pi)*Panel.btilde(alfa,b)*Panel().htilde(dens,h))
    def v1s_w(self,alfa,b,dens,h):
        return 0.5 * Panel().b_coef(alfa,b,dens,h)*( 1 + np.sqrt(1 + 8 / ( Panel().b_coef(alfa,b,dens,h)*np.sqrt((1+b)*np.pi) ) ) )
    
# ! SOLVE TRANSCENDENTAL EQUATION USING FSOLVE USING Powell's dog leg method
# ?? Powell's dog leg method is an iterative optimisation algorithm for the solution of non-linear least squares problems. 
# ??       it combines the Gauss–Newton algorithm with gradient descent, but it uses an explicit trust region. 

sol = np.zeros(len(alpha))    
v1sw_h10 = np.zeros(len(alpha))
v1sw_h12 = np.zeros(len(alpha))
v1sw_h14 = np.zeros(len(alpha))
v1sw_h29 = np.zeros(len(alpha))
i=0
for x in alpha:
    

    sol[i] = fsolve(Panel(x).transcendental_eq ,[4.0] )
    
    v1sw_h10[i] = Panel().v1s_w(x,sol[i],rho,10)
    v1sw_h12[i] = Panel().v1s_w(x,sol[i],rho,12)
    v1sw_h14[i] = Panel().v1s_w(x,sol[i],rho,14)
    v1sw_h29[i] = Panel().v1s_w(x,sol[i],rho,29)
   
    # print(alpha[i]) 
    i+=1
        
    # print(result.x)
# print(sol)




fig22 = plt.plot(figsize=(14,8))
plt.plot(alpha[0:len(alpha)-1],v1sw_h10[0:len(alpha)-1],color= "C0")
plt.plot(alpha[0:len(alpha)-1],v1sw_h12[0:len(alpha)-1],color= "C1")
plt.plot(alpha[0:len(alpha)-1],v1sw_h14[0:len(alpha)-1],color= "C2")
plt.plot(alpha[0:len(alpha)-1],v1sw_h29[0:len(alpha)-1],color= "C3")
# plt.plot(betas,Panel.b1(0.99,betas),color= "C1")
# plt.plot(betas,Panel.b2(0.99,betas),color= "C2")
plt.xlim(0.75,1.0)
plt.ylim(0,1000)
# print(Panel.b1(0.99 ,betas))

# ! SOLVE TRANSCENDENTAL EQUATION USING graphical method 
# plt.plot(betas,np.arcsinh(np.sqrt(betas))/np.sqrt(betas),color= "C0")
 
# for x in alpha:

#     r = random.random()
#     b =      random.random()
#     g = random.random()

#     color = (r, g, b)
    
#     plt.plot(betas,polinomio(x,betas) ,linestyle = ":",linewidth=1.5,color=color,label= r"$\alpha = $ %2.2f" % x)
    # for y in betas :
    #     if (abs(polinomio(x,y)-np.arcsinh(np.sqrt(y))/np.sqrt(y)) <= 0.001):
    #         # print('alpha')
    #         # print(x)
    #         # print('beta')
    #         # print( y)
    #         with open('b_a.txt', 'a',newline='\n') as f:
    #                 writer = csv.writer(f, delimiter='\t')
    #                 writer.writerows(zip([x],[y]))

# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

# plt.ylabel ( r'   ',rotation=0.0,fontsize=30)

# plt.xlabel( r' $\beta_s$ ', fontsize=30)

# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.title ( r' \textbf {Ecuación trascendente para $\beta$}' ,fontsize=40)
# plt.legend(loc=0,fontsize=30)



fig, ax1=plt.subplots(figsize=(14,8))

ax1.plot(alpha,Panel.bets_1(alpha) ,linewidth=1.5,alpha=0.6,color="C1" )

ax1.plot(alpha,Panel.bets_2(alpha) ,linewidth=1.5,alpha =0.6,color="C2")
ax1.plot(alpha,sol ,linewidth=1.5,linestyle= "dashdot",color='xkcd:crimson')

# plt.plot(betas,np.arcsinh(np.sqrt(betas))/np.sqrt(betas),color= "C0")

ax1.plot(ba['alpha'],ba['beta'],linestyle = ":",color= "C0")

ax1.errorbar(ba_sim8['alpha'], ba_sim8['beta'], yerr=ba_sim8['errorbeta'], color='C3',marker="o",linestyle="") 
ax1.errorbar(ba_sim10['alpha'], ba_sim10['beta'], yerr=ba_sim10['errorbeta'], color='C4',marker="o",linestyle="") 
ax1.errorbar(ba_sim12['alpha'], ba_sim12['beta'], yerr=ba_sim12['errorbeta'], color='C5',marker="o",linestyle="") 
ax1.errorbar(ba_sim14['alpha'], ba_sim14['beta'], yerr=ba_sim14['errorbeta'], color='C6',marker="o",linestyle="") 
ax1.errorbar(ba_sim29['alpha'], ba_sim29['beta'], yerr=ba_sim29['errorbeta'], color='C7',marker="o",linestyle="") 

ax1.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

ax1.set_ylabel( r'  $\beta_s$  ',rotation=0.0,fontsize=30)

ax1.set_xlabel( r' $\alpha$ ', fontsize=30)

ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
# ax1.legend(loc=0,fontsize=30)


# plt.xlim(0.60,1.0)

ax1.set_title( r' \textbf {Relación $\beta \leftrightarrow \alpha$}' ,fontsize=40)


# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
ax2 = plt.axes([0,0,1,1])


# Manually set the position and relative size of the inset axes within ax1


ip = InsetPosition(ax1, [0.45,0.45,0.5,0.5])


ax2.set_axes_locator(ip)

ax2.plot(alpha,Panel.bets_1(alpha) ,linewidth=1.5,color="C1",alpha=0.6, label= r"$\beta_s$  $\mathcal{O}(\beta^2)$ integral aprox." )

ax2.plot(alpha,Panel.bets_2(alpha) ,linewidth=1.5,color="C2",alpha =0.6,label= r"$\beta_s$  $\mathcal{O}(\beta^2)$ integral exacta" )

ax2.plot(ba['alpha'],ba['beta'],linewidth=1.5,linestyle = ":",color= "C0",label = r"$\beta_s$  exacta (sol. gráfica)")

ax2.plot(alpha,sol ,linewidth=1,linestyle= "dashdot",color='xkcd:crimson',label = r"$\beta_s$  integral exacta" )

ax2.errorbar(ba_sim8['alpha'], ba_sim8['beta'], yerr=ba_sim8['errorbeta'], color='C3',marker="o",linestyle="",label=r"$\beta_s$, $h= 8\sigma$ (MD)") 
ax2.errorbar(ba_sim10['alpha'], ba_sim10['beta'], yerr=ba_sim10['errorbeta'], color='C4',marker="o",linestyle="",label=r"$\beta_s$, $h= 10\sigma$ (MD)") 
ax2.errorbar(ba_sim12['alpha'], ba_sim12['beta'], yerr=ba_sim12['errorbeta'], color='C5',marker="o",linestyle="",label=r"$\beta_s$, $h= 12\sigma$ (MD)") 
ax2.errorbar(ba_sim14['alpha'], ba_sim14['beta'], yerr=ba_sim14['errorbeta'], color='C6',marker="o",linestyle="",label=r"$\beta_s$, $h= 14\sigma$ (MD)") 
ax2.errorbar(ba_sim29['alpha'], ba_sim29['beta'], yerr=ba_sim29['errorbeta'], color='C7',marker="o",linestyle="",label=r"$\beta_s$, $h= 29\sigma$ (MD)") 

ax2.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
ax2.legend(loc=0)

# Some ad hoc tweaks.

ax2.set_xlim(0.77,1.0)
ax2.set_ylim(0.0,0.77)




plt.show()