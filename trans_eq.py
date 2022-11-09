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
import random
# import scipy.special as special

# Import math Library
import math 

from sympy import *

import csv

ba= pd.read_csv("b_a.txt",header=None,sep='\s+',names=["alpha","beta"])
ba_sim = pd.read_csv("dataBeta/beta_H29.dat",header=0,sep='\s+')



betas = np.linspace(0.01,10.0,1000)
alpha = np.linspace(0.010,1.000,1000)


def polinomio(a,b):
    return np.sqrt(b+1)*(4 * ( 1 + a ) * b ** 2 + ( 16 * a - 8 ) * b + ( 9 - 3 * a ) )/( ( 30 - 18 * a) * b + ( 9 -3 * a ) )
def  bets_1(a):   
        return  10*(1-a)/(1+3*a)
        # return 32*(1-a)/(5+9*a)
def  bets_2(a):   
        # return  10*(1-a)/(1+3*a)
        return 32*(1-a)/(5+9*a)
    


    

    
# fig1=plt.figure()


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



plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

plt.ylabel ( r'   ',rotation=0.0,fontsize=30)

plt.xlabel( r' $\beta_s$ ', fontsize=30)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.title ( r' \textbf {Ecuación trascendente para $\beta$}' ,fontsize=40)
plt.legend(loc=0,fontsize=30)



fig2=plt.figure()

plt.plot(alpha,bets_1(alpha) ,linewidth=1.5,color="C1", label= r"$\beta_s$  $\mathcal{O}(\beta^2)$ integral aprox." )

plt.plot(alpha,bets_2(alpha) ,linewidth=1.5,color="C2",label= r"$\beta_s$  $\mathcal{O}(\beta^2)$ integral exacta" )

plt.plot(ba['alpha'],ba['beta'],"o",color= "C0")

plt.errorbar(ba_sim['alpha'], ba_sim['beta'], yerr=ba_sim['errorbeta'], color='C2',marker="o",linestyle="",label=r"$\beta$ (MD)") 

plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

plt.ylabel ( r'  $\beta_s$  ',rotation=0.0,fontsize=30)

plt.xlabel( r' $\alpha$ ', fontsize=30)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(loc=0,fontsize=30)


# plt.xlim(0.60,1.0)

plt.title ( r' \textbf {Relación $\beta \leftrightarrow \alpha$}' ,fontsize=40)








plt.show()