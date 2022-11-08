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
# import scipy.linal


def a2(a,b):
    return (1+a)*0.5*(np.sqrt(b+1)*(4*b**3+16*b**2-3*b)+np.sqrt(b)*(18*np.arcsinh(np.sqrt(b))*b+3*np.arcsinh(np.sqrt(b))))/(6)
def a1(b):
    return (np.sqrt(b+1)*(2*b**2-b)+np.sqrt(b)*(4*np.arcsinh(np.sqrt(b))*b+np.arcsinh(np.sqrt(b))))

alfa = np.linspace(0.60,0.99,10)
betas = np.linspace(0.0001, 1.000, 10000)


fig1=plt.figure()

plt.plot(betas,a1(betas),linewidth=1.5,color='C0')

for x in alfa:
    
    r = random.random()

    b =      random.random()

    g = random.random()

    color = (r, g, b)
    # plt.plot(kk,functions.eigenvalue1(functions.Panel(1,sigma,epsilon,vp).mu(x),functions.Panel(1,sigma,epsilon,vp).kapa(x),q,kk),color='C2',linestyle="--",label="$\lambda_1$ ")
    # plt.plot(kk,np.real(functions.eigenvalue1(functions.Panel(1,sigma,epsilon,vp).mu(x),functions.Panel(1,sigma,epsilon,vp).kapa(x),abs(functions.lamda1(x,epsilon)),kk)),linewidth=1.5,color=color,label=" $\lambda_1$ %1.3f" % x)
    # plt.plot(kk,np.real(functions.eigenvalue2(functions.Panel(1,sigma,epsilon,vp).mu(x),functions.Panel(1,sigma,epsilon,vp).kapa(x),abs(functions.lamda1(x,epsilon)),kk)),linewidth=1.5,linestyle="--",color=color,label=" $\lambda_2$ %1.3f" % x)
    
    # plt.plot(alpha,(Panel(1.0,w,1.0,x).vsw(dens,alpha)*w)**2,linewidth=1.5,color=color,label= " $H = $ %2.1f  $\sigma $" % x)
    
    plt.plot(betas,a2(x,betas),linestyle = ":",linewidth=1.5,color=color,label= r" $\alpha = $ %2.3f  " % x)

plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

plt.ylabel ( r'   ',rotation=0.0,fontsize=30)

plt.xlabel( r' $\beta_s$ ', fontsize=30)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.title ( r' \textbf {Ecuación trascendente}' ,fontsize=40)
plt.legend(loc=0,fontsize=30)

plt.show()