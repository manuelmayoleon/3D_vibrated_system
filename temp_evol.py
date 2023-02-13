from ctypes import pointer
import numpy as np
import pandas as pd
from numba import jit
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
import functions as func
import random
import os
from scipy.integrate import solve_ivp 


# ! modificar el path si cambias el numero de particulas, la densidad y el alpha 
path = os.getcwd()
path2 = path+'/MeanTemp_vs_time/'




temp= pd.read_csv(path2+"ave_alpha_p40_h29.dat",header=0,sep='\t\t+' )



temp.columns = temp.columns.str.strip()

print(temp.columns)

n=512

sigma=1.0

# altura = 29
altura = 8 

alfa = 0.400

w=0.001

rho = 0.02



txy_ini = (temp["Tx"].iloc[0] +temp["Ty"].iloc[0])
tz_ini = 2 * temp["Tz"].iloc[0]



txy =    (temp["Tx"] +temp["Ty"] )





tz =  2 *  temp["Tz"]

time = temp['# time']





fig22 = plt.subplots(figsize=(10,8))
timet0 =time*np.sqrt(txy_ini)
# #?? Representacion en funcion del tiempo

plt.plot(timet0,txy,color='C2',marker="o",linestyle="--",label="$T_{xy}$ (MD)")

plt.plot(timet0,tz,color='C3',marker="o" ,linestyle="--",label="$T_{z}$ (MD)")


# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r' $t(T_0/m\sigma^2)^{1/2}$ ', fontsize=40)
# plt.ylabel ( r' $\frac{T}{T_0}$ ',rotation=0.0,labelpad=30,fontsize=40)

plt.xlim(17000,26000)


plt.ylim(0.2e-5,0.8e-5)


plt.tick_params(axis='x', labelsize=25)
plt.tick_params(axis='y', labelsize=25)

plt.legend(loc=0,fontsize=22)


plt.tight_layout()



plt.show()
