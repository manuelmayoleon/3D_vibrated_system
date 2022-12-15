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

import functions as f


import os
import glob
  


save_temp_figures = True
save_beta_figure = True 
  
path = os.getcwd()
path = path+"/Txymean_vs_alpha/"
csv_files = glob.glob(os.path.join(path, "*.dat"))
csv_files = sorted(csv_files, key=len)
print(len(csv_files))

names = list()
elements=list()
# loop over the list of csv files
for ii in csv_files:
        
        # read the csv file
        # df= pd.read_csv(f)
       
        df= pd.read_csv(ii,header=0,sep='\t\t+',engine='python')
        elements.append(df)
        # print the location, filename and name
        print('Location:', ii)
        file_name = ii.split("/")[-1]
        print('File Name:', file_name )
        name = file_name.split(".")[0]
        print('Name:', name)
        names.append(name)
          
        # print the content
        print('Content:')
        print(df)  
        
H8 =  elements[0]  
H10 = elements[-1]   
H12 = elements[2]
H14 = elements[1]
H29 = elements[3]

print(H29)

# H8.columns = H8.iloc[0]
# ! Para quitar los espacios al principio y al final de la cabecera
H8.columns = H8.columns.str.strip()
H8 = H8.rename(columns={'# alpha': 'alpha'})
H10.columns = H10.columns.str.strip()
H10 = H10.rename(columns={'# alpha': 'alpha'})
H12.columns = H12.columns.str.strip()
H12 = H12.rename(columns={'# alpha': 'alpha'})
H14.columns = H14.columns.str.strip()
H14 = H14.rename(columns={'# alpha': 'alpha'})
H29.columns = H29.columns.str.strip()
H29 = H29.rename(columns={'# alpha': 'alpha'})
# print( H8.columns )

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
w = 0.001



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



    
# ! SOLVE TRANSCENDENTAL EQUATION USING FSOLVE. USING Powell's dog leg method
# ?? Powell's dog leg method is an iterative optimisation algorithm for the solution of non-linear least squares problems. 
# ??       it combines the Gauss–Newton algorithm with gradient descent, but it uses an explicit trust region. 

sol = np.zeros(len(alpha))    
v1sw_h8 = np.zeros(len(alpha))
v1sw_h10 = np.zeros(len(alpha))
v1sw_h12 = np.zeros(len(alpha))
v1sw_h14 = np.zeros(len(alpha))
v1sw_h29 = np.zeros(len(alpha))
i=0
for x in alpha:
    

    sol[i] = fsolve(f.Panel1(x).transcendental_eq ,[4.0] )
    v1sw_h8[i] = f.Panel1().v1s_w(x,sol[i],rho,8)
    v1sw_h10[i] = f.Panel1().v1s_w(x,sol[i],rho,10)
    v1sw_h12[i] = f.Panel1().v1s_w(x,sol[i],rho,12)
    v1sw_h14[i] = f.Panel1().v1s_w(x,sol[i],rho,14)
    v1sw_h29[i] = f.Panel1().v1s_w(x,sol[i],rho,29)
   
    # print(alpha[i]) 
    i+=1
        
    # print(result.x)
# print(sol)







# !! FIGURA T_xy <--> alpha

fig22 = plt.figure(figsize=(14,8))

# ! h =8 sigma

plt.plot(alpha[0:len(alpha)-1],v1sw_h8[0:len(alpha)-1]**2/2,color= "C0")

plt.errorbar(H8['alpha'], 2*H8['Txymean'], yerr=H8['Txy_error'],mfc="none",capsize=10,ms=12, color='C0',marker="o",linestyle="",label= r" $H =  %2i  \sigma $ "   % 8) 

plt.plot(alpha,(f.Panel2(1.0,w,1.0,8).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C0")

# ! h =10 sigma

plt.plot(alpha[0:len(alpha)-1],v1sw_h10[0:len(alpha)-1]**2/2,color= "C1")

plt.errorbar(H10['alpha'], 2*H10['Txymean'], yerr=H10['Txy_error'],mfc="none",capsize=10,ms=12, color='C1',marker="s",linestyle="",label= r" $H =  %2i  \sigma $ "   % 10) 

plt.plot(alpha,(f.Panel2(1.0,w,1.0,10).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C1")

# ! h =12 sigma

plt.plot(alpha[0:len(alpha)-1],v1sw_h12[0:len(alpha)-1]**2/2,color= "C2")

plt.errorbar(H12['alpha'], 2*H12['Txymean'], yerr=H12['Txy_error'],mfc="none",capsize=10,ms=12, color='C2',marker="d",linestyle="",label= r" $H =  %2i  \sigma $ "   % 12) 

plt.plot(alpha,(f.Panel2(1.0,w,1.0,12).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C2")

# ! h =14 sigma

plt.plot(alpha[0:len(alpha)-1],v1sw_h14[0:len(alpha)-1]**2/2,color= "C3")

plt.errorbar(H14['alpha'], 2*H14['Txymean'], yerr=H14['Txy_error'],mfc="none",capsize=10,ms=12, color='C3',marker="*",linestyle="",label= r" $H =  %2i  \sigma $ "   % 14) 


plt.plot(alpha,(f.Panel2(1.0,w,1.0,14).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C3")

# ! h =29 sigma

plt.plot(alpha[0:len(alpha)-1],v1sw_h29[0:len(alpha)-1]**2/2,color= "C4")

plt.errorbar(H29['alpha'], 2*H29['Txymean'], yerr=H29['Txy_error'],mfc="none",capsize=10,ms=12, color='C4',marker="^",linestyle="",label= r" $H =  %2i  \sigma $ "   % 29) 

plt.plot(alpha,(f.Panel2(1.0,w,1.0,29).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C4")


plt.ylabel( r'  $\frac{T^s}{m v_p^2}$  ',rotation=0.0,labelpad=30,fontsize=40)

plt.xlabel( r' $\alpha$ ',fontsize=40)


# plt.title( r' \textbf {Relación $T_{xy} \leftrightarrow \alpha$}' ,fontsize=40)

plt.xlim(0.55,0.99)
plt.ylim(0,2000)


plt.tick_params(axis='x', labelsize=25)
plt.tick_params(axis='y', labelsize=25)

# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.legend(loc=0,fontsize=22)

plt.tight_layout()

if (save_temp_figures == True) :

    plt.savefig('txy_vs_alpha.pdf', dpi = 1200)
# plt.savefig('txy_vs_alpha.eps', dpi = 500, format="eps")

# !! FIGURA T_z <--> alpha

v2sw_h8 = v1sw_h8**2/2*(1+f.Panel1.bets_2(alpha))
v2sw_h10 = v1sw_h10**2/2*(1+f.Panel1.bets_2(alpha))
v2sw_h12 = v1sw_h12**2/2*(1+f.Panel1.bets_2(alpha))
v2sw_h14 = v1sw_h14**2/2*(1+f.Panel1.bets_2(alpha))
v2sw_h29 = v1sw_h29**2/2*(1+f.Panel1.bets_2(alpha))

fig23 = plt.figure(figsize=(14,8))

# ! h = 8 sigma


plt.plot(alpha[0:len(alpha)-1],v2sw_h8[0:len(alpha)-1],color= "C0")

plt.errorbar(H8['alpha'], 2*H8['Tzmean'], yerr=H8['Tz_error'],mfc="none",capsize=10,ms=12, color='C0',marker="o",linestyle="",label= r" $H =  %2i  \sigma $ "   % 8) 

plt.plot(alpha,(1+f.Panel1.bets_1(alpha))*(f.Panel2(1.0,w,1.0,8).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C0")

# ! h = 10 sigma



plt.plot(alpha[0:len(alpha)-1],v2sw_h10[0:len(alpha)-1],color= "C1")

plt.errorbar(H10['alpha'], 2*H10['Tzmean'], yerr=H10['Tz_error'],mfc="none",capsize=10,ms=12, color='C1',marker="s",linestyle="",label= r" $H =  %2i  \sigma $ "   % 10) 

plt.plot(alpha,(1+f.Panel1.bets_1(alpha))*(f.Panel2(1.0,w,1.0,10).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C1")

# ! h =12 sigma



plt.plot(alpha[0:len(alpha)-1],v2sw_h12[0:len(alpha)-1],color= "C2")

plt.errorbar(H12['alpha'], 2*H12['Tzmean'], yerr=H12['Tz_error'],mfc="none",capsize=10,ms=12, color='C2',marker="d",linestyle="",label= r" $H =  %2i  \sigma $ "   % 12) 

plt.plot(alpha,(1+f.Panel1.bets_1(alpha))*(f.Panel2(1.0,w,1.0,12).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C2")

# ! h =14 sigma


plt.plot(alpha[0:len(alpha)-1],v2sw_h14[0:len(alpha)-1],color= "C3")

plt.errorbar(H14['alpha'], 2*H14['Tzmean'], yerr=H14['Tz_error'],mfc="none",capsize=10,ms=12, color='C3',marker="*",linestyle="",label= r" $H =  %2i  \sigma $ "   % 14) 

plt.plot(alpha,(1+f.Panel1.bets_1(alpha))*(f.Panel2(1.0,w,1.0,14).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C3")

# ! h =29 sigma

plt.plot(alpha[0:len(alpha)-1],v2sw_h29[0:len(alpha)-1],color= "C4")

plt.errorbar(H29['alpha'], 2*H29['Tzmean'], yerr=H29['Tz_error'],mfc="none",capsize=10,ms=12, color='C4',marker="^",linestyle="",label= r" $H =  %2i  \sigma $ "   % 29) 

plt.plot(alpha,(1+f.Panel1.bets_1(alpha))*(f.Panel2(1.0,w,1.0,29).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C4")

plt.ylabel( r'  $\frac{T^s_{z}}{m v_p^2}$  ',rotation=0.0,labelpad=30,fontsize=40)

plt.xlabel( r' $\alpha$ ',fontsize=40)

plt.tick_params(axis='x', labelsize=25)
plt.tick_params(axis='y', labelsize=25)

# plt.title( r' \textbf {Relación $T_{z} \leftrightarrow \alpha$}' ,fontsize=40)

plt.xlim(0.55,0.99)
plt.ylim(0,3000)

# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.legend(loc=0,fontsize=22)


plt.tight_layout()

if (save_temp_figures == True) :
    
    plt.savefig('tz_vs_alpha.pdf', dpi = 1200)
# plt.savefig('tz_vs_alpha.eps', dpi = 500, format="eps")




fig, ax1=plt.subplots(figsize=(14,8))


ax1.plot(alpha,f.Panel1.bets_1(alpha) ,linewidth=1.5,color="k")

# ax1.plot(alpha,f.Panel1.bets_2(alpha) ,linewidth=1.5,color="C2",alpha =0.6,label= r"$\beta_s$  $\mathcal{O}(\beta^2)$ integral exacta" )

# ax1.plot(ba['alpha'],ba['beta'],linewidth=1.5,linestyle = ":",color= "C0",label = r"$\beta_s$  exacta (sol. gráfica)")

# ax1.plot(alpha,sol ,linewidth=1,linestyle= "dashdot",color='xkcd:crimson' )

# ax1.plot(alpha,sol ,linewidth=1.5,linestyle= "dashdot",color='#00FF00' )

ax1.plot(alpha,sol ,linewidth=1.5,linestyle= "dashdot",color='b' )

# ax1.errorbar(ba_sim8['alpha'], ba_sim8['beta'], yerr=ba_sim8['errorbeta'],mfc="none",capsize=10,ms=12, color='C0',marker="o",linestyle="",label=r" $H= 8\sigma$ ") 
# ax1.errorbar(ba_sim10['alpha'], ba_sim10['beta'], yerr=ba_sim10['errorbeta'],mfc="none",capsize=10,ms=12, color='C1',marker="s",linestyle="",label=r" $H= 10\sigma$ ") 
# ax1.errorbar(ba_sim12['alpha'], ba_sim12['beta'], yerr=ba_sim12['errorbeta'],mfc="none",capsize=10,ms=12, color='C2',marker="d",linestyle="",label=r" $H= 12\sigma$ ") 
# ax1.errorbar(ba_sim14['alpha'], ba_sim14['beta'], yerr=ba_sim14['errorbeta'],mfc="none",capsize=10,ms=12, color='C3',marker="*",linestyle="",label=r"$H= 14\sigma$ ") 
# ax1.errorbar(ba_sim29['alpha'], ba_sim29['beta'], yerr=ba_sim29['errorbeta'],mfc="none",capsize=10,ms=12, color='C4',marker="^",linestyle="",label=r" $H= 29\sigma$ ") 







# ax1.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

ax1.set_ylabel( r'  $\beta_s$  ',rotation=0.0,labelpad=30,fontsize=40)

ax1.set_xlabel( r' $\alpha$ ', fontsize=40)

ax1.tick_params(axis='x', labelsize=25)
ax1.tick_params(axis='y', labelsize=25)
# ax1.legend(loc="lower left",fontsize=19)


# ax1.set_xlim(0.55,1.0)
# ax1.set_ylim(0.0,2.0)
# plt.xlim(0.60,1.0)

# ax1.set_title( r' \textbf {Relación $\beta \leftrightarrow \alpha$}' ,fontsize=40)


# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
ax2 = plt.axes([0,0,1,1])


# Manually set the position and relative size of the inset axes within ax1


ip = InsetPosition(ax1, [0.45,0.45,0.5,0.5])


ax2.set_axes_locator(ip)

ax2.plot(alpha,f.Panel1.bets_1(alpha) ,linewidth=1.5,color="k" )

# ax1.plot(alpha,f.Panel1.bets_2(alpha) ,linewidth=1.5,alpha =0.6,color="C2")
# ax2.plot(alpha,sol ,linewidth=1.5,linestyle= "dashdot",color='xkcd:crimson')
ax2.plot(alpha,sol ,linewidth=1.5,linestyle= "dashdot",color='b')

# plt.plot(betas,np.arcsinh(np.sqrt(betas))/np.sqrt(betas),color= "C0")

# ax1.plot(ba['alpha'],ba['beta'],linestyle = ":",color= "C0")

# ax2.errorbar(ba_sim8['alpha'], ba_sim8['beta'], yerr=ba_sim8['errorbeta'],mfc="none",capsize=10,ms=12, color='C0',marker="o",linestyle="") 
# ax2.errorbar(ba_sim10['alpha'], ba_sim10['beta'], yerr=ba_sim10['errorbeta'],mfc="none",capsize=10,ms=12, color='C1',marker="s",linestyle="") 
# ax2.errorbar(ba_sim12['alpha'], ba_sim12['beta'], yerr=ba_sim12['errorbeta'],mfc="none",capsize=10,ms=12, color='C2',marker="d",linestyle="") 
# ax2.errorbar(ba_sim14['alpha'], ba_sim14['beta'], yerr=ba_sim14['errorbeta'],mfc="none",capsize=10,ms=12, color='C3',marker="*",linestyle="") 
# ax2.errorbar(ba_sim29['alpha'], ba_sim29['beta'], yerr=ba_sim29['errorbeta'],mfc="none",capsize=10,ms=12, color='C4',marker="^",linestyle="") 

# ax2.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
# ax2.legend(loc=0)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
# Some ad hoc tweaks.

ax2.set_xlim(0.55,1.0)
ax2.set_ylim(0.0,2.00)


plt.tight_layout()

if (save_beta_figure == True) :
#  plt.savefig('beta_vs_alpha.pdf', dpi = 1200)
    plt.savefig('beta_vs_alpha_teo.pdf',dpi=1200)


plt.show()