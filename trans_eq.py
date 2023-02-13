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
  


save_temp_figures = False
save_beta_figure = False
  
path = os.getcwd()
path = path+"/Txymean_vs_alpha/"
csv_files = glob.glob(os.path.join(path, "*.dat"))
csv_files = sorted(csv_files, key=len)
# print(len(csv_files))

names = list()
elements=list()
# loop over the list of csv files
for ii in csv_files:
        
        # read the csv file
        # df= pd.read_csv(f)
       
        df= pd.read_csv(ii,header=0,sep='\t\t+',engine='python')
        elements.append(df)
        # print the location, filename and name
        # print('Location:', ii)
        file_name = ii.split("/")[-1]
        # print('File Name:', file_name )
        name = file_name.split(".")[0]
        # print('Name:', name)
        names.append(name)
          
        # print the content
        # print('Content:')
        # print(df)  
        
# print(names)
H2 =  elements[names.index('beta_Txyz_vs_alphas_H2Sigma')]   
H4 =  elements[names.index('beta_Txyz_vs_alphas_H4Sigma')]   
H6 =  elements[names.index('beta_Txyz_vs_alphas_H6Sigma')]  
H8 =  elements[names.index('beta_Txyz_vs_alphas_H8Sigma')]  
H10 = elements[names.index('beta_Txyz_vs_alphas_H10Sigma')]   
H12 = elements[names.index('beta_Txyz_vs_alphas_H12Sigma')]
H14 = elements[names.index('beta_Txyz_vs_alphas_H14Sigma')]
H29 = elements[names.index('beta_Txyz_vs_alphas_H29Sigma')]

# print(H29)

# H8.columns = H8.iloc[0]
# ! Para quitar los espacios al principio y al final de la cabecera
H2.columns = H2.columns.str.strip()
H2 = H2.rename(columns={'# alpha': 'alpha'})
H4.columns = H4.columns.str.strip()
H4 = H4.rename(columns={'# alpha': 'alpha'})
H6.columns = H6.columns.str.strip()
H6 = H6.rename(columns={'# alpha': 'alpha'})
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
ba_sim2 = pd.read_csv("dataBeta/beta_H2.dat",header=0,sep='\s+')
ba_sim4 = pd.read_csv("dataBeta/beta_H4.dat",header=0,sep='\s+')
ba_sim6 = pd.read_csv("dataBeta/beta_H6.dat",header=0,sep='\s+')
ba_sim8 = pd.read_csv("dataBeta/beta_H8.dat",header=0,sep='\s+')
ba_sim10 = pd.read_csv("dataBeta/beta_H10.dat",header=0,sep='\s+')
ba_sim12 = pd.read_csv("dataBeta/beta_H12.dat",header=0,sep='\s+')
ba_sim14 = pd.read_csv("dataBeta/beta_H14.dat",header=0,sep='\s+')
ba_sim29 = pd.read_csv("dataBeta/beta_H29.dat",header=0,sep='\s+')
ba_sim29p = pd.read_csv("dataBeta/beta_H29p.dat",header=0,sep='\s+')


betas = np.linspace(0.01,10.0,1000)
alpha = np.linspace(0.200,1.00,1000)

sigma = 1.0
rho = 0.02*sigma
w = 0.001



# ! SOLVE TRANSCENDENTAL EQUATION USING graphical method 

# beta_exacta=list()
# alpha_asoc=list()
# plt.plot(betas,np.arcsinh(np.sqrt(betas))/np.sqrt(betas*(betas+1)),color= "C0")

# for x in alpha:

#     r = random.random()
#     b =      random.random()
#     g = random.random()

#     color = (r, g, b)
    
#     plt.plot(betas,f.Panel4(1.0,x,7).trans1(x,betas) ,linestyle = ":",linewidth=1.5,color=color,label= r"$\alpha = $ %2.2f" % x)
#     for y in betas :
#         if (abs(f.Panel4(1.0,x,7).trans1(x,y)-np.arcsinh(np.sqrt(y))/np.sqrt(y*(y+1))) <= 0.001):
#             # print('alpha')
#             # print(x)
#             # print('beta')
#             # print( y)
#             # beta_exacta.append(y)
#             # alpha_asoc.append(x) 
#             with open('b_a.txt', 'a',newline='\n') as ff:
#                     writer = csv.writer(ff, delimiter='\t')
#                     writer.writerows(zip([ x],[ y]))
        

# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

# plt.ylabel ( r'   ',rotation=0.0,fontsize=30)

# plt.xlabel( r' $\beta_s$ ', fontsize=30)

# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.title ( r' \textbf {Ecuación trascendente para $\beta$}' ,fontsize=40)
# # plt.legend(loc=0,fontsize=30)



    
# ! SOLVE TRANSCENDENTAL EQUATION USING FSOLVE. USING Powell's dog leg method
# ?? Powell's dog leg method is an iterative optimisation algorithm for the solution of non-linear least squares problems. 
# ??       it combines the Gauss–Newton algorithm with gradient descent, but it uses an explicit trust region. 
sol = np.zeros(len(alpha))    




v1sw_h4 = np.zeros(len(alpha))
v1sw_h8 = np.zeros(len(alpha))
v1sw_h10 = np.zeros(len(alpha))
v1sw_h12 = np.zeros(len(alpha))
v1sw_h14 = np.zeros(len(alpha))
v1sw_h29 = np.zeros(len(alpha))
i=0


alturas = [ 7, 9, 11, 13, 15, 30 ]

sol_gen2  = np.zeros(len(alpha))
sol_gen4  = np.zeros(len(alpha))    
sol_gen6  = np.zeros(len(alpha))    
sol_gen8  = np.zeros(len(alpha))    
sol_gen29  = np.zeros(len(alpha))    
sol_gen40  = np.zeros(len(alpha))    




ts_h2= np.zeros(len(alpha))
ts_h4= np.zeros(len(alpha))
ts_h6 = np.zeros(len(alpha))
ts_h8 = np.zeros(len(alpha))
ts_h29 = np.zeros(len(alpha))

for x in alpha:
    

    sol[i] = fsolve(f.Panel1(x).transcendental_eq ,[4.0] ) 
    
    v1sw_h4[i] = f.Panel1().v1s_w(x,sol[i],rho,4)
    v1sw_h8[i] = f.Panel1().v1s_w(x,sol[i],rho,8)
    v1sw_h10[i] = f.Panel1().v1s_w(x,sol[i],rho,10)
    v1sw_h12[i] = f.Panel1().v1s_w(x,sol[i],rho,12)
    v1sw_h14[i] = f.Panel1().v1s_w(x,sol[i],rho,14)
    v1sw_h29[i] = f.Panel1().v1s_w(x,sol[i],rho,29)
    # j = 0 
    # for y in alturas:
    sol_gen2[i] = fsolve(f.Panel4(1,x,3).trans4,[0.002] )         
    sol_gen4[i] = fsolve(f.Panel4(1,x,5).trans4,[0.002] )         
    sol_gen6[i] = fsolve(f.Panel4(1,x,7).trans4,[0.002] ) 
    sol_gen8[i] = fsolve(f.Panel4(1,x,9).trans4,[0.002] ) 
    sol_gen29[i] = fsolve(f.Panel4(1,x,30).trans4,[0.001] )
    sol_gen40[i] = fsolve(f.Panel4(1,x,40).trans4,[0.001] )
    
    ts_h2[i] = f.Stat_Temp(1,x,3).ts_w_uds(sol_gen2[i],rho)
    ts_h4[i] = f.Stat_Temp(1,x,5).ts_w_uds(sol_gen4[i],rho)
    ts_h6[i] = f.Stat_Temp(1,x,7).ts_w_uds(sol_gen6[i],rho)
    ts_h8[i] = f.Stat_Temp(1,x,9).ts_w_uds(sol_gen8[i],rho)
    ts_h29[i] = f.Stat_Temp(1,x,30).ts_w_uds(sol_gen29[i],rho)
  
    #    j +=1
    # print(alpha[i]) 
    i+=1
        
    # print(result.x)
# print(sol)

# ! d =   Dimension espacial
d = 3 


print("g tilde function")
print (f.Panel3(1.0,w,1.0,9).hp(d))

# !! FIGURA T_xy <--> alpha

fig22 = plt.figure(figsize=(10,8))

# !  h=2sigma
# plt.errorbar(H2['alpha'], H2['Txymean']/2, yerr=H2['Txy_error']/2,mfc="none",capsize=10,ms=12, color='g',marker="o",linestyle="",label= r" $H =  %2i  \sigma $ "   % 2) 

# # plt.plot(alpha[0:len(alpha)-20],ts_h2[0:len(alpha)-20],color= "g",label= r" $H = 2 \sigma $ ")

# plt.plot(alpha,f.Panel3(1.0,w,1.0,3).ts_w_uds_d(alpha,rho,d),linestyle = "dashdot",linewidth=1.5,color="g")

# plt.plot(alpha,f.Panel3(1.0,w,1.0,3).ts_w_uds_2sigma(alpha,rho),linestyle = ":",linewidth=1.5,color="g")

# ! h =4 sigma

# plt.plot(alpha[0:len(alpha)-1],v1sw_h4[0:len(alpha)-1]**2/2,linestyle = "dashdot",color= "b",label= r" $H \gg \sigma $ "  )


plt.errorbar(H4['alpha'], 2* H4['Txymean'], yerr=H4['Txy_error'],mfc="none",capsize=10,ms=12, color='C6',marker="o",linestyle="",label= r" $H =  %2i  \sigma $ "   % 4) 

# plt.plot(alpha,(f.Panel2(1.0,w,1.0,4).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C6")


plt.plot(alpha,f.Panel3(1.0,w,1.0,5).ts_w_uds(alpha,rho),linestyle = ":",linewidth=1.5,color="C6")

plt.plot(alpha[0:len(alpha)-1],ts_h4[0:len(alpha)-1],color= "C6")

# ! h =6 sigma

# plt.plot(alpha[0:len(alpha)-1],v1sw_h8[0:len(alpha)-1]**2/2,color= "C0")


plt.errorbar(H6['alpha'], 2*H6['Txymean'], yerr=H6['Txy_error'],mfc="none",capsize=10,ms=12, color='C0',marker="o",linestyle="",label= r" $H =  %2i  \sigma $ "   % 6) 

# # plt.plot(alpha,(f.Panel2(1.0,w,1.0,6).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C0")


plt.plot(alpha,f.Panel3(1.0,w,1.0,7).ts_w_uds(alpha,rho),linestyle = ":",linewidth=1.5,color="C0")

plt.plot(alpha[0:len(alpha)-1],ts_h6[0:len(alpha)-1],color= "C0")

# ! h =8 sigma

# plt.plot(alpha[0:len(alpha)-1],v1sw_h8[0:len(alpha)-1]**2/2,color= "C0")

# plt.errorbar(H8['alpha'], 2*H8['Txymean'], yerr=H8['Txy_error'],mfc="none",capsize=10,ms=12, color='C5',marker="o",linestyle="",label= r" $H =  %2i  \sigma $ "   % 8) 

# # plt.plot(alpha,(f.Panel2(1.0,w,1.0,8).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C5")


# plt.plot(alpha,f.Panel3(1.0,w,1.0,9).ts_w_uds(alpha,rho),linestyle = ":",linewidth=1.5,color="C5")

# plt.plot(alpha[0:len(alpha)-1],ts_h8[0:len(alpha)-1],color= "C5")


# plt.plot(alpha,f.Panel3(1.0,w,1.0,9).ts_w_uds_d(alpha,rho,d),linestyle = "dashdot",linewidth=1.5,color="C0")

# # ! h =10 sigma

# # plt.plot(alpha[0:len(alpha)-1],v1sw_h10[0:len(alpha)-1]**2/2,color= "C1")

# plt.errorbar(H10['alpha'], 2*H10['Txymean'], yerr=H10['Txy_error'],mfc="none",capsize=10,ms=12, color='C1',marker="s",linestyle="",label= r" $H =  %2i  \sigma $ "   % 10) 

# plt.plot(alpha,(f.Panel2(1.0,w,1.0,10).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C1")

# plt.plot(alpha,f.Panel3(1.0,w,1.0,11).ts_w_uds(alpha,rho),linestyle = ":",linewidth=1.5,color="C1")

# # ! h =12 sigma

# # plt.plot(alpha[0:len(alpha)-1],v1sw_h12[0:len(alpha)-1]**2/2,color= "C2")

# plt.errorbar(H12['alpha'], 2*H12['Txymean'], yerr=H12['Txy_error'],mfc="none",capsize=10,ms=12, color='C2',marker="d",linestyle="",label= r" $H =  %2i  \sigma $ "   % 12) 

# plt.plot(alpha,(f.Panel2(1.0,w,1.0,12).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C2")

# plt.plot(alpha,f.Panel3(1.0,w,1.0,13).ts_w_uds(alpha,rho),linestyle = ":",linewidth=1.5,color="C2")
# # ! h =14 sigma

# # plt.plot(alpha[0:len(alpha)-1],v1sw_h14[0:len(alpha)-1]**2/2,color= "C3")

# plt.errorbar(H14['alpha'], 2*H14['Txymean'], yerr=H14['Txy_error'],mfc="none",capsize=10,ms=12, color='C3',marker="*",linestyle="",label= r" $H =  %2i  \sigma $ "   % 14) 


# plt.plot(alpha,(f.Panel2(1.0,w,1.0,14).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C3")

# plt.plot(alpha,f.Panel3(1.0,w,1.0,15).ts_w_uds(alpha,rho),linestyle = ":",linewidth=1.5,color="C3")

# ! h =29 sigma

# plt.plot(alpha[0:len(alpha)-1],v1sw_h29[0:len(alpha)-1]**2/2,color= "C4")

# plt.errorbar(H29['alpha'], 2*H29['Txymean'], yerr=H29['Txy_error'],mfc="none",capsize=10,ms=12, color='C4',marker="^",linestyle="",label= r" $H =  %2i  \sigma $ "   % 29) 

# # plt.plot(alpha,(f.Panel2(1.0,w,1.0,29).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C4")

# plt.plot(alpha,f.Panel3(1.0,w,1.0,30).ts_w_uds(alpha,rho),linestyle = ":",linewidth=1.5,color="C4")

# plt.plot(alpha[0:len(alpha)-1],ts_h29[0:len(alpha)-1],color= "C4")

plt.ylabel( r'  $\frac{T^s}{m v_p^2}$  ',rotation=0.0,labelpad=30,fontsize=40)

plt.xlabel( r' $\alpha$ ',fontsize=40)


# plt.title( r' \textbf {Relación $T_{xy} \leftrightarrow \alpha$}' ,fontsize=40)

# plt.xlim(0.55,0.96)
plt.ylim(0,3500)
# plt.ylim(0,100000)

plt.tick_params(axis='x', labelsize=25)
plt.tick_params(axis='y', labelsize=25)

# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.legend(loc=0,fontsize=22)

plt.tight_layout()

if (save_temp_figures == True) :

    plt.savefig('txy_vs_alpha.pdf', dpi = 1200)
# plt.savefig('txy_vs_alpha.eps', dpi = 500, format="eps")

# !! FIGURA T_z <--> alpha

# v2sw_h8 = v1sw_h8**2/2*(1+f.Panel1.bets_2(alpha))
# v2sw_h10 = v1sw_h10**2/2*(1+f.Panel1.bets_2(alpha))
# v2sw_h12 = v1sw_h12**2/2*(1+f.Panel1.bets_2(alpha))
# v2sw_h14 = v1sw_h14**2/2*(1+f.Panel1.bets_2(alpha))
# v2sw_h29 = v1sw_h29**2/2*(1+f.Panel1.bets_2(alpha))

# fig23 = plt.figure(figsize=(10,8))

# # ! h = 8 sigma


# plt.plot(alpha[0:len(alpha)-1],v2sw_h8[0:len(alpha)-1],color= "C0")

# plt.errorbar(H8['alpha'], 2*H8['Tzmean'], yerr=H8['Tz_error'],mfc="none",capsize=10,ms=12, color='C0',marker="o",linestyle="",label= r" $H =  %2i  \sigma $ "   % 8) 

# plt.plot(alpha,(1+f.Panel1.bets_1(alpha))*(f.Panel2(1.0,w,1.0,8).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C0")

# # ! h = 10 sigma



# plt.plot(alpha[0:len(alpha)-1],v2sw_h10[0:len(alpha)-1],color= "C1")

# plt.errorbar(H10['alpha'], 2*H10['Tzmean'], yerr=H10['Tz_error'],mfc="none",capsize=10,ms=12, color='C1',marker="s",linestyle="",label= r" $H =  %2i  \sigma $ "   % 10) 

# plt.plot(alpha,(1+f.Panel1.bets_1(alpha))*(f.Panel2(1.0,w,1.0,10).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C1")

# # ! h =12 sigma



# plt.plot(alpha[0:len(alpha)-1],v2sw_h12[0:len(alpha)-1],color= "C2")

# plt.errorbar(H12['alpha'], 2*H12['Tzmean'], yerr=H12['Tz_error'],mfc="none",capsize=10,ms=12, color='C2',marker="d",linestyle="",label= r" $H =  %2i  \sigma $ "   % 12) 

# plt.plot(alpha,(1+f.Panel1.bets_1(alpha))*(f.Panel2(1.0,w,1.0,12).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C2")

# # ! h =14 sigma


# plt.plot(alpha[0:len(alpha)-1],v2sw_h14[0:len(alpha)-1],color= "C3")

# plt.errorbar(H14['alpha'], 2*H14['Tzmean'], yerr=H14['Tz_error'],mfc="none",capsize=10,ms=12, color='C3',marker="*",linestyle="",label= r" $H =  %2i  \sigma $ "   % 14) 

# plt.plot(alpha,(1+f.Panel1.bets_1(alpha))*(f.Panel2(1.0,w,1.0,14).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C3")

# # ! h =29 sigma

# plt.plot(alpha[0:len(alpha)-1],v2sw_h29[0:len(alpha)-1],color= "C4")

# plt.errorbar(H29['alpha'], 2*H29['Tzmean'], yerr=H29['Tz_error'],mfc="none",capsize=10,ms=12, color='C4',marker="^",linestyle="",label= r" $H =  %2i  \sigma $ "   % 29) 

# plt.plot(alpha,(1+f.Panel1.bets_1(alpha))*(f.Panel2(1.0,w,1.0,29).vsw_mano(rho,alpha))**2/2,linestyle = "dashdot",linewidth=1.5,color="C4")

# plt.ylabel( r'  $\frac{T^s_{z}}{m v_p^2}$  ',rotation=0.0,labelpad=30,fontsize=40)

# plt.xlabel( r' $\alpha$ ',fontsize=40)

# plt.tick_params(axis='x', labelsize=25)
# plt.tick_params(axis='y', labelsize=25)

# # plt.title( r' \textbf {Relación $T_{z} \leftrightarrow \alpha$}' ,fontsize=40)

# plt.xlim(0.55,0.99)
# plt.ylim(0,3000)

# # plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
# plt.legend(loc=0,fontsize=22)


# plt.tight_layout()

# if (save_temp_figures == True) :
    
#     plt.savefig('tz_vs_alpha.pdf', dpi = 1200)
# # plt.savefig('tz_vs_alpha.eps', dpi = 500, format="eps")




fig, ax1=plt.subplots(figsize=(14,8))


# ax1.plot(alpha,f.Panel1.bets_1(alpha) ,linewidth=1.5,color="k")
# ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 30).betas_h2(alpha) ,linewidth=1.5,color="k")
# ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 2).betas_d(alpha,d) ,linewidth=3.5,alpha = 0.5 ,linestyle= "dashdot" ,color="k")


# print(sol_gen)


# ax1.plot(ba['alpha'],ba['beta'] ,linewidth=1.5,linestyle= "dashdot",color='b' )
ax1.plot(alpha,sol_gen4,linewidth=1.5,linestyle= "dashdot",color='C6' )
ax1.plot(alpha,sol_gen6,linewidth=1.5,linestyle= "dashdot",color='C2' )
ax1.plot(alpha,sol_gen8,linewidth=1.5,linestyle= "dashdot",color='C0' )
ax1.plot(alpha,sol_gen29,linewidth=1.5,linestyle= "dashdot",color='C4' )
ax1.plot(alpha,sol_gen40,linewidth=1.5,linestyle= "dashdot",color='g' )
# ax1.plot(alpha,sol ,linewidth=1.5,linestyle= "dashdot",color='b',label=r" $H \gg \sigma$ "  )
ax1.plot(alpha,sol_gen2,linewidth=1.5,linestyle= "dashdot",color='b',label=r" $H = 2 \sigma$ " )


# ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 6).betas(alpha),linestyle = ":",color='C2')
ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 5).betas_d(alpha,d),linestyle = ":",color='C6')
# ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 7).betas_d(alpha,d),linestyle = ":",color='C2')
# ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 9).betas(alpha),linestyle = ":",color='C0')
# ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 30).betas(alpha),linestyle = ":",color='C4')

ax1.errorbar(ba_sim4['alpha'], ba_sim4['beta'], yerr=ba_sim4['errorbeta'],mfc="none",capsize=10,ms=12, color='C6',marker="o",linestyle="",label=r" $H= 4\sigma$ ") 
ax1.errorbar(ba_sim2['alpha'], ba_sim2['beta']/2, yerr=ba_sim2['errorbeta'],mfc="none",capsize=10,ms=12, color='b',marker="o",linestyle="",label=r" $H= 2\sigma$ ") 


# ax1.errorbar(ba_sim6['alpha'], ba_sim6['beta'], yerr=ba_sim6['errorbeta'],mfc="none",capsize=10,ms=12, color='C2',marker="o",linestyle="",label=r" $H= 6\sigma$ ") 

# ax1.errorbar(ba_sim8['alpha'], ba_sim8['beta'], yerr=ba_sim8['errorbeta'],mfc="none",capsize=10,ms=12, color='C0',marker="o",linestyle="",label=r" $H= 8\sigma$ ") 
# ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 9).betas_d(alpha,d),linestyle = "dashdot",color='C0')

# ax1.errorbar(ba_sim10['alpha'], ba_sim10['beta'], yerr=ba_sim10['errorbeta'],mfc="none",capsize=10,ms=12, color='C1',marker="s",linestyle="",label=r" $H= 10\sigma$ ") 
# # # ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 10).betas(alpha),linestyle = ":",color='C1')
# ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 11).betas_d(alpha,d),linestyle = "dashdot",color='C1')

# ax1.errorbar(ba_sim12['alpha'], ba_sim12['beta'], yerr=ba_sim12['errorbeta'],mfc="none",capsize=10,ms=12, color='C2',marker="d",linestyle="",label=r" $H= 12\sigma$ ") 
# # ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 12).betas(alpha),linestyle = ":",color='C2')

# ax1.errorbar(ba_sim14['alpha'], ba_sim14['beta'], yerr=ba_sim14['errorbeta'],mfc="none",capsize=10,ms=12, color='C3',marker="*",linestyle="",label=r"$H= 14\sigma$ ") 
# # ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 14).betas(alpha),linestyle = ":",color='C3')


ax1.errorbar(ba_sim29['alpha'], ba_sim29['beta'], yerr=ba_sim29['errorbeta'],mfc="none",capsize=10,ms=12, color='C4',marker="^",linestyle="",label=r" $H= 29\sigma$ ") 

ax1.errorbar(ba_sim29p['alpha'], ba_sim29p['beta'], yerr=ba_sim29p['errorbeta'],mfc="none",capsize=10,ms=12, color='C8',marker="^",linestyle="",label=r" $H= 29\sigma $ p ") 


# ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 29).betas(alpha),linestyle = ":",color='C4')
# ax1.plot(alpha,f.Panel3(1.0,0.001,1.0, 30).betas_d(alpha,d),linestyle = "dashdot",color='C4')





# ax1.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

ax1.set_ylabel( r'  $\beta_s$  ',rotation=0.0,labelpad=30,fontsize=40)

ax1.set_xlabel( r' $\alpha$ ', fontsize=40)

ax1.tick_params(axis='x', labelsize=25)
ax1.tick_params(axis='y', labelsize=25)
ax1.legend(loc="lower left",fontsize=19)


# ax1.set_xlim(0.55,1.0)
ax1.set_ylim(0.0,2.0)
# plt.xlim(0.60,1.0)

# ax1.set_title( r' \textbf {Relación $\beta \leftrightarrow \alpha$}' ,fontsize=40)


# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
ax2 = plt.axes([0,0,1,1])


# Manually set the position and relative size of the inset axes within ax1


ip = InsetPosition(ax1, [0.45,0.45,0.5,0.5])


ax2.set_axes_locator(ip)

ax2.plot(alpha,f.Panel1.bets_1(alpha) ,linewidth=1.5,color="k" )


# ax2.plot(alpha,sol ,linewidth=1.5,linestyle= "dashdot",color='b')


# ax1.plot(ba['alpha'],ba['beta'],linestyle = ":",color= "C0")

ax2.plot(alpha,sol_gen4,linewidth=1.5,linestyle= "dashdot",color='C6' )
# ax2.plot(alpha,sol_gen6,linewidth=1.5,linestyle= "dashdot",color='C2' )
# ax2.plot(alpha,sol_gen8,linewidth=1.5,linestyle= "dashdot",color='C0' )
# ax2.plot(alpha,sol_gen29,linewidth=1.5,linestyle= "dashdot",color='C4' )
# ax2.plot(alpha,sol ,linewidth=1.5,linestyle= "dashdot",color='b' )


ax2.errorbar(ba_sim4['alpha'], ba_sim4['beta'], yerr=ba_sim4['errorbeta'],mfc="none",capsize=10,ms=12, color='C6',marker="o",linestyle="") 

ax2.plot(alpha,f.Panel3(1.0,0.001,1.0, 5).betas_d(alpha,d),linestyle = ":",color='C6')


# ax2.errorbar(ba_sim6['alpha'], ba_sim6['beta'], yerr=ba_sim6['errorbeta'],mfc="none",capsize=10,ms=12, color='C2',marker="o",linestyle="") 
# # ax2.plot(alpha,f.Panel3(1.0,0.001,1.0, 6).betas(alpha),linestyle = ":",color='C2')
# ax2.plot(alpha,f.Panel3(1.0,0.001,1.0, 7).betas_d(alpha,d),linestyle = ":",color='C2')


# ax2.errorbar(ba_sim8['alpha'], ba_sim8['beta'], yerr=ba_sim8['errorbeta'],mfc="none",capsize=10,ms=12, color='C0',marker="o",linestyle="") 
# # ax2.plot(alpha,f.Panel3(1.0,0.001,1.0, 8).betas(alpha),linestyle = ":",color='C0')
# ax2.plot(alpha,f.Panel3(1.0,0.001,1.0, 9).betas_d(alpha,d),linestyle = ":",color='C0')
# ax2.errorbar(ba_sim10['alpha'], ba_sim10['beta'], yerr=ba_sim10['errorbeta'],mfc="none",capsize=10,ms=12, color='C1',marker="s",linestyle="") 
# # ax2.plot(alpha,f.Panel3(1.0,0.001,1.0, 10).betas(alpha),linestyle = ":",color='C1')

# ax2.errorbar(ba_sim12['alpha'], ba_sim12['beta'], yerr=ba_sim12['errorbeta'],mfc="none",capsize=10,ms=12, color='C2',marker="d",linestyle="") 
# # ax2.plot(alpha,f.Panel3(1.0,0.001,1.0, 12).betas(alpha),linestyle = ":",color='C2')

# ax2.errorbar(ba_sim14['alpha'], ba_sim14['beta'], yerr=ba_sim14['errorbeta'],mfc="none",capsize=10,ms=12, color='C3',marker="*",linestyle="") 
# # ax2.plot(alpha,f.Panel3(1.0,0.001,1.0, 14).betas(alpha),linestyle = ":",color='C3')

# ax2.errorbar(ba_sim29['alpha'], ba_sim29['beta'], yerr=ba_sim29['errorbeta'],mfc="none",capsize=10,ms=12, color='C4',marker="^",linestyle="") 
# # # ax2.plot(alpha,f.Panel3(1.0,0.001,1.0, 29).betas(alpha),linestyle = ":",color='C4')
# ax2.plot(alpha,f.Panel3(1.0,0.001,1.0, 30).betas_d(alpha,d),linestyle = ":",color='C4')
# ax2.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
ax2.legend(loc=0)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
# Some ad hoc tweaks.

ax2.set_xlim(0.84,1.0)
ax2.set_ylim(0.0,0.50)


plt.tight_layout()

if (save_beta_figure == True) :
#  plt.savefig('beta_vs_alpha.pdf', dpi = 1200)
    plt.savefig('beta_vs_alpha_teo.pdf',dpi=1200)


plt.show()