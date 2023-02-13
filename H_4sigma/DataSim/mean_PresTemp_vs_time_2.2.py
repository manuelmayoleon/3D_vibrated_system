 # -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 21:10:41 2018

@author: Juan Petit
"""

import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit
import os

path1 = 'average_PresTemp'

# ! Acordarse de cambiarlo 
R = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,14,15,16,17,18,19,20]
# R = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# alpha = [60, 65, 70, 75, 80, 85, 90, 95]

alpha = [60, 65, 70, 75, 80]

if not os.path.exists(path1): os.makedirs(path1)

#-----------------------------------------------------------
# CHEKING WHAT DATA FILE HAS THE MINIMUM LENGTH

# for q in alpha:
    
#     D = []

#     for m in R:
        
#         Tx = []
    
#         data = open('alpha_p%d/pressures/Data_R%d.dat' %(q,m))
            
#         lines = data.readlines()[1:]
        
#         for line in lines:
#             pos = line.split()
#             if pos != []:
#                 Tx.append(float(pos[4]))
        
#         D.append(len(Tx))
    
#     minD = min(D)

#-----------------------------------------------------------

error = []
beta_mean = []

for q in alpha:
    
    #-----------------------------------------------------------
    # CHEKING WHAT DATA FILE HAS THE MINIMUM LENGTH
    
    D = []

    for m in R:
        
        Tx = []
    
        data = open('alpha_p%d/pressures/Data_R%d.dat' %(q,m))
            
        lines = data.readlines()[1:]
        
        for line in lines:
            pos = line.split()
            if pos != []:
                Tx.append(float(pos[4]))
        
        D.append(len(Tx))
    
    minD = min(D)
    
    #-----------------------------------------------------------

    
    Pcon_matrix = []
    Pkin_matrix = []
    Gtemp_matrix = []
    Tx_matrix = []
    Ty_matrix = []
    Tz_matrix = []
    beta = []
    
    for m in R:
    
        data = open('alpha_p%d/pressures/Data_R%d.dat' %(q,m))
            
        lines = data.readlines()[1:]
        
        t = []
        NU = []
        Pkin = []
        Gtemp = []
        Tx = []
        Ty = []
        Tz = []
        
        for line in lines:
            pos = line.split()
            if pos != []:
                t.append(float(pos[0]))
                NU.append(float(pos[1]))
                Pkin.append(float(pos[2]))
                Gtemp.append(float(pos[3]))
                Tx.append(float(pos[4]))
                Ty.append(float(pos[5]))
                Tz.append(float(pos[6]))
                
        Pkin_matrix.append(Pkin)
        Gtemp_matrix.append(Gtemp)
        Tx_matrix.append(Tx)
        Ty_matrix.append(Ty)
        Tz_matrix.append(Tz)
        
        beta.append(((2*Tz[-1])/((Tx[-1] + Ty[-1]))) - 1)
    
    error_sum = 0.0
    beta_mean.append(np.mean(beta))
    
    # ERROR CALCULATION 
    
    for k in range(len(beta)):
        error_sum = error_sum + (beta[k] - np.mean(beta))**2
    
    error.append(np.sqrt(error_sum/(len(beta) - 1)))
    
    Pkin_ave = []
    Gtemp_ave = []
    Tx_ave = []
    Ty_ave = []
    Tz_ave = []
    
    n = len(Tz_matrix)
    
    for j in range(0, minD):
        
        pkinsum = 0.0
        gtempsum = 0.0
        txsum = 0.0
        tysum = 0.0
        tzsum = 0.0
        
        for k in range(n):
            
            pkinsum = pkinsum + Pkin_matrix[k][j]
            gtempsum = gtempsum + Gtemp_matrix[k][j]
            txsum = txsum + Tx_matrix[k][j] 
            tysum = tysum + Ty_matrix[k][j] 
            tzsum = tzsum + Tz_matrix[k][j] 
            
        Pkin_ave.append(pkinsum/n)
        Gtemp_ave.append(gtempsum/n)
        Tx_ave.append(txsum/n)
        Ty_ave.append(tysum/n)
        Tz_ave.append(tzsum/n)
    
    data = np.column_stack((t[0:minD], NU[0:minD], Pkin_ave, Gtemp_ave, Tx_ave, Ty_ave, Tz_ave))
    np.savetxt(os.path.join(path1, 'ave_alpha_p%d.dat' %(q)), data, fmt='%.10f', delimiter='\t\t ', header='time \t\t NU \t\t Pkin \t\t\t Gtemp \t\t Tx \t\t Ty \t\t Tz')

alphas = []
for p in range(len(alpha)):
    alphas.append(alpha[p]/100.0)

data = np.column_stack((alphas, beta_mean,  error))
np.savetxt(os.path.join(path1, 'beta_vs_alphas.dat'), data, fmt='%.10f', delimiter='\t\t ', header='alpha \t\t beta \t\t error_beta')


# data = np.column_stack((alpha/100.0, beta_mean,  error))
# np.savetxt(os.path.join(path1, 'beta_vs_alpha_p%d.dat' %(alpha)), data, fmt='%.10f', delimiter='\t\t ', header='alpha \t\t beta \t\t\t error')

    
    