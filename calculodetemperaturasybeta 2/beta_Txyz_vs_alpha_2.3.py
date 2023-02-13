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


R = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
alpha = [60, 65, 70, 75, 80, 85, 90, 95]

omega = 0.001

if not os.path.exists(path1): os.makedirs(path1)

beta_error = []
beta_mean = []
Txy_mean = []
Tz_mean = []
Txy_error = []
Tz_error = []

for q in alpha:
    
    Txy_matrix = []
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
                
        Txy_matrix.append((Tx[-1] + Ty[-1])/(2*omega**2))  # "removing units" 
        Tz_matrix.append(Tz[-1]/omega**2)
        
        beta.append(((Tz[-1])/(0.5*(Tx[-1] + Ty[-1]))) - 1)
    
    beta_error_sum = 0.0
    Txy_error_sum = 0.0
    Tz_error_sum = 0.0
    
    beta_mean.append(np.mean(beta))
    Txy_mean.append(np.mean(Txy_matrix))
    Tz_mean.append(np.mean(Tz_matrix))
    
    # ERROR CALCULATION 
    
    for k in range(len(beta)):
        beta_error_sum = beta_error_sum + (beta[k] - np.mean(beta))**2
        Txy_error_sum = Txy_error_sum + (Txy_matrix[k] - np.mean(Txy_matrix))**2
        Tz_error_sum = Tz_error_sum + (Tz_matrix[k] - np.mean(Tz_matrix))**2

    beta_error.append(np.sqrt(beta_error_sum/(len(beta) - 1)))
    Txy_error.append(np.sqrt(Txy_error_sum/(len(Txy_matrix) - 1)))
    Tz_error.append(np.sqrt(Tz_error_sum/(len(Tz_matrix) - 1)))
    
   
alphas = []
for p in range(len(alpha)):
    alphas.append(alpha[p]/100.0)

data = np.column_stack((alphas, beta_mean,  beta_error, Txy_mean, Txy_error, Tz_mean, Tz_error))
np.savetxt(os.path.join(path1, 'beta_Txyz_vs_alphas.dat'), data, fmt='%.10f', delimiter='\t\t ', header='alpha \t\t betamean \t\t error_beta \t\t Txymean \t\t Txy_error \t\t Tzmean \t\t Tz_error')

    
    