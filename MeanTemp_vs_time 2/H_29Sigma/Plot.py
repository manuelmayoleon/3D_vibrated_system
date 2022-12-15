 # -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 21:10:41 2018

@author: Juan Petit
"""

import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit
import os

alpha = [10, 20, 30, 40, 50, 60, 65, 70, 75, 80, 85, 90, 95]

# if not os.path.exists(path1): os.makedirs(path1)
# if not os.path.exists('RattlersFree/SRp%d' %m): os.makedirs('RattlersFree/SRp%d' %m)

beta = []

for m in alpha:

    data = open('ave_alpha_p%d.dat' %(m))
        
    lines = data.readlines()[1:]
    
    t = []
    Gtemp = []
    Tx = []
    Ty = []
    Tz = []
    
    for line in lines:
        pos = line.split()
        if pos != []:
            t.append(float(pos[0]))
            Gtemp.append(float(pos[3]))
            Tx.append(float(pos[4]))
            Ty.append(float(pos[5]))
            Tz.append(float(pos[6]))
    
    plt.figure(1, figsize = (12.5,12))
  
               
    plt.subplot(3, 1, 1)
    
    if m == 50:
        plt.plot(t, Tx, '-', mfc='none', color = 'purple', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))  
    if m == 60:
        plt.plot(t, Tx, '-', mfc='none', color = 'cyan', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))  
    if m == 70:
        plt.plot(t, Tx, '-', mfc='none', color = 'magenta', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))    
    if m == 80:
        plt.plot(t, Tx, '-', mfc='none', color = 'red', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))
    if m == 85:
        plt.plot(t, Tx, '-', mfc='none', color = 'blue', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))
    if m == 90:
        plt.plot(t, Tx, '-', mfc='none', color = 'green', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))
    if m == 95:
        plt.plot(t, Tx, '-', mfc='none', color = 'k', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))
        
    plt.xlim(0.0, 5e5)
    plt.ylim(0, 0.001)
    locs,labels = plt.yticks()
    plt.yticks(locs, map(lambda x: "%.4f" % x, locs*1.0), fontsize = 20)
    locs,labels = plt.xticks()
    plt.xticks(locs, map(lambda x: "%.0f" % x, locs*1.0), fontsize = 20)
    # plt.ylabel(r"P", fontsize = 30)
    plt.legend(loc = 'upper right', fontsize = 14)
    # plt.title(r'$x_{\mathrm{S}} = 0.1$', fontsize  = 24)
    plt.text(x = 95, y = 2, s = r'$\phi = %.2f$' %(m/100), fontsize = 18)
    plt.ylabel(r'$T_{x}$', fontsize = 26)
    # plt.xlabel(r"$t$", fontsize = 30)
    # plt.tight_layout()
    # plt.show()
        
        # xvalues = np.asarray(t) 
        
        # def fitfun(X, a, b, c):
        #     return  a - b*(X**c)

        # init_vals = [0.5, 0.5, 1]    
        # best_vals, covar1 = curve_fit(fitfun, xvalues, Ptot, p0=init_vals, maxfev = 100000000)
        # fig4, = plt.plot(xvalues,fitfun(xvalues,best_vals[0], best_vals[1], best_vals[2]), 'k', lw = 2.5)
        
        # print('parameters are a, b, c', best_vals)
        
    plt.subplot(3, 1, 2)
    
    if m == 50:
        plt.plot(t, Ty, '-', mfc='none', color = 'purple', markersize = 12, label = r'$\alpha = %.2f$' %(m/100)) 
    if m == 60:
        plt.plot(t, Ty, '-', mfc='none', color = 'cyan', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))  
    if m == 70:
        plt.plot(t, Ty, '-', mfc='none', color = 'magenta', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))    
    if m == 80:
        plt.plot(t, Ty, '-', mfc='none', color = 'red', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))
    if m == 85:
        plt.plot(t, Ty, '-', mfc='none', color = 'blue', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))
    if m == 90:
        plt.plot(t, Ty, '-', mfc='none', color = 'green', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))
    if m == 95:
        plt.plot(t, Ty, '-', mfc='none', color = 'k', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))

        
    plt.xlim(0.0, 5e5)
    plt.ylim(0, 0.001)
    locs,labels = plt.yticks()
    plt.yticks(locs, map(lambda x: "%.4f" % x, locs*1.0), fontsize = 20)
    locs,labels = plt.xticks()
    plt.xticks(locs, map(lambda x: "%.0f" % x, locs*1.0), fontsize = 20)
    # plt.ylabel(r"P", fontsize = 30)
    plt.legend(loc = 'upper right', fontsize = 14)
    # plt.title(r'$x_{\mathrm{S}} = 0.1$', fontsize  = 24)
    plt.text(x = 95, y = 2, s = r'$\phi = %.2f$' %(m/100), fontsize = 18)
    plt.ylabel(r'$T_{y}$', fontsize = 26)
    # plt.xlabel(r"$t$", fontsize = 30)
        
    plt.subplot(3, 1, 3)
    
    if m == 50:
        plt.plot(t, Tz, '-', mfc='none', color = 'purple', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))    
    if m == 60:
        plt.plot(t, Tz, '-', mfc='none', color = 'cyan', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))    
    if m == 70:
        plt.plot(t, Tz, '-', mfc='none', color = 'magenta', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))    
    if m == 80:
        plt.plot(t, Tz, '-', mfc='none', color = 'red', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))
    if m == 85:
        plt.plot(t, Tz, '-', mfc='none', color = 'blue', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))
    if m == 90:
        plt.plot(t, Tz, '-', mfc='none', color = 'green', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))
    if m == 95:
        plt.plot(t, Tz, '-', mfc='none', color = 'k', markersize = 12, label = r'$\alpha = %.2f$' %(m/100))

        # plt.plot(t, Gtemp, '-', mfc='none', color = 'green', markersize = 12, label = r'$T_{g}$', lw = 2)
        # plt.plot(t, Tx, '-', mfc='none', color = 'red', markersize = 12, label = r'$T_{x}$', lw = 2)
        # plt.plot(t, Ty, '-', mfc='none', color = 'blue', markersize = 12, label = r'$T_{y}$', lw = 2)
        # plt.plot(t, Tz, '-', mfc='none', color = 'k', markersize = 12, label = r'$T_{z}$', lw = 2)

        
    plt.xlim(0.0, 5e5)
    plt.ylim(0, 0.001)
    locs,labels = plt.yticks()
    plt.yticks(locs, map(lambda x: "%.4f" % x, locs*1.0), fontsize = 20)
    locs,labels = plt.xticks()
    plt.xticks(locs, map(lambda x: "%.0f" % x, locs*1.0), fontsize = 20)
    # plt.ylabel(r"P", fontsize = 30)
    plt.legend(loc = 'upper right', fontsize = 14)
    # plt.title(r'$x_{\mathrm{S}} = 0.1$', fontsize  = 24)
    plt.text(x = 95, y = 2, s = r'$\phi = %.2f$' %(m/100), fontsize = 18)
    plt.ylabel(r'$T_{z}$', fontsize = 26)
    plt.xlabel(r"$t$", fontsize = 30)
    
    plt.savefig('T_vs_t_movewall.pdf', dpi = 500)
    # plt.close()
    
    
    
data = open('beta_vs_alphas.dat')
    
lines = data.readlines()[1:]

alpha = []
beta = []
errorbeta = []

for line in lines:
    pos = line.split()
    if pos != []:
        alpha.append(float(pos[0]))
        beta.append(float(pos[1]))
        errorbeta.append(float(pos[2]))
    
beta_theory = []
a = np.linspace(0, 1, 100)
for i in range(len(a)):
    beta_theory.append((10*(1 - a[i]))/(1 + 3*a[i]))

plt.figure(2, figsize = (8,6))

if  m == 80:
    plt.errorbar(alpha, beta, yerr = errorbeta, mfc='none', fmt='o', color = 'red', markersize = 12, label = r'$\frac{T^{s}_{2}}{T^{s}_{1}} -1$')        
else:
    plt.errorbar(alpha, beta, yerr = errorbeta, mfc='none', fmt='o', color = 'red', markersize = 12)
   
    plt.plot(a, beta_theory, '-', mfc='none', color = 'red')
    plt.xlim(0, 1)
    plt.ylim(0, 10)
    locs,labels = plt.yticks()
    plt.yticks(locs, map(lambda x: "%.1f" % x, locs*1.0), fontsize = 20)
    locs,labels = plt.xticks()
    plt.xticks(locs, map(lambda x: "%.2f" % x, locs*1.0), fontsize = 20)
    # plt.ylabel(r"P", fontsize = 30)
    plt.legend(loc = 'upper right', fontsize = 18)
    # plt.title(r'$x_{\mathrm{S}} = 0.1$', fontsize  = 24)
    # plt.ylabel(r'$\beta$', fontsize = 26)
    plt.xlabel(r"$\alpha$", fontsize = 30)
    plt.tight_layout()

plt.savefig('beta_vs_alpha.pdf', dpi = 500)



  