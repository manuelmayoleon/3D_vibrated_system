 # -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 21:10:41 2018

@author: Juan Petit
"""

import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit
import os


# alpha = [10, 20, 30, 40, 50, 60, 65, 70, 75, 80, 85, 90, 95]
# alpha = [60, 65, 70, 75, 80, 85, 90, 95]
alpha = [ 65, 70, 75]
L = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

mass = 1

for q in alpha:

    for m in L:
              
        if not os.path.exists('alpha_p%d/pressures' %q): os.makedirs('alpha_p%d/pressures' %q)
        
        count = 0
        # Iterate directory
        for path in os.listdir('alpha_p%d/R%d' %(q,m)):
            # check if current path is a file
            if os.path.isfile(os.path.join('alpha_p%d/R%d' %(q,m), path)):
                count += 1
                
        print('File count:', count)
        
        ini = 0
        fin = count - 1
        
        NU_rat = []
        time = []
        pressure_contact = []
        pressure_kin = []
        Gtemp = []
        Tx = []
        Ty = []
        Tz = []
    
    
        for k in range(ini, fin):
      
            data = open('alpha_p%d/R%d/3Dgas_%d.dat' %(q,m,k))
    
            firstline = data.readlines()[1:2]
    
            for p in firstline:
                FL = p.split()
                if FL != []:
                    time.append(float(FL[0]))
                    N = float(FL[1])
                    Lx = float(FL[2])
                    Ly = float(FL[3])
                    Lz = float(FL[4])
                    R = float(FL[5])

    
            data = open('alpha_p%d/R%d/3Dgas_%d.dat' %(q,m,k))
    
            lines = data.readlines()[3:]
    
            x = []
            y = []
            z = []
            vx = []
            vy = []
            vz = []
    
            for line in lines:
                pos = line.split()
                if pos != []:
                    x.append(float(pos[0]))
                    y.append(float(pos[1]))
                    z.append(float(pos[2]))
                    vx.append(float(pos[3]))
                    vy.append(float(pos[4]))
                    vz.append(float(pos[5]))
                
            ux = np.mean(vx)
            uy = np.mean(vy)
            uz = np.mean(vz)
            pkin = 0.0 
            temp = 0.0
            tempx = 0.0
            tempy = 0.0
            tempz = 0.0
        
            
            for i in range(0, len(x)):
                pkin = pkin + mass*((vx[i] - ux)**2 + (vy[i] - uy)**2 + (vz[i] - uz)**2)  
                temp = temp + 0.5*mass*( (vx[i] - ux)**2 + (vy[i] - uy)**2 + (vz[i] - uz)**2)  
                tempx = tempx + 0.5*mass*( (vx[i] - ux)**2 )  
                tempy = tempy + 0.5*mass*( (vy[i] - uy)**2 )  
                tempz = tempz + 0.5*mass*( (vz[i] - uz)**2 )  
    
                
            pressure_kin.append(pkin/(3*Lx*Ly*Lz))                  
            Gtemp.append(temp/len(x))   
            Tx.append(tempx/len(x))
            Ty.append(tempy/len(x))
            Tz.append(tempz/len(x))
                
            Vp = 0.0
            for i in range(0,len(x)):
                Vp = Vp + (4.0/3.0)*np.pi*R**3  # Three DIMENSIONS
            NU_rat.append(Vp/(Lx*Ly*Lz))
    
            
               
            print('Analysed DATA: alpha Rep, state  ---> %d %d %d' %(q,m,k))
                
            data = np.column_stack((time, NU_rat, pressure_kin, Gtemp, Tx, Ty, Tz))
            np.savetxt(os.path.join('alpha_p%d/pressures' %q, 'Data_R%d.dat' %(m)), data, fmt='%.10f', delimiter='\t\t ', header='time \t\t NU \t\t Pkin \t\t\t Gtemp \t\t Tx \t\t Ty \t\t Tz')
