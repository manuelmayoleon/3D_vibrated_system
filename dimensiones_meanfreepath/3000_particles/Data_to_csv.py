#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 17:31:12 2022

@author: juanpetit
"""

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import os
import pandas as pd
import os
import glob
  
path1 = 'data_csv'
path = os.getcwd()
path = path + '/R1/'

dat_files = glob.glob(os.path.join(path,"3Dgas*"+"*.dat"))
dat_files = sorted(dat_files, key=len)
print(len(dat_files))
ini = 0
fin = len(dat_files)

for i in range(ini, fin):
    
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []
    
    if not os.path.exists(path1): os.makedirs(path1)
    
    # print(path)
    data = open(path+'3Dgas_%d.dat' %i)
    # data = dat_files[i]
    lines = data.readlines()[1:2]
    for line in lines:
        pos = line.split()
        if pos != []:
            t = float(pos[0])
            N = float(pos[1])
            Lx = float(pos[2])
            Ly = float(pos[3])
            Lz = float(pos[4])
            R = float(pos[5])
                       
    
    data = open(path+'3Dgas_%d.dat' %i)
    lines = data.readlines()[3:]
    for line in lines:
        pos = line.split()
        if pos != []:
            x.append(float(pos[0]))
            y.append(float(pos[1]))
            z.append(float(pos[2]))
            vx.append(float(pos[3]))
            vy.append(float(pos[4]))
            vz.append(float(pos[5]))
            
    Rad = R*np.ones(len(x))
    
    X1 = []
    for k in range(len(x)):
        X1.append([x[k], y[k], z[k], Rad[k]])
        
    arr = np.asarray(X1)
    pd.DataFrame(arr).to_csv(os.path.join(path1, 'state.csv.%d' %i), index_label = "Index", header  = ['x','y','z', 'R'])

    print('save fig ---> %d' %i)

