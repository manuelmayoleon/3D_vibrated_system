# remove
from numpy.core.function_base import linspace
import pandas as pd
import numpy as np
import scipy . stats as ss
import math
import matplotlib . mlab as mlab
import matplotlib . pyplot as plt
from scipy import optimize
import matplotlib


import os
import glob
  
  
# use glob to get all the csv files 
# in the folder
path = os.getcwd()
path = path+"/Txymean_vs_alpha/"
csv_files = glob.glob(os.path.join(path, "*.dat"))
csv_files = sorted(csv_files, key=len)
print(len(csv_files))

names = list()
elements=list()
# loop over the list of csv files
for f in csv_files:
        
        # read the csv file
        # df= pd.read_csv(f)
       
        df= pd.read_csv(f,header=0,sep='\t\t+',engine='python')
        elements.append(df)
        # print the location, filename and name
        print('Location:', f)
        file_name = f.split("/")[-1]
        print('File Name:', file_name )
        name = file_name.split(".")[0]
        print('Name:', name)
        names.append(name)
          
        # print the content
        print('Content:')
        print(df)  
        
        
print(names)
print(elements[1])

H8 =  elements[0]  
H10 = elements[-1]   
H12 = elements[2]
H14 = elements[1]
H29 = elements[3]

# H8.columns = H8.iloc[0]

H8.columns = H8.columns.str.strip()

print( H8.columns )








