from numpy.core.function_base import linspace
import pandas as pd
import numpy as np
import scipy . stats as ss
import math
import matplotlib . mlab as mlab
import matplotlib . pyplot as plt
import matplotlib.animation as manimation
from scipy import optimize
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
import os 
# PARA USAR A LA HORA DE GUARDAR DOS COLUMNAS EN UN ARCHIVO
import csv
import glob
from re import search


import cv2

# r=np.genfromtxt("pos.txt",names=["rx","ry"])
# v = np.genfromtxt("velocidad.txt", names=["vx","vy"])
# rinicial=np.genfromtxt("posiciones_init.txt",names=["rx","ry"])
# vinit=np.genfromtxt("velocidad_init.txt", names=["vx","vy"])
# temp=np.genfromtxt("temperaturas.txt" ,names= ["y","z"])
# tiempo=np.genfromtxt("tiemposdecol.txt",names=["t"])


# ! modificar el path si cambias el numero de particulas, la densidad y el alpha 

path = os.getcwd()


# ! Find all .dat files and put into array
csv_files = glob.glob(os.path.join(path, "*.dat"))

print(csv_files)

# temp1= pd.read_csv(csv_files[0],header=0,sep='\s+' )
# temp2= pd.read_csv(csv_files[1],header=0,sep='\s+' )


names = list()

# loop over the list of csv files
for ii in csv_files:
        
        # read the csv file
        # df= pd.read_csv(f)
         # print the location, filename and name
        print('Location:', ii)
        file_name = ii.split("/")[-1]
        # print('File Name:', file_name )
        name = file_name.split(".")[0]
        print('Name:', name)
        names.append(name)
        len(name)
        if search("T1_", name):
           temp1= pd.read_csv(ii,header=0,sep='\t\t+')
        if search("T2_", name):
          temp2= pd.read_csv(ii,header=0,sep='\t\t+' )

# ! Para quitar los espacios al principio y al final de la cabecera
temp1.columns = temp1.columns.str.strip()
temp1 = temp1.rename(columns={'# z': 'z'})
temp2.columns = temp2.columns.str.strip()
temp2 = temp2.rename(columns={'# z': 'z'})


# for filename in os.listdir(path):
#     if filename.startswith('pos_0.995_tiempo'):
#         ri[i] = pd.read_csv(filename,header=None,sep='\s+',names=["ryy","rzz"])
#     i=i+1   

# r= pd.read_csv("pos_0.995_tiempo_50000000.txt",header=None,sep='\s+',names=["rx","ry"])
# r= pd.read_csv("pos_0.995_tiempo_50000000.txt",header=None,sep='\s+',names=["rx","ry"])
# r= pd.read_csv(path+"pos_0.980.txt",header=None,sep='\s+',names=["rx","ry"])
# v =  pd.read_csv(path+"vel_0.980.txt" ,header=None,sep='\s+' , names=["vx","vy"])
# rinicial= pd.read_csv(path+"posiciones_init.txt",header=None,sep='\s+',names=["y","ry"])
# vinit= pd.read_csv(path+"velocidad_init.txt",header=None,sep='\s+', names=["vx","vy"])
# temp= pd.read_csv(path+"temperaturas_0.980.dat",header=0,sep='\s+' )
# tiempo= pd.read_csv(path+"tiemposdecol_0.980.txt",header=None,sep='\s+',names=["t"])


guardar_beta= False
guardar_temps = False

print(temp1)


# ! Variables

h =29
sigma =1.0
alpha = 0.9
vp = 0.001


T1 = temp1["T1_ave"] / (10*2*vp**2)
errorT1 = temp1["error"] /(10*2* vp**2)

T2 = 2*temp2["T2_ave"] / vp**2
errorT2 =2* temp2["error"] / vp**2 


z =  temp1["z"]

# print(T1)


print(z)


# # change default range so that new circles will work
figure, axis = plt.subplots(1, 2 ,figsize=(30,10) ) 

# For Sine Function
axis[0].errorbar(z, T1, yerr=errorT1,mfc="none",capsize=10,ms=12, color='C2',marker="o",linestyle="") 

axis[0].set_xlabel("Sine Function")
# For Cosine Function
axis[1].errorbar(z, T2, yerr=errorT2,mfc="none",capsize=10,ms=12, color='C2',marker="o",linestyle="") 
# # plt.plot(r["rx"][0:len(r["rx"])-1],r["ry"][0:len(r["rx"])-1], "o",markersize=10,c='#00FF00')

# plt.scatter( r["rx"][0:len(r["rx"])-1],r["ry"][0:len(r["rx"])-1], s=1200 ,  edgecolors='#00B000', facecolors='#00FF00' )  
# #? s = 1e5 para bolas de radio 1 
# # plt.scatter( r["rx"][0:len(r["rx"])-1],r["ry"][0:len(r["rx"])-1], s=points_radius**2 ,  facecolors='none', edgecolors='#00FF00' ) 
# # teta = np.linspace(0,2*np.pi,1000)
# # for (x,y) in zip(r['rx'],r['ry']):
# #     # cc = plt.Circle(( x , y ), 0.5 ) 
 
# #     # axes.set_aspect( 1 ) 
# #     xx =  (0.5*np.cos(teta)+x)/100
# #     yy= (0.5 *np.sin(teta)+y)/100
# #     axes.plot( xx, yy ) 
# # axes.axis('equal') # set aspect ratio to equal


# # plt.xlabel ( r' $x$ ', fontsize=40)
# # plt.ylabel ( r' $z$ ',rotation=0.0,labelpad=30,fontsize=40)
# # Turn off tick labels


# plt.xlim(-80,-40)
# plt.ylim(0,h+1)

# plt.xticks([])
# plt.yticks([])

# # plt.xticks(fontsize=25)
# # plt.yticks(fontsize=25)

# plt.tight_layout()



# plt.savefig ('tfield_h'+str(h)+'alfa_'+str(alfa)+'.pdf',format ='pdf')

plt.show()


