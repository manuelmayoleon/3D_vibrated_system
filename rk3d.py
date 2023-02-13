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
path2 = path+'/MeanTemp_vs_time/H29/'


# temp= pd.read_csv("ave_alpha_p90-error.dat",header=0,sep='\t\t+' )
temp= pd.read_csv(path2+"ave_alpha_p90_h4.dat",header=0,sep='\t\t+' )

# temp8sigma= pd.read_csv(path2+"ave_alpha_p90.dat",header=0,sep='\t\t+' )



temp.columns = temp.columns.str.strip()
# temp8sigma.columns = temp8sigma.columns.str.strip()

print(temp.columns)

n=512

sigma=1.0

# altura = 29
altura = 5

alfa = 0.90
# alfa = 0.750

w=0.001

rho = 0.02

txy_ini =  (temp["Tx"].iloc[0]+ temp["Ty"].iloc[0])

# txy_ini = temp["Txy"].iloc[0]
tz_ini = 2 * temp["Tz"].iloc[0] 





# txy_ini = (temp8sigma["Tx"].iloc[0] +temp8sigma["Ty"].iloc[0])
# tz_ini = 2* temp8sigma["Tz"].iloc[0]


# txy =   temp["Txy"][0:20]/txy_ini

txy =    (temp["Tx"][0:20] +temp["Ty"][0:20]) 



# txyerror =   temp["errorTxy"][0:20]/txy_ini

tz = 2 * temp["Tz"][0:20] 


print("beta sim")


print(tz[18]/txy[18]-1)

# tz =  2 * temp8sigma["Tz"][0:20]/txy_ini

# tzerror =  2 * temp["errorTz"][0:20]/txy_ini

time = temp['time'] [0:20]

# time = temp['# time'] [0:20]


# #?? resolucion de las ecuaciones en funcion del tiempo
x0= txy_ini
y0=tz_ini
a=0
b=time[len(time)-1]



h=0.01

# #??Resolucion de las ecuaciones en colpp
# # x0= temp['y'][0]/vp**2
# # y0=temp['z'][0]/vp**2
# # a=0
# # b=int(cols)



# # print(temp)
# # print(tiempo)

class Panel:
    def __init__(self,alfa=0.97,altura= 29,rho=0.02 ,w=0.001):
        self.rho = rho
        self.w=  w
        self.alfa = alfa
        self.altura = altura
    def f(self,t1,t2,t):
        
        return (1+self.alfa)*self.rho*np.sqrt(t1*np.pi)*( -4*(1-self.alfa)*t1 + ( 2 + 6 *self.alfa)*(t2 -t1 )/5)/3
 

    def g(self,t1,t2,t):
    
        return  2 * (1+self.alfa)*self.rho*np.sqrt(t1*np.pi)*( - 2*( 1 - self.alfa ) * t1 -( 17 - 9 * self.alfa ) * ( t2 - t1 ) / 5 ) / 3 +( 4.0*self.w*np.sqrt(2*t2)/self.altura ) * ( self.w / np.sqrt( np.pi ) + np.sqrt( t2 / 2 ) )

        #  return  2 * (1+self.alfa)*self.rho*np.sqrt(t1*np.pi)*( - 2*( 1 - self.alfa ) * t1 -( 17 - 9 * self.alfa ) * ( t2 - t1 ) / 5 ) / 3 +( 4.0*self.w*t2/self.altura ) 

        #  return  2 * (1+self.alfa)*self.rho*np.sqrt(t1*np.pi)*( - 2*( 1 - self.alfa ) * t1 -( 17 - 9 * self.alfa ) * ( t2 - t1 ) / 5 ) / 3 +( 4.0*self.w*t2/self.altura ) 
    
    def coupled(t,z,alfa,rho,w,altura):
    
        t1, t2 = z

        return [(1+alfa)*rho*np.sqrt(t1*np.pi)*( -4*(1-alfa)*t1 + ( 2 + 6 *alfa)*(t2 -t1 )/5)/3,2 * (1+alfa)*rho*np.sqrt(t1*np.pi)*( - 2*( 1 - alfa ) * t1 -( 17 - 9 * alfa ) * ( t2 - t1 ) / 5 ) / 3 +( 4.0*w*np.sqrt(2*t2)/altura ) * (w / np.sqrt( np.pi ) + np.sqrt( t2 / 2 )) ]



sol = solve_ivp(Panel.coupled, [a, b], [x0, y0], args=(alfa, rho, w, altura),
                    dense_output=True)

# print(sol.y[0])
# print(sol.x)
# print(sol.t)

t = np.linspace(0, b, 300)




z = sol.sol(t)




txysol= z[0]
tzsol= z[1] 




#     # return -4.0*epsilon**3.0*rho*np.sqrt(tx)*( ty -tx )/(3.0*np.sqrt(np.pi))
# #?? Representacion en funcion de las colisiones por particula
# # colpp=np.linspace(0,cols,len(temp["z"]))
# # plt.plot(colpp,temp["z"]/vp**2,color='C0',label="$T_z$ (MD)")
# # plt.plot(colpp,temp["y"]/vp**2,color='C1',label="$T_y$ (MD)")      
# # plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
# # plt.xlabel ( r' $s$ ', fontsize=30)
# # plt.ylabel ( r' $\tilde{T}$ ',rotation=0.0,fontsize=30)

fig22 = plt.subplots(figsize=(10,8))
timet0 =time*np.sqrt(txy_ini)



plt.plot(t*np.sqrt(txy_ini),txysol,color='C2',linestyle=":",label="$T_{xy}$ ")

plt.plot(t*np.sqrt(txy_ini),tzsol,color='C3',linestyle=":",label="$T_{z}$ ")


# #?? Representacion en funcion del tiempo
# plt.errorbar(timet0,txy,yerr=txyerror,color='C0',mfc="none",capsize=10,ms=12,marker="o",linestyle="",label="$T$ ")

plt.plot(timet0,txy,color='C2',marker="o",linestyle="",label="$T_{xy}$ (MD)")

plt.plot(timet0,tz,color='C3',marker="o",linestyle="",label="$T_{z}$ (MD)")


# plt.errorbar(timet0,tz,yerr=tzerror,color='C1',mfc="none",capsize=10,ms=12,marker="s",linestyle="",label="$T_z$ ")    
  
# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r' $t(T_0/m\sigma^2)^{1/2}$ ', fontsize=40)
# plt.ylabel ( r' $\frac{T}{T_0}$ ',rotation=0.0,labelpad=30,fontsize=40)

# plt.xlim(-2.5,340)


# plt.ylim(0,max(txy)+0.2)


plt.tick_params(axis='x', labelsize=25)
plt.tick_params(axis='y', labelsize=25)

plt.legend(loc=0,fontsize=22)


# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)

# plt.title ( r' \textbf {Solución de las ecuaciones para las temperaturas}  ',fontsize=40)


@jit(nopython=False)
def runge_kutta_system(f, g, x0, y0, a, b, h):
    t = np.arange(a, b + h, h)
    n = len(t)
    x = np.zeros(n)
    y = np.zeros(n)
    x[0] = x0
    y[0] = y0
    for i in range(n - 1):
        k1 = h * f(x[i], y[i], t[i])
        l1 = h * g(x[i], y[i], t[i])
        k2 = h * f(x[i] + k1 / 2, y[i] + l1 / 2, t[i] + h / 2)
        l2 = h * g(x[i] + k1 / 2, y[i] + l1 / 2, t[i] + h / 2)
        k3 = h * f(x[i] + k2 / 2, y[i] + l2 / 2, t[i] + h / 2)
        l3 = h * g(x[i] + k2 / 2, y[i] + l2 / 2, t[i] + h / 2)
        k4 = h * f(x[i] + k3, y[i] + l3, t[i] + h)
        l4 = h * g(x[i] + k3, y[i] + l3, t[i] + h)
        x[i + 1] = x[i] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + 2 * k4)
        y[i + 1] = y[i] + (1 / 6) * (l1 + 2 * l2 + 2 * l3 + 2 * l4)
    
    color=["#"+''.join([random.choice('0123456789ABCDEF') for i in range(6)])
       for j in range(2)]
    
    plt.plot(t*np.sqrt(x0), x ,color=color[0] , label ="txy"  )
    plt.plot(t*np.sqrt(x0), y ,color=color[1]  , label ="tz" )
    print(min(x))
    print("beta")
    print(y[-1]/x[-1]-1)
    plt.legend(loc=0,fontsize=30)
    # plt.show()
# # # np.seterr('raise')
# # # #?? Aplicación RK  a ecuacion en escala de colpp
# # # runge_kutta_system(Panel(alfa,epsilon,rho,vp).f_colpp,Panel(alfa,epsilon,rho,vp).g_colpp,x0,y0,a,b,h)
# # #?? Aplicación RK  a ecuacion en escala temporal 
# runge_kutta_system(Panel(alfa,30,rho,w).f,Panel(alfa,30,rho,w).g,x0,y0,a,b,h)
runge_kutta_system(func.Evol_eq(sigma,rho,w,alfa,altura,3).fpd_dim,func.Evol_eq(sigma,rho,w,alfa,altura,3).gpd_dim,x0,y0,a,b,h)


plt.tight_layout()

# plt.savefig('temp_dynamics.pdf', dpi = 1200)

plt.show()
