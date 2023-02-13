from numpy.core.function_base import linspace
import pandas as pd
import numpy as np
import scipy . stats as ss
import math
import matplotlib . mlab as mlab
import matplotlib . pyplot as plt
from scipy import optimize
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
from sympy import *
import functions as fun
from numba import jit




alfa= Symbol("alfa",positive=True)
h=Symbol("h",positive=True)
x = symbols('x')
np.seterr(divide='ignore', invalid='ignore')

d = 2

# #  ! d_dimension
# M=Matrix([ [- fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = d ).factor(d)*fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = d ).a11(alfa),-  fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = d ).factor(d)* fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = d ).a12(alfa) ],
#      [ - fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = d ).factor(d)* fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = d ).a21(alfa) ,- fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = d ).factor(d)*fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h).a22(alfa)] ] )
# # print(M)



# p = M.charpoly(x)

# print(p)


# poly_d = p.coeffs()[0] * x**2 +  p.coeffs()[1] * x +  p.coeffs()[2] 
# quadratic_equation = Eq(poly_d, 0)


# solucion=solve(quadratic_equation, x)

# print(solucion[0])
# print(solucion[1])

#  ! 3_dimension



# M=Matrix([ [ fun.hom_eigen3d(sigma=1.0,alfa=alfa, h  = h ).a11_3d(alfa),fun.hom_eigen3d(sigma=1.0,alfa=alfa, h  = h ).a12_3d(alfa) ],
#      [fun.hom_eigen3d(sigma=1.0,alfa=alfa, h  = h ).a21_3d(alfa) ,fun.hom_eigen3d(sigma=1.0,alfa=alfa, h  = h ).a22_3d(alfa)] ] )

# ?  M Matrix 
 

# print(M)

# ? Charasteristic polinomic of M
 
# p = M.charpoly(x)


# ? show  Charasteristic polinomic of M

# print(p)

# ? show  all coefficients of  Charasteristic polinomic

# print(p.coeffs())


# poly_3d = p.coeffs()[0] * x**2 +  p.coeffs()[1] * x +  p.coeffs()[2] 

# quadratic_equation = Eq(poly_3d, 0)


# solucion=solve(quadratic_equation, x)

# print(solucion[0])
# print(solucion[1])

# print( fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = 3 ).A(alfa))

# print( fun.hom_eigen3d(sigma=1.0,alfa=alfa, h  = h ).A_3d(alfa))


# print( fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = 3 ).B(alfa))

# print( fun.hom_eigen3d(sigma=1.0,alfa=alfa, h  = h ).B_3d(alfa))


# print( fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = 3 ).Ap(alfa))

# print( fun.hom_eigen3d(sigma=1.0,alfa=alfa, h  = h ).Ap_3d(alfa))


# print( fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = 3 ).Bp(alfa))

# print( fun.hom_eigen3d(sigma=1.0,alfa=alfa, h  = h ).Bp_3d(alfa))

# print( fun.hom_eigen(sigma=1.0,alfa=alfa, h  = h,d = 3 ).factor(3))
# print (4/15)
# print( fun.hom_eigen3d(sigma=1.0,alfa=alfa, h  = 8 ).betas(0.90))




def lamda1(alfa,h):
    
    return 1.33333333333333e-13*(1.664e+15*alfa**2*h**3 - 6.93200000000001e+15*alfa**2*h**2 + 9.547e+15*alfa**2*h - 4.34931249999999e+15*alfa**2 - 4.224e+15*alfa*h**3 + 1.7672e+16*alfa*h**2 - 2.4372e+16*alfa*h + 1.1064625e+16*alfa + 4.096e+15*h**3 - 1.8228e+16*h**2 + 2.6993e+16*h - 8.33483942856729e+16*np.sqrt(-0.000468151021503043*alfa**4*h**6 + 0.00440232947361771*alfa**4*h**5 - 0.0170497380360827*alfa**4*h**4 + 0.0348536424843862*alfa**4*h**3 - 0.0397044142721164*alfa**4*h**2 + 0.0239185038540117*alfa**4*h - 0.00595715407992565*alfa**4 + 0.00416736954657873*alfa**3*h**6 - 0.0386234910930032*alfa**3*h**5 + 0.148145762450567*alfa**3*h**4 - 0.301133059188179*alfa**3*h**3 + 0.342241907163804*alfa**3*h**2 - 0.206263954396548*alfa**3*h + 0.0515157559307716*alfa**3 - 0.0111059101272435*alfa**2*h**6 + 0.102906109726353*alfa**2*h**5 - 0.395393248871895*alfa**2*h**4 + 0.806516014154075*alfa**2*h**3 - 0.921278584878893*alfa**2*h**2 + 0.558864472780802*alfa**2*h - 0.140672850139514*alfa**2 + 0.011752123627958*alfa*h**6 - 0.109383574615815*alfa*h**5 + 0.422711324869806*alfa*h**4 - 0.868291213576237*alfa*h**3 + alfa*h**2 - 0.61232256864264*alfa*h + 0.155757948441184*alfa - 0.0040058161714001*h**6 + 0.0373873719285433*h**5 - 0.144962128679912*h**4 + 0.29890867737224*h**3 - 0.345737295156704*h**2 + 0.212714498047917*h - 0.0543904162226421+0j) - 1.33063125e+16)/(448.0*alfa*h**3 - 1744.0*alfa*h**2 + 2219.0*alfa*h - 923.0*alfa - 704.0*h**3 + 2832.0*h**2 - 3727.0*h + 1599.0)
def lamda2(alfa,h):            
   
    return 1.33333333333333e-13*(1.664e+15*alfa**2*h**3 - 6.93200000000001e+15*alfa**2*h**2 + 9.547e+15*alfa**2*h - 4.34931249999999e+15*alfa**2 - 4.224e+15*alfa*h**3 + 1.7672e+16*alfa*h**2 - 2.4372e+16*alfa*h + 1.1064625e+16*alfa + 4.096e+15*h**3 - 1.8228e+16*h**2 + 2.6993e+16*h + 8.33483942856729e+16*np.sqrt(-0.000468151021503043*alfa**4*h**6 + 0.00440232947361771*alfa**4*h**5 - 0.0170497380360827*alfa**4*h**4 + 0.0348536424843862*alfa**4*h**3 - 0.0397044142721164*alfa**4*h**2 + 0.0239185038540117*alfa**4*h - 0.00595715407992565*alfa**4 + 0.00416736954657873*alfa**3*h**6 - 0.0386234910930032*alfa**3*h**5 + 0.148145762450567*alfa**3*h**4 - 0.301133059188179*alfa**3*h**3 + 0.342241907163804*alfa**3*h**2 - 0.206263954396548*alfa**3*h + 0.0515157559307716*alfa**3 - 0.0111059101272435*alfa**2*h**6 + 0.102906109726353*alfa**2*h**5 - 0.395393248871895*alfa**2*h**4 + 0.806516014154075*alfa**2*h**3 - 0.921278584878893*alfa**2*h**2 + 0.558864472780802*alfa**2*h - 0.140672850139514*alfa**2 + 0.011752123627958*alfa*h**6 - 0.109383574615815*alfa*h**5 + 0.422711324869806*alfa*h**4 - 0.868291213576237*alfa*h**3 + alfa*h**2 - 0.61232256864264*alfa*h + 0.155757948441184*alfa - 0.0040058161714001*h**6 + 0.0373873719285433*h**5 - 0.144962128679912*h**4 + 0.29890867737224*h**3 - 0.345737295156704*h**2 + 0.212714498047917*h - 0.0543904162226421+0j) - 1.33063125e+16)/(448.0*alfa*h**3 - 1744.0*alfa*h**2 + 2219.0*alfa*h - 923.0*alfa - 704.0*h**3 + 2832.0*h**2 - 3727.0*h + 1599.0)
def lamda1_3d(alfa,h):
    
    return 1.66666666666667e-14*(1.3312e+16*alfa**2*h**3 - 5.54559999999999e+16*alfa**2*h**2 + 7.63759999999999e+16*alfa**2*h - 3.47945e+16*alfa**2 - 3.3792e+16*alfa*h**3 + 1.41376e+17*alfa*h**2 - 1.94976e+17*alfa*h + 8.85170000000001e+16*alfa + 3.2768e+16*h**3 - 1.45824e+17*h**2 + 2.15944e+17*h - 6.66787154285384e+17*np.sqrt(-0.000468151021503043*alfa**4*h**6 + 0.00440232947361774*alfa**4*h**5 - 0.0170497380360828*alfa**4*h**4 + 0.0348536424843863*alfa**4*h**3 - 0.0397044142721164*alfa**4*h**2 + 0.0239185038540117*alfa**4*h - 0.00595715407992568*alfa**4 + 0.00416736954657873*alfa**3*h**6 - 0.0386234910930031*alfa**3*h**5 + 0.148145762450567*alfa**3*h**4 - 0.30113305918818*alfa**3*h**3 + 0.342241907163805*alfa**3*h**2 - 0.206263954396548*alfa**3*h + 0.0515157559307717*alfa**3 - 0.0111059101272435*alfa**2*h**6 + 0.102906109726353*alfa**2*h**5 - 0.395393248871895*alfa**2*h**4 + 0.806516014154076*alfa**2*h**3 - 0.921278584878895*alfa**2*h**2 + 0.558864472780804*alfa**2*h - 0.140672850139514*alfa**2 + 0.011752123627958*alfa*h**6 - 0.109383574615815*alfa*h**5 + 0.422711324869807*alfa*h**4 - 0.868291213576237*alfa*h**3 + alfa*h**2 - 0.612322568642641*alfa*h + 0.155757948441184*alfa - 0.00400581617140009*h**6 + 0.0373873719285433*h**5 - 0.144962128679913*h**4 + 0.298908677372241*h**3 - 0.345737295156704*h**2 + 0.212714498047916*h - 0.0543904162226421+0j) - 1.064505e+17)/(448.0*alfa*h**3 - 1744.0*alfa*h**2 + 2219.0*alfa*h - 923.0*alfa - 704.0*h**3 + 2832.0*h**2 - 3727.0*h + 1599.0)
def lamda2_3d(alfa,h):
   
    return 1.66666666666667e-14*(1.3312e+16*alfa**2*h**3 - 5.54559999999999e+16*alfa**2*h**2 + 7.63759999999999e+16*alfa**2*h - 3.47945e+16*alfa**2 - 3.3792e+16*alfa*h**3 + 1.41376e+17*alfa*h**2 - 1.94976e+17*alfa*h + 8.85170000000001e+16*alfa + 3.2768e+16*h**3 - 1.45824e+17*h**2 + 2.15944e+17*h + 6.66787154285384e+17*np.sqrt(-0.000468151021503043*alfa**4*h**6 + 0.00440232947361774*alfa**4*h**5 - 0.0170497380360828*alfa**4*h**4 + 0.0348536424843863*alfa**4*h**3 - 0.0397044142721164*alfa**4*h**2 + 0.0239185038540117*alfa**4*h - 0.00595715407992568*alfa**4 + 0.00416736954657873*alfa**3*h**6 - 0.0386234910930031*alfa**3*h**5 + 0.148145762450567*alfa**3*h**4 - 0.30113305918818*alfa**3*h**3 + 0.342241907163805*alfa**3*h**2 - 0.206263954396548*alfa**3*h + 0.0515157559307717*alfa**3 - 0.0111059101272435*alfa**2*h**6 + 0.102906109726353*alfa**2*h**5 - 0.395393248871895*alfa**2*h**4 + 0.806516014154076*alfa**2*h**3 - 0.921278584878895*alfa**2*h**2 + 0.558864472780804*alfa**2*h - 0.140672850139514*alfa**2 + 0.011752123627958*alfa*h**6 - 0.109383574615815*alfa*h**5 + 0.422711324869807*alfa*h**4 - 0.868291213576237*alfa*h**3 + alfa*h**2 - 0.612322568642641*alfa*h + 0.155757948441184*alfa - 0.00400581617140009*h**6 + 0.0373873719285433*h**5 - 0.144962128679913*h**4 + 0.298908677372241*h**3 - 0.345737295156704*h**2 + 0.212714498047916*h - 0.0543904162226421+0j) - 1.064505e+17)/(448.0*alfa*h**3 - 1744.0*alfa*h**2 + 2219.0*alfa*h - 923.0*alfa - 704.0*h**3 + 2832.0*h**2 - 3727.0*h + 1599.0)
def lamda1_2d(alfa,h):
    
    return 4.0e-13*(1.19380520836412e+15*alfa**2*h**3 - 4.97963972051546e+15*alfa**2*h**2 + 6.86542371606075e+15*alfa**2*h - 3.12958920390942e+15*alfa**2 - 3.0159289474462e+15*alfa*h**3 + 1.26122350331848e+16*alfa*h**2 - 1.73585084536063e+16*alfa*h + 7.8447913399881e+15*alfa + 3.33008821280518e+15*h**3 - 1.50110961900438e+16*h**2 + 2.25273065237373e+16*h - 8.04815498329378e+16*np.sqrt(-0.000919111957516917*alfa**4*h**6 + 0.00850207678793603*alfa**4*h**5 - 0.0325343020445037*alfa**4*h**4 + 0.0659463842669901*alfa**4*h**3 - 0.0747044063789934*alfa**4*h**2 + 0.0448563213590659*alfa**4*h - 0.0111569971096413*alfa**4 + 0.00547810230382099*alfa**3*h**6 - 0.0509948719708735*alfa**3*h**5 + 0.196587060942328*alfa**3*h**4 - 0.401836648317936*alfa**3*h**3 + 0.459458279768645*alfa**3*h**2 - 0.278691895349622*alfa**3*h + 0.0700763821001939*alfa**3 - 0.0115705712211546*alfa**2*h**6 + 0.108388969407236*alfa**2*h**5 - 0.420832425382051*alfa**2*h**4 + 0.867045160171535*alfa**2*h**3 - alfa**2*h**2 + 0.612273468235399*alfa**2*h - 0.155505760566657*alfa**2 + 0.0106368314871255*alfa*h**6 - 0.100299286928392*alfa*h**5 + 0.39230905223747*alfa*h**4 - 0.814915454843055*alfa*h**3 + 0.948340749234587*alfa*h**2 - 0.586331980161359*alfa*h + 0.150493694278258*alfa - 0.00343021359211221*h**6 + 0.0324888490607037*h**5 - 0.127700718105717*h**4 + 0.266684605327385*h**3 - 0.312143056533324*h**2 + 0.194183428450345*h - 0.0501687962111496+0j) - 1.12576577780142e+16)/(448.0*alfa*h**3 - 1744.0*alfa*h**2 + 2219.0*alfa*h - 923.0*alfa - 704.0*h**3 + 2832.0*h**2 - 3727.0*h + 1599.0)


def lamda2_2d(alfa,h):
   
    return 4.0e-13*(1.19380520836412e+15*alfa**2*h**3 - 4.97963972051546e+15*alfa**2*h**2 + 6.86542371606075e+15*alfa**2*h - 3.12958920390942e+15*alfa**2 - 3.0159289474462e+15*alfa*h**3 + 1.26122350331848e+16*alfa*h**2 - 1.73585084536063e+16*alfa*h + 7.8447913399881e+15*alfa + 3.33008821280518e+15*h**3 - 1.50110961900438e+16*h**2 + 2.25273065237373e+16*h + 8.04815498329378e+16*np.sqrt(-0.000919111957516917*alfa**4*h**6 + 0.00850207678793603*alfa**4*h**5 - 0.0325343020445037*alfa**4*h**4 + 0.0659463842669901*alfa**4*h**3 - 0.0747044063789934*alfa**4*h**2 + 0.0448563213590659*alfa**4*h - 0.0111569971096413*alfa**4 + 0.00547810230382099*alfa**3*h**6 - 0.0509948719708735*alfa**3*h**5 + 0.196587060942328*alfa**3*h**4 - 0.401836648317936*alfa**3*h**3 + 0.459458279768645*alfa**3*h**2 - 0.278691895349622*alfa**3*h + 0.0700763821001939*alfa**3 - 0.0115705712211546*alfa**2*h**6 + 0.108388969407236*alfa**2*h**5 - 0.420832425382051*alfa**2*h**4 + 0.867045160171535*alfa**2*h**3 - alfa**2*h**2 + 0.612273468235399*alfa**2*h - 0.155505760566657*alfa**2 + 0.0106368314871255*alfa*h**6 - 0.100299286928392*alfa*h**5 + 0.39230905223747*alfa*h**4 - 0.814915454843055*alfa*h**3 + 0.948340749234587*alfa*h**2 - 0.586331980161359*alfa*h + 0.150493694278258*alfa - 0.00343021359211221*h**6 + 0.0324888490607037*h**5 - 0.127700718105717*h**4 + 0.266684605327385*h**3 - 0.312143056533324*h**2 + 0.194183428450345*h - 0.0501687962111496+0j) - 1.12576577780142e+16)/(448.0*alfa*h**3 - 1744.0*alfa*h**2 + 2219.0*alfa*h - 923.0*alfa - 704.0*h**3 + 2832.0*h**2 - 3727.0*h + 1599.0)


alfa=np.linspace(0.40,0.99,500)

hh=4

fig, ax = plt.subplots(1,1,figsize=(14,8))

# plt.plot(alfa,lamda1(alfa,h),color="C0",label="$\lambda_1, h=4 \sigma$")
# plt.plot(alfa,lamda2(alfa,h),color="C0",linestyle="--",label="$\lambda_2  ,  h=4 \sigma $")


# plt.plot(alfa,np.real(lamda1(alfa,h)),color="C1",label="$\lambda_1, h=4 \sigma$")
# plt.plot(alfa,np.real(lamda2(alfa,h)),color="C1",linestyle="--",label="$\lambda_2  ,  h=4 \sigma $")



plt.plot(alfa,np.real(lamda1_3d(alfa,2)),color="C0",linewidth= 1.5,label="$\mathcal{ R } (\lambda_1), h=2 \sigma$")
plt.plot(alfa,np.real(lamda2_3d(alfa,2)),color="C0",linewidth= 1.5,linestyle="--",label=r"$\mathcal{ R}(\lambda_2)  ,  h=2 \sigma $")


plt.plot(alfa,np.imag(lamda1(alfa,2)),alpha=0.5,linewidth= 1.5,color="C0",label=r"$\mathcal{ I}(\lambda_1), h=2 \sigma$")
plt.plot(alfa,np.imag(lamda2(alfa,2)),alpha=0.5,linewidth= 1.5,color="C0",linestyle="--",label=r"$\mathcal{ I}(\lambda_2), h=2 \sigma$")



plt.plot(alfa,np.real(lamda1_3d(alfa,4)),color="C1",linewidth= 1.5,label="$\mathcal{ R } (\lambda_1), h=4 \sigma$")
plt.plot(alfa,np.real(lamda2_3d(alfa,4)),color="C1",linewidth= 1.5,linestyle="--",label=r"$\mathcal{ R}(\lambda_2)  ,  h=4 \sigma $")


plt.plot(alfa,np.imag(lamda1(alfa,4)),alpha=0.5,linewidth= 1.5,color="C1",label=r"$\mathcal{ I}(\lambda_1), h=4 \sigma$")
plt.plot(alfa,np.imag(lamda2(alfa,4)),alpha=0.5,linewidth= 1.5,color="C1",linestyle="--",label=r"$\mathcal{ I}(\lambda_2), h=4 \sigma$")


plt.plot(alfa,np.real(lamda1_3d(alfa,8)),color="C2",linewidth= 1.5,label="$\mathcal{ R } (\lambda_1), h=8 \sigma$")
plt.plot(alfa,np.real(lamda2_3d(alfa,8)),color="C2",linewidth= 1.5,linestyle="--",label=r"$\mathcal{ R}(\lambda_2)  ,  h=8 \sigma $")


plt.plot(alfa,np.imag(lamda1(alfa,8)),alpha=0.5,linewidth= 1.5,color="C2",label=r"$\mathcal{ I}(\lambda_1), h=8 \sigma$")
plt.plot(alfa,np.imag(lamda2(alfa,8)),alpha=0.5,linewidth= 1.5,color="C2",linestyle="--",label=r"$\mathcal{ I}(\lambda_2), h=8\sigma$")



plt.plot(alfa,np.real(lamda1_3d(alfa,29)),color="C3",linewidth= 1.5,label="$\mathcal{ R } (\lambda_1), h=29 \sigma$")
plt.plot(alfa,np.real(lamda2_3d(alfa,29)),color="C3",linewidth= 1.5,linestyle="--",label=r"$\mathcal{ R}(\lambda_2)  ,  h=29 \sigma $")


plt.plot(alfa,np.imag(lamda1(alfa,29)),alpha=0.5,linewidth= 1.5,color="C3",label=r"$\mathcal{ I}(\lambda_1), h=29 \sigma$")
plt.plot(alfa,np.imag(lamda2(alfa,29)),alpha=0.5,linewidth= 1.5,color="C3",linestyle="--",label=r"$\mathcal{ I}(\lambda_2), h=29\sigma$")





# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r' $\alpha$ ', fontsize=30)
plt.ylabel ( r'  ',fontsize=30)
plt.legend(loc=0,fontsize=30)
plt.title ( r' \textbf {Autovalores de la matriz $M$ $(d=3)$}  ',fontsize=40)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)



# height = np.linspace(2,29,1000)


# altura = list()
# alpha = list()


# for i in height:
#     count=0
#     for j in alfa:
#         if (np.imag(lamda1(j,i))==0 and count == 0):
#             altura.append(i)
#             alpha.append(j)
#             count+=1

# fig22, ax22 = plt.subplots(1,1,figsize=(14,8))


# print(alpha[len(alpha)-1])

# print(altura[len(altura)-1])


# plt.plot(alpha,altura,alpha=0.5,linewidth= 1.5,color="C3")






# # plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
# plt.xlabel ( r' $\alpha$ ', fontsize=30)
# plt.ylabel ( r'  ',fontsize=30)
# # plt.legend(loc=0,fontsize=30)
# plt.title ( r' \textbf {Dependencia entre h y $\alpha$ $(d=3)$}  ',fontsize=40)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)


fig2, ax2 = plt.subplots(1,1,figsize=(14,8))

# plt.plot(alfa,lamda1(alfa,h),color="C0",label="$\lambda_1, h=4 \sigma$")
# plt.plot(alfa,lamda2(alfa,h),color="C0",linestyle="--",label="$\lambda_2  ,  h=4 \sigma $")


# plt.plot(alfa,np.real(lamda1(alfa,h)),color="C1",label="$\lambda_1, h=4 \sigma$")
# plt.plot(alfa,np.real(lamda2(alfa,h)),color="C1",linestyle="--",label="$\lambda_2  ,  h=4 \sigma $")



plt.plot(alfa,np.real(lamda1_2d(alfa,2)),color="C0",linewidth= 1.5,label="$\mathcal{ R } (\lambda_1), h=2 \sigma$")
plt.plot(alfa,np.real(lamda2_2d(alfa,2)),color="C0",linewidth= 1.5,linestyle="--",label=r"$\mathcal{ R}(\lambda_2)  ,  h=2 \sigma $")


plt.plot(alfa,np.imag(lamda1_2d(alfa,2)),alpha=0.5,linewidth= 1.5,color="C0",label=r"$\mathcal{ I}(\lambda_1), h=2 \sigma$")
plt.plot(alfa,np.imag(lamda2_2d(alfa,2)),alpha=0.5,linewidth= 1.5,color="C0",linestyle="--",label=r"$\mathcal{ I}(\lambda_2), h=2 \sigma$")



plt.plot(alfa,np.real(lamda1_2d(alfa,4)),color="C1",linewidth= 1.5,label="$\mathcal{ R } (\lambda_1), h=4 \sigma$")
plt.plot(alfa,np.real(lamda2_2d(alfa,4)),color="C1",linewidth= 1.5,linestyle="--",label=r"$\mathcal{ R}(\lambda_2)  ,  h=4 \sigma $")


plt.plot(alfa,np.imag(lamda1_2d(alfa,4)),alpha=0.5,linewidth= 1.5,color="C1",label=r"$\mathcal{ I}(\lambda_1), h=4 \sigma$")
plt.plot(alfa,np.imag(lamda2_2d(alfa,4)),alpha=0.5,linewidth= 1.5,color="C1",linestyle="--",label=r"$\mathcal{ I}(\lambda_2), h=4 \sigma$")


plt.plot(alfa,np.real(lamda1_2d(alfa,8)),color="C2",linewidth= 1.5,label="$\mathcal{ R } (\lambda_1), h=8 \sigma$")
plt.plot(alfa,np.real(lamda2_2d(alfa,8)),color="C2",linewidth= 1.5,linestyle="--",label=r"$\mathcal{ R}(\lambda_2)  ,  h=8 \sigma $")


plt.plot(alfa,np.imag(lamda1_2d(alfa,8)),alpha=0.5,linewidth= 1.5,color="C2",label=r"$\mathcal{ I}(\lambda_1), h=8 \sigma$")
plt.plot(alfa,np.imag(lamda2_2d(alfa,8)),alpha=0.5,linewidth= 1.5,color="C2",linestyle="--",label=r"$\mathcal{ I}(\lambda_2), h=8\sigma$")



plt.plot(alfa,np.real(lamda1_2d(alfa,29)),color="C3",linewidth= 1.5,label="$\mathcal{ R } (\lambda_1), h=29 \sigma$")
plt.plot(alfa,np.real(lamda2_2d(alfa,29)),color="C3",linewidth= 1.5,linestyle="--",label=r"$\mathcal{ R}(\lambda_2)  ,  h=29 \sigma $")


plt.plot(alfa,np.imag(lamda1_2d(alfa,29)),alpha=0.5,linewidth= 1.5,color="C3",label=r"$\mathcal{ I}(\lambda_1), h=29 \sigma$")
plt.plot(alfa,np.imag(lamda2_2d(alfa,29)),alpha=0.5,linewidth= 1.5,color="C3",linestyle="--",label=r"$\mathcal{ I}(\lambda_2), h=29\sigma$")



# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r' $\alpha$ ', fontsize=30)
plt.ylabel ( r'  ',fontsize=30)
# plt.legend(loc=0,fontsize=30)
plt.title ( r' \textbf {Autovalores de la matriz $M$ $(d=2)$}  ',fontsize=40)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

fig23, ax23 = plt.subplots(1,1,figsize=(14,8))

# plt.plot(alfa,lamda1(alfa,h),color="C0",label="$\lambda_1, h=4 \sigma$")
# plt.plot(alfa,lamda2(alfa,h),color="C0",linestyle="--",label="$\lambda_2  ,  h=4 \sigma $")


# plt.plot(alfa,np.real(lamda1(alfa,h)),color="C1",label="$\lambda_1, h=4 \sigma$")
# plt.plot(alfa,np.real(lamda2(alfa,h)),color="C1",linestyle="--",label="$\lambda_2  ,  h=4 \sigma $")

time = np.linspace(0,40,1000)

plt.plot(time,np.exp(lamda1_3d(0.6,29)*time)*np.cos(np.imag(lamda1_3d(0.6,29))*time)+np.exp(lamda2_3d(0.6,29)*time)*np.cos(np.imag(lamda2_3d(0.6,29))*time),color="C0",linewidth= 1.5,label="$\alpha = 0.4, h=4 \sigma$")

# plt.plot(time,np.exp(lamda1_3d(0.8,4)*time)+np.exp(lamda2_3d(0.8,4)*time),color="C1",linewidth= 1.5,label="$\alpha = 0.8, h=4 \sigma$")

# plt.plot(alfa,np.real(lamda2_2d(alfa,2)),color="C0",linewidth= 1.5,linestyle="--",label=r"$\mathcal{ R}(\lambda_2)  ,  h=29 \sigma $")





# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r' $t$ ', fontsize=30)
plt.ylabel ( r'  ',fontsize=30)
# plt.legend(loc=0,fontsize=30)
# plt.title ( r' \textbf {Autovalores de la matriz $M$ $(d=2)$}  ',fontsize=40)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)







plt.tight_layout()
 
 
plt.show()
