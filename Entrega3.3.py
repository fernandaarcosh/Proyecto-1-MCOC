from matplotlib.pylab import *
from numpy import random

L=1.       #Largo dominio
n=100      #n° intervalos
dx= L / n
dt=1.      #s
#vector x
x=linspace(0, L, n+1)

#condicion inicial
def fun_u0(x):
    return 10*exp(-(x-0.5)**2/0.1**2)

u0 = fun_u0(x)

#Creacion de listas de temperaturas aleatoreas entre 0 y 20 para las condiciones de bordes iniciales y finales
u_ki=range(5)   #Condicion de borde inicial o en x=0
u_kf=range(5)   #Condicion de borde final o en x=n   
for i in range(len(u_ki)):
    u_ki[i]=random.randint(0,20)
    u_kf[i]=random.randint(0,20)

#Valores usando Hierro
K=79.5    #W/(m*K)
c=450.    #J/(Kg*K)
rho=7800. #Kg/m^3

Alpha=K*dt/(c*rho*dx**2) 

plot(x,u0, "k--")
#Loop para usar los distintos pares de condiciones de bordes
for j in range(len(u_ki)):
    #vector U en un tiempo t = k
    u_k = u0.copy()
    #vector U en un tiempo t = k+1
    u_km1=u_k.copy()
    #Loop en el tiempo
    k=0
    for k in range(1000):
        t=dt*k
        #Definicion de las condiciones de bordes
        u_k[0]=u_ki[j]
        u_k[n]=u_kf[j]
        #Loop en el espacio
        for i in range(1,n): 
            u_km1[i] = u_k[i] + Alpha*(u_k[i+1]-2*u_k[i]+u_k[i-1])
        u_k = u_km1
        
        if k % 50 == 0:
            plot(x,u_k)
            axis([0, 1, 0, 20])
        title('Para u_ki = {} C, u_kf = {} C, t = {} s'.format(u_ki[j], u_kf[j], k*dt))
        xlabel('Largo barra (m)')
        ylabel('Temperatura (C)')
    show()