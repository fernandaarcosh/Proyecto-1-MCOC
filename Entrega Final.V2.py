from matplotlib.pylab import *
from scipy import interpolate
import numpy
import Datos_Amb    #Se utiliza otro archivo python para optimizar la velocidad de este archivo

#Por algun motivo, funciona sin problemas en python 2.7 pero no en python 3.7

#Primero se definen todas las funciones a utilizar
#Fórmula para el calculo del q(t), conociendo el grado de hidratacion
def Q(t):
    if t==0:
        return 0
    else:
        tau = 9.8
        b = 0.98
        au = 0.7
        E = 27.1
        R = 8.31
        tr = 20
        tc = 27
        return 38591.46031*0.400463*(tau/(t*60))**b*(b/(t*60))*au*np.exp(-(tau/(t*60))**b)*np.exp(E/R*(1/(273+tr)- 1/(273+tc)))


#En caso de requerir este metodo en vez de utilizar el llamado de listado con numeros decimales (probablemente no se utilizara)
#==============================================================================
# def promedioX(u,i,j,Z):
#     return (u[i,j,Z] + u[i+1,j,Z])/2
# def promedioY(u,i,j,Z):
#     return (u[i,j,Z] + u[i,j+1,Z])/2
# def promedioZ(u,i,j,Z):
#     return (u[i,j,Z] + u[i,j,Z+1])/2
#=======================================
# Dimensiones del bloque de hormigon

a = 0.5        #Alto del dominio y en m.
b = 1.         #Largo del dominio x en m.
C = 0.57       #Ancho del dominio z en m.
Nx = 25        #Numero de intervalos en x
Ny = 12.5        #Numero de intervalos en Y
Nz = 14.25        #Numero de intervalos en Z
#==============================================================================
# Se utilizaron estos Nx, Ny y Nz para localizar facilmente la separación de 4 cm de los bordes, 
# al mismo tiempo optimizar el tiempo de actividad del programa y mantener la igualdad de dx, dy y dz.
#==============================================================================
dx = b / Nx     #Discretizacion espacial en X
dy = a / Ny     #Discretizacion espacial en Y
dz = C / Nz     #Discretizacion espacial en Z

h = dx    # = dy = dz 

if dx != dy:
    print("ERRROR!!!!! dx != dy")
    exit(-1)   #-1 le dice al SO que el programa fallo.....


#Funcion de conveniencia para calcular coordenadas del punto (i,j,Z)
 
coords = lambda i, j, Z : (dx*i, dy*j, dz*Z)
x, y, z = coords(1,1,1) 
 
#print "x = ", x
#print "y = ", y
#print "z = ", z
 
u_k = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)
u_km1 = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)

#Parametros del problema (hormigon)
dt = 1.0        # min
K = 9.87        # m^2 / s   
c = 0.983       # J / kg C
rho = 2425.0    # kg / m^3
#HT = 0.579      # kj / g
#alpha = K*dt/(c*rho*dx**2)

alpha_bueno = 0.0001
#dt = alpha_bueno*(c*rho*dx**2)/K
#alpha = K*dt/(c*rho*dx**2) 
 
#Informar cosas interesantes
#print "dt = ", dt
#print "dx = ", dx
#print "K = ", K
#print "c = ", c
#print "rho = ", rho
#print "alpha = ", alpha

#Se denominaron los puntos en base a foto entregada 

Punto1 = []         #Costado arriba 1
Punto2 = []         #Costado medio 2
Punto3 = []         #Costado bajo 3
Punto4 = []         #Centro arriba 4
Punto5 = []         #Centro medio 5
Punto6 = []         #Centro bajo 6
Punto7 = []         #Esquina arriba 7
Punto8 = []         #Esquina medio 8
Punto9 = []         #Esquina bajo 9
time=[]             #Graficar el tiempo en Dias, en vez de minutos o segundos
k = 0
u_k[:,:,:]=20
#Contadores para distribuir el tiempo transcurrido en Dias, Hrs., Min. 
Hora = 0
#Loop en el tiempo 
for k in range(len(ejeTiempo)):  #se usa ese rango para equiparar el largo de la lista de los valores de la temperatura ambiente entregados en el archivo excel
    t = dt*k              #Por como se definio el tiempo, se avanza por cada minuto en vez de cada segundo.
    if k%60 == 0 and k>0:     #Contador para especificar bien la hora y que no supere las 24 horas.
        Hora +=1
    if Hora == 24:
        Hora = 0
    #print "k = ", k, "t = ", t
    
    #CB esencial
    u_k[:,:,-1] = Tambiente[k]
    #print "temperatura ambiente:", Tambiente[k]
    #Loop en el espacio   i = 1 ... Nx-1, j = 1 ... Ny-1, Z = 1 ... Nz-1
    for i in range(1,int(Nx)):
        for j in range(1,int(Ny)):
            for Z in range(1,int(Nz)):
                #Algoritmo de diferencias finitas 3-D para difusion
 
                #Laplaciano 3-D
                nabla_u_k = (u_k[i-1,j,Z] + u_k[i+1,j,Z] + u_k[i,j-1,Z] + u_k[i,j+1,Z] + u_k[i,j,Z-1] + u_k[i,j, Z+1] - 6*u_k[i,j,Z])/h**2

                #Forward euler..
                u_km1[i,j,Z] = u_k[i,j,Z] + alpha_bueno*nabla_u_k + Q(t)
 
    #CB natural
    u_km1[0,:,:] = u_km1[1,:,:]             #du/dx(0,y,z,t)=0    
    u_km1[:,0,:] = u_km1[:,1,:]             #du/dy(x,0,,z,t)=0   
    u_km1[:,:,0] = u_km1[:,:,1]             #du/dx(x,y,0,t)=0    
    u_km1[-1,:,:] = u_km1[-2,:,:]           #du/dx(x,y,C,t)=0
    u_km1[:,-1,:] = u_km1[:,-2,:]           #du/dx(x,b,z,t)=0
    #Avanzar la solucion a k + 1
    u_k = u_km1
    #Por cada coordenada que se avanza en el eje x, y o z, se avanza cada 0.04, es decir: u_k[1,1,1] esta en la posicion x=0.04, y=0.04, z=0.04
    Punto_1 = u_k[12.5,11.5,1]      #punto 1
    Punto_2 = u_k[12.5,6.25,1]      #punto 2
    Punto_3 = u_k[12.5,1,1]         #punto 3
    Punto_4 = u_k[12.5,11.5,7.125]  #punto 4
    Punto_5 = u_k[12.5,6.25,7.125]  #punto 5
    Punto_6 = u_k[12.5,1,7.125]     #punto 6
    Punto_7 = u_k[24,11.5,1]        #punto 7
    Punto_8 = u_k[24,6.25,1]        #punto 8
    Punto_9 = u_k[24,1,1]           #punto 9
    
#==============================================================================
    #print "Punto 1: {}".format(Punto_1)
    #print "Punto 2: {}".format(Punto_2)
    #print "Punto 3: {}".format(Punto_3)
    #print "Punto 4: {}".format(Punto_4)
    #print "Punto 5: {}".format(Punto_5)
    #print "Punto 6: {}".format(Punto_6)
    #print "Punto 7: {}".format(Punto_7)
    #print "Punto 8: {}".format(Punto_8)
    #print "Punto 9: {}".format(Punto_9)
#==============================================================================

    print "Tmax = ", u_k.max()
    print ("tiempo transcurrido: {} Dias, {} Hrs., {} Mins.".format(((k/24)/60)%60, Hora, t%60))
    
    time.append(t/60./24.)
    Punto1.append(Punto_1)
    Punto2.append(Punto_2)
    Punto3.append(Punto_3)
    Punto4.append(Punto_4)
    Punto5.append(Punto_5)
    Punto6.append(Punto_6)
    Punto7.append(Punto_7)
    Punto8.append(Punto_8)
    Punto9.append(Punto_9)   

t=linspace(0,12,10000)
plot(t,tempAmb(t))
plot(time, Punto1, 'b')
plot(time, Punto2, 'r')
plot(time, Punto3, 'm')
plot(time, Punto4, 'k')
plot(time, Punto5, 'gray')
plot(time, Punto6, 'c')
plot(time, Punto7, 'y')
plot(time, Punto8, 'g')
plot(time, Punto9, 'violet')
xlabel('Tiempo (Dias)')
ylabel('Temperatura (C)')
legend(['T Ambiente','Costado arriba 1','Costado medio 2','Costado bajo 3','Centro arriba 4','Centro medio 5','Centro bajo 6','Esquina arriba 7','Esquina medio 8','Esquina bajo 9'], loc='best')
show()