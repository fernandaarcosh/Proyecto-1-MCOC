from matplotlib.pylab import *

#Scipy.interpolate
#Caso en donde u(0,y,t)=u(x,0,t)=20 y du/dx(b,y,t)=du/dy(x,a,t)=0
  
a = 1.        #Alto del dominio y en cm.
b = 1.         #Largo del dominio x en cm.
C = 1.       #Ancho del dominio z en cm.
Nx = 3      #Numero de intervalos en x
Ny = 3       #Numero de intervalos en Y
Nz = 3       #Numero de intervalos en Z
#==============================================================================
# Se utilizan esos Nx, Ny y Nz para poder localizar con mayor facilidad los puntos 0.04 de separacion de los sensores 
#==============================================================================
dx = b / Nx     #Discretizacion espacial en X
dy = a / Ny     #Discretizacion espacial en Y
dz = C / Nz     #Discretizacion espacial en Z
 
h = dx    # = dy = dz 

#if dx != dy:
#    print("ERRROR!!!!! dx != dy")
#    exit(-1)   #-1 le dice al SO que el programa fallo.....
 
#Funcion de conveniencia para calcular coordenadas del punto (i,j)
 
coords = lambda i, j, Z : (dx*i, dy*j, dz*Z)
x, y, z = coords(1,2,3) 
 
print "x = ", x
print "y = ", y
print "z = ", z
 
u_k = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)
u_km1 = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)
 
#Buena idea definir funciones que hagan el codigo expresivo
def printbien(u):
    print u.T[Nx::-1,:]

#Con la separaciòn actual de Nz, no se lgora llegar al punto medio en z especificamente, por lo cual se hace un promedio entre el punto anterior y siguiente al punto medio exacto
def promedio(u,i,j,Z):
    return (u[i,j,Z] + u[i,j,Z+1])/2
    
#print u_k               
#printbien(u_k)                  #Imprime con el eje y invertido
 
def imshowbien(u):
    imshow(u.T[Nx::-1,:])
    colorbar(extend='both',cmap='plasma')
    clim(10, 30)
 
#Parametros del problema (hierro)
dt = 1.0       # s
K = 79.5       # m^2 / s   
c = 450.       # J / kg C
rho = 7800.    # kg / m^3
alpha = K*dt/(c*rho*dx**2)

alpha_bueno = 0.0001
dt = alpha_bueno*(c*rho*dx**2)/K
alpha = K*dt/(c*rho*dx**2) 
 
#Informar cosas interesantes
print "dt = ", dt
print "dx = ", dx
print "K = ", K
print "c = ", c
print "rho = ", rho
print "alpha = ", alpha

#Puntos
#==============================================================================
#u_k[50,46,4] #punto 1
#u_k[50,25,4] #punto 2
#u_k[50,4,4]  #punto 3
# promedio(u_k[50,46,28]) #punto 4
# promedio(u_k[50,25,28]) #punto 5
# promedio(u_k[50,4,28])  #punto 6
# u_k[96,46,4] #punto 7
# u_k[96,25,4] #punto 8
# u_k[96,4,4]  #punto 9
#==============================================================================

time=[]
m=[]
k = 0
u_k[:,:,:]=20. 
#print u_k
Hora = 0
Dias = 0
#Loop en el tiempo 
dnext_t = 60   #  20.00
next_t = 0.
framenum = 0
for k in range(5000):    
    t = dt*(k+1)
    if k%60 == 0 and k>0:
        Hora +=1
    if Hora == 24:
        Hora = 0
        Dias += 1
    print "k = ", k, "t = ", t
 
    #CB esencial
    u_k[0,:,:] = u_k[1,:,:]             #du/dx(0,y,z,t)=0    Cara lateral 1
    u_k[:,0,:] = u_k[:,1,:]             #du/dy(x,0,,z,t)=0   Cara abajo
    u_k[:,:,0] = u_k[:,:,1]             #du/dx(x,y,0,t)=0    Caras del ancho 1
    u_k[:,:,-1] = u_k[:,:,-2]           #du/dy(x,y,,C,t)=0   Cara del ancho 2
    u_k[-1,:,:] = 20 + 10*sin(2*pi/(24.0/(60./60./60.))*t)     #du/dy(x,0,z,t)=0  Cara lateral 2
    u_k[:,-1,:] = 20 + 10*sin(2*pi/(24.0/(60./60./60.))*t)     #du/dy(x,0,z,t)=0  Cara de arriba
    print 20 + 10*sin(2*pi/(24.0/(60./60./60.))*t)
    #Loop en el espacio   i = 1 ... Nx-1, j = 1 ... Ny-1, Z = 1 ... Nz-1
    for i in range(1,Nx):
        for j in range(1,Ny):
            for Z in range(1,Nz):
                #Algoritmo de diferencias finitas 3-D para difusion
 
                #Laplaciano
                nabla_u_k = (u_k[i-1,j,Z] + u_k[i+1,j,Z] + u_k[i,j-1,Z] + u_k[i,j+1,Z] + u_k[i,j,Z-1] + u_k[i,j, Z+1] - 6*u_k[i,j,Z])/h**2

                #Forward euler..
                u_km1[i,j,Z] = u_k[i,j,Z] + alpha*nabla_u_k
 
    #CB natural
    u_km1[Nx,:,:] = u_km1[Nx-1,:,:]
    u_km1[:,Ny,:] = u_km1[:,Ny-1,:]   
    u_km1[:,:,Nz] = u_km1[:,:,Nz-1]
    #Avanzar la solucion a k + 1
    u_k = u_km1
    print u_k
    #CB esencial una ultima vez
    u_k[0,:,:] = u_k[1,:,:]             #du/dx(0,y,z,t)=0    Cara lateral 1
    u_k[:,0,:] = u_k[:,1,:]             #du/dy(x,0,,z,t)=0   Cara abajo
    u_k[:,:,0] = u_k[:,:,1]             #du/dx(x,y,0,t)=0    Caras del ancho 1
    u_k[:,:,Nz] = u_k[:,:,Nz-1]           #du/dy(x,y,,C,t)=0   Cara del ancho 2
    u_k[-1,:,:] = 20 + 10*sin(2*pi/(24.0/(60./60./60.))*t)     #du/dy(x,0,z,t)=0  Cara lateral 2
    u_k[:,-1,:] = 20 + 10*sin(2*pi/(24.0/(60./60./60.))*t)     #du/dy(x,0,z,t)=0  Cara de arriba
    
    print "Tmax = ", u_k.max()
    print ("t = {} Dias, {} Horas, {} Minutos".format(Dias, Hora, k%60))
    #time.append(k)
    #m.append(u_k[9,9,9])
#print time
#print m

#plot(time, m)
#show()
    #if t > next_t and i==xM:
    #    figure(1)
    #    imshowbien(u_k)
    #    title("k = {0:4.0f}   t = {1:05.2f} s".format(k, k*dt))
    #    savefig("E6_1/frame_{0:04.0f}.png".format(framenum))
    #    framenum += 1
    #    next_t += dnext_t
    #    close(1)
 
# figure(2)
# imshowbien(u_k)
# title("k = {}   t = {} s".format(k, (k+1)*dt))
 
#show()