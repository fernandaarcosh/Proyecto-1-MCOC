from matplotlib.pylab import *

#Scipy.interpolate
#Caso en donde u(0,y,t)=u(x,0,t)=20 y du/dx(b,y,t)=du/dy(x,a,t)=0
  
a = 1.          #Alto del dominio y
b = 1.          #Largo del dominio x
C = 1.          #Ancho del dominio z
Nx = 5         #Numero de intervalos en x
Ny = 5         #Numero de intervalos en Y
Nz = 5         #Numero de intervalos en Z
 
dx = b / Nx     #Discretizacion espacial en X
dy = a / Ny     #Discretizacion espacial en Y
dz = C / Nz
 
h = dx    # = dy = dz 

xM = round((Nx/2)/10.)*10
if dx != dy:
    print("ERRROR!!!!! dx != dy")
    exit(-1)   #-1 le dice al SO que el programa fallo.....
 
#Funcion de conveniencia para calcular coordenadas del punto (i,j)
 
coords = lambda i, j, Z : (dx*i, dy*j, dz*Z)
x, y, z = coords(4,2,3) 
 
print "x = ", x
print "y = ", y
print "z = ", z
 
u_k = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)
u_km1 = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)
 
#Buena idea definir funciones que hagan el codigo expresivo
#def printbien(u):
 #   print u.T[Nx::-1,:]
 
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
 
k = 0
u_k[:,:,:]=20. 
#Loop en el tiempo 
dnext_t = 60   #  20.00
next_t = 0.
framenum = 0
for k in range(int32(5.*24*60+19*60.)):
    t = dt*(k+1)
    print "k = ", k, "t = ", t
 
    #CB esencial
    u_k[0,:,:] = u_k[1,:,:]             #du/dx(0,y,z,t)=0    Cara lateral 1
    u_k[:,0,:] = u_k[:,1,:]             #du/dy(x,0,,z,t)=0   Cara abajo
    u_k[:,:,0] = u_k[:,:,1]             #du/dx(x,y,0,t)=0    Caras del ancho 1
    u_k[:,:,-1] = u_k[:,:,-2]           #du/dy(x,y,,C,t)=0   Cara del ancho 2
    u_k[-1,:,:] = 20 + 10*sin(2*pi/(24.0/(60./60./60.))*t)     #du/dy(x,0,z,t)=0  Cara lateral 2
    u_k[:,-1,:] = 20 + 10*sin(2*pi/(24.0/(60./60./60.))*t)     #du/dy(x,0,z,t)=0  Cara de arriba
    
    #Loop en el espacio   i = 1 ... Nx-1, j = 1 ... Ny-1, Z = 1 ... Nz-1
    for i in range(1,Nx):
        if t > next_t and i==xM:
           figure(1)
           imshowbien(u_k)
           title("k = {0:4.0f}   t = {1:05.2f} s".format(k, k*dt))
           savefig("E6_1/frame_{0:04.0f}.png".format(framenum))
           framenum += 1
           next_t += dnext_t
           close(1)
        for j in range(1,Ny):
            for Z in range(1,Nz):
                #Algoritmo de diferencias finitas 3-D para difusion
 
                #Laplaciano
                nabla_u_k = (u_k[i-1,j,Z] + u_k[i+1,j,Z] + u_k[i,j-1,Z] + u_k[i,j+1,Z] + u_k[i,j,Z-1] + u_k[i,j, Z+1] - 6*u_k[i,j,Z])/h**2
 
                #Forward euler..
                u_km1[i,j,z] = u_k[i,j,z] + alpha*nabla_u_k
 
    #CB natural
    u_km1[Nx,:,:] = u_km1[Nx-1,:,:]
    u_km1[:,Ny,:] = u_km1[:,Ny-1,:]
    u_km1[:,:,Nz] = u_km1[:,:,Nz-1]
    #Avanzar la solucion a k + 1
    u_k = u_km1
 
    #CB esencial una ultima vez
    u_k[0,:,:] = u_k[1,:,:]             #du/dx(0,y,z,t)=0    Cara lateral 1
    u_k[:,0,:] = u_k[:,1,:]             #du/dy(x,0,,z,t)=0   Cara abajo
    u_k[:,:,0] = u_k[:,:,1]             #du/dx(x,y,0,t)=0    Caras del ancho 1
    u_k[:,:,-1] = u_k[:,:,-2]           #du/dy(x,y,,C,t)=0   Cara del ancho 2
    u_k[-1,:,:] = 20 + 10*sin(2*pi/(24.0/(60./60./60.))*t)     #du/dy(x,0,z,t)=0  Cara lateral 2
    u_k[:,-1,:] = 20 + 10*sin(2*pi/(24.0/(60./60./60.))*t)     #du/dy(x,0,z,t)=0  Cara de arriba
    
    print "Tmax = ", u_k.max()
    if t > next_t and i==xM:
        figure(1)
        imshowbien(u_k)
        title("k = {0:4.0f}   t = {1:05.2f} s".format(k, k*dt))
        savefig("E6_1/frame_{0:04.0f}.png".format(framenum))
        framenum += 1
        next_t += dnext_t
        close(1)
 
# figure(2)
# imshowbien(u_k)
# title("k = {}   t = {} s".format(k, (k+1)*dt))
 
 
show()