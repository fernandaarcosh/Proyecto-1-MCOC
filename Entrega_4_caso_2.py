from matplotlib.pylab import *

#Caso en donde u(0,y,t)=u(x,0,t)=u(b,y,t)=20 y u(x,a,t)=20 + 10*sin(2*pi*t/(7))

a = 1.          #Ancho del dominio y
b = 1.          #Largo del dominio x
Nx = 70         #Numero de intervalos en x, si se intenta usar el rango que el profesor desea de 100 el codigo no funciona de manera optima
Ny = 70         #Numero de intervalos en Y, 70 es el mayor numero donde el codigo entrega resultados estables.
 
dx = b / Nx     #Discretizacion espacial en X
dy = a / Ny     #Discretizacion espacial en Y
 
h = dx    # = dy
 
 
if dx != dy:
    print("ERRROR!!!!! dx != dy")
    exit(-1)   #-1 le dice al SO que el programa fallo.....
 
#Funcion de conveniencia para calcular coordenadas del punto (i,j)
coords = lambda i, j : (dx*i, dy*j)
x, y = coords(4,2) 
 
print "x = ", x
print "y = ", y
 
u_k = zeros((Nx+1,Ny+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)
u_km1 = zeros((Nx+1,Ny+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)
 
#CB esencial
u_k[:,:] = 20.    #Para lograr situacion inicial del material donde esta a la misma temperatura en todas partes
#u_k[:,Ny] = 20 + 10*sin(2*pi*t/(24*3600))
 
#Buena idea definir funciones que hagan el codigo expresivo
def printbien(u):
    print u.T[Nx::-1,:]
 
#print u_k               
printbien(u_k)                  #Imprime con el eje y invertido
 
def imshowbien(u):
    imshow(u.T[Nx::-1,:])
    colorbar(extend='both',cmap='plasma')
    clim(10, 30)
 
#Parametros del problema (hierro)
dt = 60.       # s, para que avance cada 1 minuto
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
 
#Loop en el tiempo 
dnext_t = 0.05   #  20.00
next_t = 0.
framenum = 0
for k in range(int32(5./dt)):
    t = dt*(k+1)
    print "k = ", k, " t = ", t
 
    #CB esencial
    u_k[0,:] = 20.
    u_k[:,0] = 20.
    u_k[-1,:] = 20.
    u_k[:,-1] = 20 + 10*sin(2*pi*t/(7))
 
    #Loop en el espacio   i = 1 ... n-1   u_km1[0] = 0  u_km1[n] = 20
    for i in range(1,Nx):
        for j in range(1,Ny):
            #Algoritmo de diferencias finitas 2-D para difusion
            
 
            #Laplaciano
            nabla_u_k = (u_k[i-1,j] + u_k[i+1,j] + u_k[i,j-1] + u_k[i,j+1] - 4*u_k[i,j])/h**2
 
            #Forward euler..
            u_km1[i,j] = u_k[i,j] + alpha*nabla_u_k
 
    #CB natural
    u_km1[Nx,:] = u_km1[Nx-1,:]   #pendiente 
    u_km1[:,Ny] = u_km1[:,Ny-1]
 
    #Avanzar la solucion a k + 1
    u_k = u_km1
 
    #CB esencial una ultima vez
    u_k[0,:] = 20.
    u_k[:,0] = 20.
    u_k[-1,:] = 20.
    u_k[0,:] = u_k[1,:]
    u_k[Nx,:] = u_k[Nx-1,:]
    u_k[:,0] = u_k[:,1] 
    u_k[:,Ny] = 20+10*np.sin((2*np.pi/7)*t)
 
    print "Tmax = ", u_k.max()
 
    if t > next_t:
        figure(1)
        imshowbien(u_k)
        title("k = {0:4.0f}   t = {1:05.2f} s".format(k, k*dt))
        savefig("E4_caso_2/frame_{0:04.0f}.png".format(framenum))
        framenum += 1
        next_t += dnext_t
        close(1)
 
show()