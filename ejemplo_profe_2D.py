from matplotlib.pylab import *

A=1        #Largo dominio en Y
B=1.       #Largo dominio en X
Nx=6       #n° intervalos
Ny=6
dx= B / Nx
dy= A / Ny
#print dx
#vector x e y
#x=linspace(0, A, Nx)
#y=linspace(0, B, Ny)

h=dx
def coords(i,j):        #funcion de conveniencia que sirve para calcular coordenadas del punto (i,j)
    return dx*i, dy*j

    
#x, y = coords(4, 2)
#i, j = 4, 2
#x, y =dx*i, dy*j

coords = lambda i, j : (dx*i, dy*j)
x, y =coords(4, 2)

print "x = ", x
print "y = ", y
u_k=zeros((Nx+1, Ny+1), dtype=double)   #dtype es el tipo de datos (double, float, int32)
u_km1=zeros((Nx+1, Ny+1), dtype=double)

u_k[0,:]=20
u_k[:,0]=20

def printbien(u):
    print u_k.T[Nx::-1,:]
def imshowbien(u):
    imshow(u.T[Nx::-1,:])

printbien(u_k)


#Valores usando Hierro
dt=1.
K=79.5    #W/(m*K)
c=450.    #J/(Kg*K)
rho=7800. #Kg/m^3

k=0
Alpha=K*dt/(c*rho*dx**2)
figure(1) 
imshowbien(u_k)
colorbar()
title("k = {}  t = {} s".format(k, k*dt))

for k in range(1000):
    t=dt*(k+1)
    #Definicion de las condiciones de bordes
    u_k[0,:]=20
    u_k[:,0]=20
    #Loop en el espacio
    for i in range(1,Nx-1):
        for j in range(1,Ny-1):
            nabla_u_k = (u_k[i-1,j] + u_k[i+1,j] + u_k[i,j+1] + u_k[i,j-1] - 4*u_k[i,j])/h**2
            u_km1[i,j]= u_k[i,j] + Alpha*nabla_u_k
    u_km1[Nx,:] = u_km1[Nx-1,:]
    u_km1[:,Ny] = u_km1[:,Ny-1]
    #Avanza al tiempo k+1
    u_k = u_km1
u_k[0,:]=20
u_k[:,0]=20
figure(2)         
imshowbien(u_k)
title("k = {}  t = {} s".format(k, k*dt))
colorbar()