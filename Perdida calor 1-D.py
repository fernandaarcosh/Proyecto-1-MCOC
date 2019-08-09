from matplotlib.pylab import *

L=1.       #Largo dominio
n=100     #n° intervalos
dx= L / n
#print dx
#vector x
x=linspace(0, L, n+1)

#condicion inicial
def fun_u0(x):
    return 10*exp(-(x-0.5)**2/0.1**2)

u0 = fun_u0(x)

#U en el tiempo K
u_k = u0.copy()

#Condiciones de borde
u_k[0] = 0
u_k[n] = 20

#T° en tiempo k+1  t = dt*[k+1]
u_km1=u_k.copy()
dt=1.
K=79.5
c=450.
rho=7800.

Alpha=K*dt/(c*rho*dx**2) 
#print Alpha
#Loop en el tiempo
plot(x,u0, "k--")
k = 0
for k in range(1000):
    t=dt*k
    print "k=", k, "t=", t
    
    u_k[0]=0
    u_k[n]=20.
    #Loop en el espacio
    for i in range(1,n): 
        #print i        
        u_km1[i] = u_k[i] + Alpha*(u_k[i+1]-2*u_k[i]+u_k[i-1])
    u_k = u_km1
        
    if k % 50 == 0:
        plot(x,u_k)

title('k = {} t = {} s'.format(k, k*dt))

show()