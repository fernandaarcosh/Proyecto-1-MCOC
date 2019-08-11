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

#condicion de calor dependiendo de x y de t
def q_value(x,i):
    if (x[i]<L/3) or (x[i]>2*L/3):
        return 0
    else:
        return t/x[i]
    
u0 = fun_u0(x)

#U en el tiempo K
u_k = u0.copy()

#Condiciones de borde
u_k[0] = 0
u_k[n] = 20

#T° en tiempo k+1  t = dt*[k+1]
u_km1=u_k.copy()
dt=1.      #s

#Materiales a utilizar y sus respectivos valores
Materiales=["Hierro", "Madera", "Hormigon", "Acero dulce", "Poliestireno Expandido", "Pastico (Polipropile)", "Marmol", "Vidrio", "Aluminio", "Arena"]
K=[79.5, 0.13, 1.4, 49., 0.0373, 0.157, 2.09, 0.81, 209., 0.33] # W/(m*K)
c=[450., 1381., 837., 460., 1200, 1800., 879., 833., 909., 795.] #J/Kg*K
rho=[7800., 840., 2200., 7800., 25., 946., 2400., 2600., 2700., 1500.] #Kg/m^3

#Calcula los Alpha para los distintos materiales
Alpha = range(10)
for i in range(len(Alpha)):
    Alpha[i]=K[i]*dt/(c[i]*rho[i]*dx**2)

plot(x,u0, "k--")
#Loop en el tiempo para los distintos materiales
k = 0
for j in range(len(Alpha)):
    for k in range(1000):
        t=dt*k
        #print "k=", k, "t=", t
        u_k[0]=0
        u_k[n]=20
        #Loop en el espacio
        for i in range(1,n):
            q = q_value(x,i)
            u_km1[i] = q*dt + u_k[i] + Alpha[j]*(u_k[i+1]-2*u_k[i]+u_k[i-1])
        u_k = u_km1   
        if k % 50 == 0:    #Grafica cada 50 segundos
            plot(x,u_k)
    
        title('Material: {} ; t = {} s'.format(Materiales[j], k*dt))
        xlabel('Largo barra (m)')
        ylabel('Temperatura (C)')
    print 'Alpha de {} = {}'.format(Materiales[j], Alpha[j])
    #Mostrar el grafico de cada material en grafico separado
    show()
