# Proyecto-1-MCOC
Proyecto 1 MCOC

# Introducción
Utilizando simulación numérica en diferencias finitas, se espera predecir la evolución térmica en hormigones masivos durannte su proceso de maduración.
Para etso se utiliza una ecuación diferencial en derivadas parciales conocida como la ecuación de difusión, la cual describe la evolucion de la temperatura por expulsión de calor en un sólido.

# [Entrega 3]- Simulación caso 1-D (10 casos)
Para la realización de esta entrega, se ejecuto un programa en Phython, el cual entrega de manera gráfica el vector u verus el tiempo, este vector corresponde a la temperatura en el tiempo.
Se consideraron 10 diferentes materiales, a los cuales se les aplico la siguiente ecuación:
```
 u_km1[i] = q*dt + u_k[i] + Alpha[j]*(u_k[i+1]-2*u_k[i]+u_k[i-1])
```
En donde q, se determino con la función:
```
def q_value(x,i):
    if (x[i]<L/3) or (x[i]>2*L/3):
        return 0
    else:
        return (t/3600)*x[i]
```
Y alpha se determino como:
```
Alpha[i]=K[i]*dt/(c[i]*rho[i]*dx**2)
```
En donde K, c y rho son parámetros dependientes de cada material.
El programa realiza un ciclo for en donde se va evaluando el cambio de la temeperatura de cada material por difusion a traves de cierto perido de tiempo k.

# [Entrega 4]- Simulación caso 2-D (3 casos)
Para esta entrega, se creó un código en Puthon intentando simular los resultados obtenidos en tres casos realizador por el profesor. En donde el primero contaba con dos condiciones de borde en donde la temperatura se mantenia constante en 20°C (eje x e y); el segundo caso consistia en que se añadia un borde mas en temperatura constenate, dejando el borde superior con temperatura ambiente la cual se debia  modelar con la siguiente función:
```
ambiente(t,T) = 20+10*np.sin((2*np.pi/T)*t
u_k[:,-1] = 20 + 10*sin(2*pi*t/(7))

```
Finalmente, el último caso, consistio en que tres de los bordes eran aisalntes, por lo que la pendiente entre estos puntos debiese ser 0. Se definio como du/dx = 0, du/dy = 0.

```
  u_k[0,:] = u_k[1,:]             #du/dx(0,y,t)=0
  u_k[Nx,:] = u_k[Nx-1,:]         #du/dy(b,y,t)=0
  u_k[:,0] = u_k[:,1]  

```

# [Entrega 6]- Simulación caso 3-D
Para esta última entrega, se simula la gráfica obtenida por los sensores en el experimento realizado en el Laboratorio Uandes; estos sensores fueron colocados en un bloque de hormigon de dimensiones .....
Como fue investigado para la [Entrega 5] lo parametros que se utilizaron fueron:
- K = conductividad térmica de 9.87[kJ/m h °C]
- c =  calor específico de 0.983 [kJ/ kg °C]
- rho = densidad de 2.425 [kg/m3]
Y con un calor máximo, entregado por la hidratación del cemento de HT=0,597 [kJ/g].
