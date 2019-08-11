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
        return t/x[i]
```
Y alpha se determino como:
```
Alpha[i]=K[i]*dt/(c[i]*rho[i]*dx**2)
```
En donde K,c y rho son parámetros dependientes de cada material.
El programa realiza un ciclo for en donde se va evaluando el cambio de la temeperatura de cada material por difusion a traves de cierto perido de tiempo k.
