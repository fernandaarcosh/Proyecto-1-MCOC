# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:51:20 2019

@author: Administrador
"""
from matplotlib.pylab import *
from scipy import interpolate
import numpy

archivo = open('TemperaturaAmbiente.csv')
datos=[]
n=0
for i in archivo:
    datos.append((i.split(';')[0]).split(','))#[Fecha,Hora,sen,temp]
    cambio = datos[-1][-1][:-1]
    #print cambio
    datos[-1][-1]=cambio
    #print datos[-1]
    #print datos[n]
    n+=1
    
ejeTiempo=[]
Tambiente=[]
for i in datos:    #Loop para asegurar que el grafico de la temperatura ambiente esta correcto
    ejeTiempo.append(float(i[2]))
    Tambiente.append(float(i[4]))
#pyplot.plot(ejeTiempo,grafico[0])

tempAmb=interpolate.interp1d(array(ejeTiempo),Tambiente,kind='slinear')  #Interpolacion de los datos
t=linspace(0,12,10000)  #Rango de los dias de datos otorgados
plot(t,tempAmb(t))
archivo = close('TemperaturaAmbiente.csv')