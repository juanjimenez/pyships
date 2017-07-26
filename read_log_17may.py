# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:56:35 2017
Sirve para leer los log que producen los barquitos azules de Josemaria...
La parte de manipulacion que sigue a la lectura deberia estar en un fichero 
aparte, pero desgraciadamente somos asi de guarros.
@author: juan
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from scipy import stats
import dubing
import usv
###############################################################################
"""
Funcion para pasar de longitud y latitud a metros en xy esta copiada de
/home/juan/gldrive/pyships/data_exprmt/matread.py
"""
def traslada(datos,zero,d = 0):
    '''Trasforma datos en logitud latitud a x y, referidos al valor de cero 
    suministrado en zero, si d ==0, y de x y a longitud latitud
    se supones que datos se introuduce como una matriz cuya fila cero contiene
    las longitudes la fila uno las latitudes. devolvera en x las distancias
    correspondientes a las longitudes y en y las correspondientes a las
    latitudes. si se meten coordenadas en metros [x,y] de en datos entonces
    x se trasforma la longitudes e y a latitudes... (en lo barcos gordos se 
    hace justo al reves.)
    si d <> zero. 
    (el metodo me lo paso Martin y esta implementado tal cual...)'''
    deg2m = 111194.9266445587
    if d == 0:
        y = (datos[1] - zero[1]) * deg2m
        x = (datos[0] - zero[0]) * deg2m *np.cos(datos[1]*np.pi/180.)
        dev = [x,y]
    else: 
        lat = datos[1]/deg2m + zero[1]    
        lon = datos[0]/(deg2m*np.cos(lat*np.pi/180.))+zero[0]
        dev = [lat,lon]
    return dev
###############################################################################    
#Cargamos los datos ##########################################################
try:
# GPS =\
# pd.read_table('E:\JUANFRAN JIMENEZ\Google Drive\\pyships\Visualiz\\GPS_10.log',\
# sep = '\s+', encoding='latin_1',na_values = " ")   
# Actuators =\
# pd.read_table('E:\JUANFRAN JIMENEZ\Google Drive\\pyships\Visualiz\\Actuators_10.log',\
# sep = '\s+',encoding='latin_1',na_values = " ")
# Compass =\
# pd.read_table('E:\JUANFRAN JIMENEZ\Google Drive\\pyships\Visualiz\\Compass_10.log',\
# sep = '\s+',encoding='latin_1',na_values = " ")
# ZodiacState =\
# pd.read_table('E:\JUANFRAN JIMENEZ\Google Drive\\pyships\Visualiz\\ZodiacState_10.log',\
# sep = '\s+',encoding='latin_1',na_values = " ")
# Formation =\
# pd.read_table('E:\JUANFRAN JIMENEZ\Google Drive\\pyships\Visualiz\\Formation_10.log',\
# sep = '\s+',encoding='latin_1',na_values = " ")
 GPS =\
 pd.read_table('C:\Users\juan\Google Drive\\pyships\Visualiz\\GPS_11.log',sep = '\s',\
 encoding='latin_1',na_values = " ")   
 Actuators =\
 pd.read_table('C:\Users\juan\Google Drive\\pyships\Visualiz\\Actuators_11.log',sep = '\s+',\
 encoding='latin_1',na_values = " ")
 Compass =\
 pd.read_table('C:\Users\juan\Google Drive\\pyships\Visualiz\\Compass_11.log',sep = '\s+',\
 encoding='latin_1',na_values = " ")
 ZodiacState =\
 pd.read_table('C:\Users\juan\Google Drive\\pyships\Visualiz\\ZodiacState_11.log',sep = '\s+',\
 encoding='latin_1',na_values = " ")
 Formation =\
 pd.read_table('C:\Users\juan\Google Drive\\pyships\Visualiz\\Formation_11.log',sep = '\s+',\
 encoding='latin_1',na_values = " ")
except:
 GPS=\
 pd.read_table('/home/juan/gldrive/pyships/Visualiz/GPS_11.log',sep = '\s+',\
 encoding='latin_1',na_values = " ")   
 Actuators =\
 pd.read_table('/home/juan/gldrive/pyships/Visualiz/Actuators_11.log',sep = '\s+',\
 encoding='latin_1',na_values = " ")
 Compass =\
 pd.read_table('/home/juan/gldrive/pyships/Visualiz/Compass_11.log',sep = '\s+',\
 encoding='latin_1',na_values = " ")
 ZodiacState =\
 pd.read_table('/home/juan/gldrive/pyships/Visualiz/ZodiacState_11.log',sep = '\s+',\
 encoding='latin_1',na_values = " ")
 Formation =\
 pd.read_table('/home/juan/gldrive/pyships/Visualiz/Formation_11.log',sep = '\s+',\
 encoding='latin_1',na_values = " ")

#adding a column with real times. Mola pero es inutil porque los datos se 
#repiten segun los muestreo de cada variable que ni siquiera coinciden...
#Es mas útil el promedio bestia hecho despues...
GPS.loc[:,'time_r']=GPS.time_sec +GPS.time_usec * 10**-6 -GPS.time_sec[0]
Actuators.loc[:,'time_r']=Actuators.time_sec +Actuators.time_usec * 10**-6-\
Actuators.time_sec[0]
Compass.loc[:,'time_r']=Compass.time_sec +Compass.time_usec * 10**-6-\
Compass.time_sec[0]
ZodiacState.loc[:,'time_r']=ZodiacState.time_sec-\
ZodiacState.time_usec * 10**-6 -ZodiacState.time_sec[0]
Formation.loc[:,'time_r']=Formation.time_sec +Formation.time_usec * 10**-6-\
Formation.time_sec[0]
##############################################################################
#hacemos un promedio a lo bruto de los valores por segundos....

GPStemp = GPS.groupby(['time_sec']).mean()
GPStemp = GPStemp.reset_index()
Compasstemp = Compass.groupby(['time_sec']).mean()
Compasstemp = Compasstemp.reset_index()

#tomamos solo los puntos de la traza GPS que nos intesan
ini = 652
fin = 2700
lonlat = GPStemp[['longitude','latitude']]\
[(GPStemp.time_sec>=GPS.time_sec[ini])&(GPStemp.time_sec<=GPS.time_sec[fin])].get_values()
compassxy = Compasstemp[['compass_x','compass_y']]\
[(Compasstemp.time_sec>=GPS.time_sec[ini])&(Compasstemp.time_sec<=GPS.time_sec[fin])].get_values()
 
rlonlat = Formation[['referenceLon','referenceLat']][ini:fin].get_values()
#los pasamos a metros, poniendo el cero en el primer puntos. ojo ya veremos

datos_m = traslada(lonlat.T,lonlat[0])
datosr_m  = traslada(rlonlat.T,lonlat[0])
#pl.plot(datos_m[0],datos_m[1],'k')
#pl.plot(datosr_m[0],datosr_m[1],'k')


pl.figure(1)
bexp = usv.barco(1)
bexp.ls = 1.
color = 'r'
heading = np.arctan2(compassxy.T[1],compassxy.T[0])
flip = 0 
for i in zip(heading,np.array(datos_m).T):
    bexp.pb = i[1]
    if flip > 184:
        color = 'r'
    bexp.theta = i[0]
    bexp.dibujar(color)
    flip += 1
pl.xlabel('m') 
pl.ylabel('m')   
axis('equal')

#dubings pa completar la fiesta
wp = [np.array([0.,0.,0,-1]),\
np.array([-48.5,69,14.,-1]),\
np.array([-92.5,35.5,14,-1]),\
np.array([-48.5,69,14.,-1])]

###dubin secondary waypoints:
wp2 = dubing.secnwayp(wp)
wp2[-1].append(wp2[1][1])
dubing.pintadubing(wp2)      
#vel = GPStemp['velocity']\
#[(GPStemp.time_sec>=GPS.time_sec[2585])&(GPStemp.time_sec<=GPS.time_sec[2906])].get_values()

