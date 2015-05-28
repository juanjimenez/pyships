# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 20:10:35 2015
Calcula los coeficientes de un polinomio de Bézier a partir de los puntos
de control... 
Tambien de la opción de pintarlos para regocijo y alegría del usuario...
@author: juan
"""

import numpy as np
from scipy.misc import comb
import matplotlib.pyplot as pl

def bezier4p(p1,p2,v1,v2,tf,theta1,theta2,n):
    '''Esta es una función ad-hoc para calcularse puntos de una curva de
    Becier de grado 3 (cuatro puntos de paso que pasa por los puntos 
    p1 =[[x],[y]]y p2), Ojo a la estructura estan metidos como vectores columnas
    t es el valor del parámetro entre 0 y 1 que define en que punto de la
    curva estoy... t debe ser un valor real o un array de puntos'''
    
    t = np.linspace(0,1,num = n)*tf
        
    #definir los puntos de control 
    p = np.array([p1, p1 + np.array([[v1*tf/3*np.cos(theta1)],[v1*tf/3*np.sin(theta1)]]),\
    p2 - np.array([[v2*tf/3*np.cos(theta2)],[v2*tf/3*np.sin(theta2)]]),p2])
    
    #obtencion puntos de control de la odografa
    d = (p[1:]-p[:-1])        
    
    r = np.zeros((2,n))
    v = np.zeros((2,n))    
            
    for i in range(4):    
        r += tf**-3 * comb(3,i) * t ** i * (tf - t) ** (3-i) * p[i]
    
    
    for i in range(3):
        v += tf**-3 * comb(2,i) * t ** i * (tf - t) ** (2-i) * 3*d[i]
        
    return r,v,t

def pintar_bezier(r,v = None):
    
    pl.plot(r[0],r[1])
    pl.hold(True)
    if v is not None:
        for i,j,k,l in zip(r[0],r[1],v[0],v[1]):
            pl.arrow(i,j,k,l)
        