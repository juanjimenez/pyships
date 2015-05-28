# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 20:51:26 2014
Curvitas de bezier de grado tres (4 puntos de control) Para interpolar un 
conjunto de puntos.
¿Sera esto lo que los expertos llaman B-splines?
@author: juan
"""
import numpy as np
import scipy as sc
import matplotlib as mpl

def Bzi(puntos,alpha,beta):
    """ Se espera que Bzi sea un array de puntos sobre los que interpolar
    Se espera que tenga tantas filas como puntos[[coordedas x,coordenadas y],...]
    La idea es cosntruir una curvita de Bezier entre cada dos puntos suministrados
    alpha marca la dirección de entrada. beta la direccion de salida. suponemos
    de momento v inicial = 1 v final = 1
    en en array"""

    #creamos un array para guardar los puntos de control. entre cada dos puntos 
    #suministrados debemos crear dos nuevos...
    
    p_control =np.zeros([3*(puntos.shape[0]-1)+1,2])
    p_control[0] = puntos[0]
    #el primer punto de control asegura la direccion de entrada del movil en
    #la trayectoria
    p_control[1] = np.array([np.cos(alpha),np.sin(alpha)])/3.0 + p_control[0]
    #el ultimo punto de control es el ultimo punto de paso
    p_control[-1] = puntos[-1]
    #fijamos el penultimo punto para asegurar la dirección de salida del movil
    # de la trayetoria
    p_control[-2] = p_control[-1]- np.array([np.cos(beta),np.sin(beta)])/3
    
    #de los 2 puntos de paso intermedio, el primero  lo fija la continuidad
    # de la trayectoria. Podemos fijar el segundo como el punto medio
    #entre el primero y el punto final del tramo
    for p in range(1,puntos.shape[0]-1):
        p_control[3*p] = puntos[p] #son puntos de paso
        
        #segundo punto
        p_control[3*p-1] = p_control[3*p-2]+(p_control[3*p] - p_control[3*p-2])/2
        #primer punto de la siguiente
        p_control[3*p+1] = 2 * p_control[3*p] -p_control[3*p-1]
        

    return p_control
    
def cypbz(p_control,paso,pinta):
    '''Esta función permite construir los polinomios y pintarlo. p_control
    es un array de puntos de control, paso define el numero de puntos de
    paso para el parametro t de los polinomios de bernstein y por ultimo
    pinta es un parametro que puede ser 0 o 1 para pintar o no el resultado
        
    Como todos los polinomios son de grado 3. los polinomios de bezier
    tienen una estructura fija. definimos los coeficientes directamente 
    en un array.'''
        
    B = np.array([1,3,3,1])
    #el numero de puntos generado nos lo da paso y el numero de puntos
    # paso divido por 4
        
    t = np.arange(0,1,paso)
    pasos = np.zeros([(p_control.shape[0]-1)*t.shape[0]/3-1,2])
        
    for i in range(0,p_control.shape[0]-3,3):
        #contruimos el polinomio de berstein para cada cuatro puntos
        for j in range(t.shape[0]):
            print i
            print j
            print i*(t.shape[0]-1)/3+j
            pasos[i*(t.shape[0]-1)/3+j] = B[0]*p_control[i]*(1-t[j])**3\
            + B[1]*p_control[i+1]*(1-t[j])**2*t[j]\
            + B[2]*p_control[i+2]*(1-t[j])*t[j]**2\
            + B[3]*p_control[i+3]*t[j]**3
            #print pasos
            
    #añadimos el último a mano
    pasos[-1] = p_control[-1]
            
    
    return pasos
        
        
        
        
    
   
    
