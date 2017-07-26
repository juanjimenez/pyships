# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 22:43:48 2013
Este fichero define un objeto llamado cadena. El objeto contiene una cadena de
logitud variable, y tiene especificadas una serie de propiedades mecanicas de 
las cadenas.
Ademas, incluye funciones para pintar cadenas, moverlas etc.
Nota: se distingue de los ficheros  barchain_matrix_cbl_sp, en que estos ultimos
utilizan matrices sparse
Esta version incluye masas a~nadidas, resitencias cuadraticas al avance y los
extremos de la cadena estan unidos a los barcos mediante cables.
Se ha excluido el caso de arratre de tres barcos

DEPRECATED. See boom.py. strain matrices are at present included in cadena
There is a reason to preserve these files: barchain_matrix_x they include code
to include a three (central) ship towing the boom 
objects. 06.03.2017
 
@author: juan
"""
import numpy as np
import numpy.linalg.linalg as npl
import boom_nv

#creamos una matriz pero esto depende de la situación
def Matriz_t_ini(cad = boom_nv.cadena()):
    ''' inicio de la matriz de tensiones de la cadena hay dos posibilidades
    1. Hay dos barcos en los extemos unidos al boom por cables
    2. Solo hay barco a la derecha de la cadena el otro extremo esta fijado a 
       Tierra
    '''     
    M = np.zeros((2 * cad.esl + 2, 2 * cad.esl + 2))
    B = np.zeros(2 * cad.esl + 2)
    T = np.zeros(2 * cad.esl + 2)
    return M, B, T 
        
        
def Matriz_t(M, B, T, cad, barcoi = None, barcod = None):
    ''' Calculo de la matriz de coeficientes para calcular las tensiones
    soportadas por los eslabones de la cadena. los extremos de la cadena
    pueden estan libres por defecto, además se puede enganchar un objeto
    barcoi al extremo izquierdo de la cadena, un objeto barcod al extremo
    derecho de la cadena y un objeto barcoc a un eslabon intermedio
    cualquiera. En principio, el eslabon concreto al que se engancha barcoc
    viene determinado por el atributo barcoc.link.'''
    #print 'en cadena', barcoi.Fm
    #print 'en cadena', barcod.Fm
#######################Revisado################################################
#calculamos factores que se repiten frecuentemente en el calculo de los
#elementos de la matriz
    mcA = 1./(cad.m + cad.mA * cad.normal) #inversas masas a~ndidas eslabones
    T1xa = mcA -cad.L**2 * cad.normal**2 / cad.I
    T1ya = - cad.L**2 * cad.normal[0] * cad .normal[1] / cad.I
    
    #resistencia al avance Fr
    factor1 = mcA *\
    ((cad.s * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.vr * cad.para,0))) / cad.vmodr+\
    (cad.s2 * np.abs(np.sum(cad.vr * cad.normal,0)) +\
    cad.q2 * abs(np.sum(cad.vr * cad.para,0)))) * cad.vr
    #codiciones de cierre y resistencia al giro    
    factor2 =  cad.L * (cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal /cad.I) 
 

        
    T2xa = T1xa[:,:-1] -2 * mcA[:,:-1] + T1xa[:,1:] -2 * mcA[:,1:] 
    T2ya = T1ya[:-1] + T1ya[1:]  
        
      
    for i in range(0,cad.normal.shape[1]-1): 
    #modificamos las filas de la matriz m, correspondientes a eslabones
    #interiores, los terminos de los extrremos de la cadena se tratan a parte
        M[2*i+2:2*i+4,2*i:2*i+6] = np.array(\
        [[T1xa[0,i],T1ya[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        T1ya[i+1]],[T1ya[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        T1ya[i+1], T1xa[1,i+1]]])
    
        #calculo de los terminos independientes B
        #factores para construir los terminos independientes
        B[range(2,B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
        + factor2[0,1:] +factor2[0,0:-1]\
        
        B[range(3,B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
        + factor2[1,1:] +factor2[1,0:-1]
        
###############################################################################    
    #añadimos la contribución de los barcos si existen o el 'cierre' de las
    # matrices de coeficientes si no existen...
###############################################################################
    if barcoi:
        #calculamos la distancia entre primer eslabon y barco. Asumimos que 
        #siempre hay cable si hay barco
        vcbli =  barcoi.pb + barcoi.csbt([-barcoi.ls/2,0])\
        -cad.para[:,0] * cad.L - cad.cms[:,0] 
        disti = npl.norm(vcbli)
        dt = disti >= cad.cbl[0]
        
        #El cable esta tenso, si dt = 1
        #hacemos vcbl unitario que es tanto como guardar seno y coseno del cable
        vcbli = vcbli/(disti + (disti == 0.))  
        
        M[0,0] = 1
        M[1,1] = 1
        #metemos fuerza proporcional a la distancia        
        B[0] = -10000 * (disti -cad.cbl[0]) * dt * vcbli[0]
        B[1] = -10000 * (disti -cad.cbl[0]) * dt * vcbli[1]
                                         

        
        
            

    else:
        #que mas puede pasar si no hay barco a la izquierda
        if cad.Fi[0] == -1: #metemos esta condición cuando el extremo izquierdo
            #esta fijo al muelle... asi es que -1
            #quitamos las dos primeras ecuaciones, no hay que calcular tensiones
            #a la izquierda del primer eslabon porque es fijo...
            M[2,2:4] = M[2,2:4] + M[2,0:2]
            M[3,2:4] = M[3,2:4] + M[3,0:2]
            M = M[2:,2:]
            B = B[2:]
            T = T[2:]
            #recalculamos los coeficientes de la primera columna de A anadiendo
            #la condicion de que las tensiones a izquierda y derecha del primer
            #eslabon deben ser iguales...
                        
        else:
            #el extremo esta sujeto a una fuerza arbitraria cad.Fi
            M[0,0] = 1
            M[1,1] = 1
            B[0] = cad.Fi[0]
            B[1] = cad.Fi[1]
###############################################################################        
        
###############################################################################        
        
    if barcod:
       vcbld =  cad.cms[:,-1] - cad.para[:,-1] * cad.L \
       - barcod.csbt([-barcod.ls/2,0]) - barcod.pb 
       distd = npl.norm(vcbld)
       dt = distd >= cad.cbl[1]
       
       #El cable esta tenso, si dt = 1
       #hacemos vcbl unitario que es tanto como guardar seno y coseno del cable
       vcbld = vcbld/(distd + (distd == 0.))
       

       M[-1,-1] = 1
       M[-2,-2] = 1
       B[-2] = -10000 * (distd -cad.cbl[1]) * dt * vcbld[0]
       B[-1] = -10000 * (distd -cad.cbl[1]) * dt * vcbld[1]
    else:
       M[-1,-1] = 1
       M[-2,-2] = 1
       B[-2] = cad.Fd[0]
       B[-1] = cad.Fd[1]
       
      
    #print 'izq', disti, 'dcha', distd ,'\n'      
    T = np.linalg.linalg.solve(M,B)
#    print('dentro', T)
#    if barcoi:
#        T[0] = T[1] * vcbli[0]
#        T[1] = T[1] * vcbli[1]
#    if barcod:
#        T[-2] = T[-1] * vcbld[0]
#        T[-1] = T[-1] * vcbld[1]
#    print 'i', disti, 'd', distd
    return M,B,T       