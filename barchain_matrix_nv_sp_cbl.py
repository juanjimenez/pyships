# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 22:43:48 2013
Este fichero contiene las matrices de tensiones de los eslabones de la cadena 
calculadas para distintas situaciones cadenas, barcos en las puntas extremos,
libres,
DEPRECATED. See boom.py. strain matrices are at present included in cadena
There is a reason to preserve these files: barchain_matrix_x they include code
to include a three (central) ship towing the boom 
@author: juan
"""
import numpy as np
import numpy.linalg.linalg as npl
import boom_nv
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

#creamos una matriz pero esto depende de la situación
def Matriz_t_ini(cad = boom_nv.cadena()):
    ''' inicio de la matriz de tensiones de la cadena, si barciint es distinto
    de cero, se está inicializando la matriz suponiendo que hay un barco en el
    centro, con lo cual hay que calcular una tensión más. El resto de
    situaciones ¡¡¡dan igual???'''     
    

    #construimos una matriz sparse (csr) vacia de dimension el 2*numero de
    #eslabomes mas dos. Usamos csr (compressed sparse row format) porque es
    #el que mejor se adecua a la estructura de la matriz de coeficientes
    #de los formatos con lo que trabaja spsolve...
   
    M = csr_matrix((2 * cad.esl +2,2 * cad.esl +2)) 
    # Costruimos el puntero de elementos no nulos de la matriz
    M.indptr[1:3] = [4,8]
    M.indptr[3:-2] = np.arange(14,14+6*(M.indptr.size-5),6)
    M.indptr[-2:] = [M.indptr[-3] + 4, M.indptr[-3] + 8]      
    
    #cosntruimos los indices de columna de los elementos no nulos de la 
    #matriz
    M.indices= np.array([0,1,2,3,0,1,2,3,0,1,2,3,4,5,0,1,2,3,4,5],\
    dtype = 'int32')
    for i in range(4,M.indptr.size-4,2):
        M.indices= np.concatenate((M.indices,i-2 + M.indices[8:20]))
    M.indices= np.concatenate((M.indices, i + M.indices[0:8]))
    M.data = np.zeros(M.indices.size)
    B = np.zeros(2 * cad.esl +2)
    T = np.zeros(2 * cad.esl +2)
       
    return M, B, T 
        
        
def Matriz_t(M, B, T, cad, barcoi = None, barcod = None):
    ''' Calculo de la matriz de coeficientes para calcular las tensiones
    soportadas por los eslabones de la cadena. los extremos de la cadena
    pueden estan libres por defecto, además se puede enganchar un objeto
    barcoi al extremo izquierdo de la cadena, un objeto barcod al extremo
    derecho de la cadena y un objeto barcoc a un eslabon intermedio
    cualquiera. En principio, el eslabon concreto al que se engancha barcoc
    viene determinado por el atributo barcoc.link.'''
########################Revisado###############################################    
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
    
  
    for i in range(0,cad.normal.shape[1]-1): #modificamos las filas de
    #la matriz m 
        M.data[M.indptr[2*i+2]:M.indptr[2*i+4]] = \
        [T1xa[0,i],T1ya[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        T1ya[i+1],T1ya[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        T1ya[i+1], T1xa[1,i+1]]

    #calculo de los terminos independientes B
    #factores para construir los terminos independientes
    B[range(2,B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
    + factor2[0,1:] +factor2[0,0:-1]\
    
    B[range(3,B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
    + factor2[1,1:] +factor2[1,0:-1]
        
    
    #añadimos la contribución de los barcos si existen o el 'cierre' de las
    # matrices de coeficientes si no existen...(no me gusta la condicional)
    #pero de momento hasta que me aclare es lo que hay...
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
        if cad.Fi[0] == -1:
            
        #metemos esta condicion cuano el extremo izquierdo
        #esta fijo al muelle. En este caso, podemos eliminar las dos primeras 
        #ecuaciones de la matriz, ya que no hay que calcular tensiones a la
        #izquierda de la izquierda del primer eslabon... 
        #compio aqui de momento las condiciones impuestas a la matriz full que 
        #se emplea  en barchain_matrix_nv... esto esta en construccion hay que
        #ver como manejarlo...
            M[2,2:4] = M[2,2:4] + M[2,0:2]
            M[3,2:4] = M[3,2:4] + M[3,0:2]
            M = M[2:,2:]
            B = B[2:]
            T = T[2:]
        #recalculamos los coeficientes de la primera columna de A anadiendo
        #la condicion de que las tensiones a izquierda y derecha del primer
        #eslabon deben ser iguales...

        else:
            #este supuesto es que se aplica una fuerza externa directamente al
            #extremo de la cadena... posiblemente no es buena idea tener tanto
            #if then en un codigo pensado para iterar...
            M[0,0] = 1
            M[1,1] = 1
            B[0] = cad.Fi[0]
            B[1] = cad.Fi[1]
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
    
      
           
    T = spsolve(M,B)
#    print('dentro', T)        
    return M,B,T       