# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 22:43:48 2013
Este fichero define un objeto llamado cadena. El objeto contiene una cadena de
logitud variable, y tiene especificadas una serie de propiedades mecanicas de 
las cadenas.
Ademas, incluye funciones para pintar cadenas, moverlas etc.
Nota: se distingue de los ficheros  barchain_matrix_sp, en que estos ultimos
utilizan matrices sparse
Nota2: Esta nueva version (nv) incluye masas a~nadidas y resistencias cuadraticas

DEPRECATED. See boom.py. strain matrices are at present included in cadena
There is a reason to preserve these files: barchain_matrix_x they include code
to include a three (central) ship towing the boom 
@author: juan
"""
import numpy as np
import boom

#creamos una matriz pero esto depende de la situación
def Matriz_t_ini(cad = boom.cadena() , barcoint = 0):
    ''' inicio de la matriz de tensiones de la cadena Hay de momento dos
    posibilidades que una cubre el caso de que no haya barcos en los extremos'''     
    if barcoint == 0:    
        M = np.zeros((2 * cad.esl +2, 2 * cad.esl +2))
        B = np.zeros(2 * cad.esl +2)
        T = np.zeros(2 * cad.esl +2)
        #recalculamos aqui porque si ha cambiado la masa del eslabon
        #tenemos problemas
        #cad.mL2 = cad. m * cad.L**2
        #cad.I = 0.5 * cad.mL2 
    else:
        M = np.zeros((2 * cad.esl +4, 2 * cad.esl +4))
        B = np.zeros(2 * cad.esl +4)
        T = np.zeros(2 * cad.esl +4)
        #recalculamos aqui porque si ha cambiado la masa y/o la longitud del 
        #eslabon tenemos problemas
        #En cualquier caso no se si estos valores tiene ya mucho sentido
        #definirlos, la masa cambia dinamicamente al emplear masas a~nadidas
        #cad.mL2 = cad.m * cad.L**2
        #cad.I = cad.mL2/12
        #cad.I = 0.5 * cad.mL2
    return M, B, T 
        
        
def Matriz_t(M, B, T, cad, barcoi = None, barcod = None, barcoc = None):
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
 
    if barcoc:
       
        T2xa =np.concatenate((T1xa[:,:barcoc.link-2]-2*mcA[:,:barcoc.link-2]+\
        T1xa[:,1:barcoc.link-1] - 2 * mcA[:,1:barcoc.link-1],\
        [[0],[0]],[[0],[0]],\
        T1xa[:,barcoc.link-1:-1] -2*mcA[:,barcoc.link-1:-1] +\
        T1xa[:,barcoc.link:] -2*mcA[:,barcoc.link:]),axis = 1)
        T2ya = np.concatenate((T1ya[:barcoc.link-2]\
        + T1ya[1:barcoc.link-1],\
        [0],[0],\
        T1ya[barcoc.link-1:-1]\
        + T1ya[barcoc.link:]), axis = 1)
        
      
        for i in range(0,barcoc.link-2): #modificamos las filas de
        #la matriz m 
            M[2*i+2:2*i+4,2*i:2*i+6] = np.array(\
            [[T1xa[0,i],T1ya[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
            T1ya[i+1]],[T1ya[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
            T1ya[i+1], T1xa[1,i+1]]])
        for i in range(barcoc.link,cad.normal.shape[1]): #modificamos las filas de
        #la matriz m 
            M[2*i+2:2*i+4,2*i:2*i+6] = np.array(\
            [[T1xa[0,i-1],T1ya[i-1],T2xa[0,i],T2ya[i],T1xa[0,i],\
            T1ya[i]],[T1ya[i-1],T1xa[1,i-1],T2ya[i],T2xa[1,i],\
            T1ya[i],T1xa[1,i]]])
        #calculo de los terminos independientes B
        #factores para construir los terminos independientes
        
        B[range(2,2*barcoc.link-3,2)] =\
        + factor1[0,1:barcoc.link-1] - factor1[0,0:barcoc.link-2]\
        + factor2[0,1:barcoc.link-1] + factor2[0,0:barcoc.link-2]\
        
        B[range(3,2*barcoc.link-2,2)] =\
        + factor1[1,1:barcoc.link-1] - factor1[1,0:barcoc.link-2]\
        + factor2[1,1:barcoc.link-1] + factor2[1,0:barcoc.link-2]\
        
        B[range(2*barcoc.link+2,B.shape[0]-3,2)] =\
        + factor1[0,barcoc.link:] - factor1[0,barcoc.link-1:-1]\
        + factor2[0,barcoc.link:] + factor2[0,barcoc.link-1:-1]\
        
        B[range(2*barcoc.link+3,B.shape[0]-2,2)] = factor1[1,barcoc.link:]\
        -factor1[1,barcoc.link-1:-1]\
        + factor2[1,barcoc.link:] +factor2[1,barcoc.link-1:-1]
        
    else:
        
        T2xa = T1xa[:,:-1] -2 * mcA[:,:-1] + T1xa[:,1:] -2 * mcA[:,1:] 
        T2ya = T1ya[:-1] + T1ya[1:]  
        
      
        for i in range(0,cad.normal.shape[1]-1): #modificamos las filas de
        #la matriz m 
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
        
    
    #añadimos la contribución de los barcos si existen o el 'cierre' de las
    # matrices de coeficientes si no existen...(no me gusta la condicional)
    #pero de momento hasta que me aclare es lo que hay...
    if barcoi:
        #calculamos las masas a~nadidas y otros terminos que se repiten
          
        comps = np.array([np.cos(barcoi.theta),np.sin(barcoi.theta)])
        compt = np.array([comps[1],-comps[0]])
        mbA =np.array(\
        [barcoi.mb + barcoi.mA[0] * comps[0] - barcoi.mA[1] * comps[1],
        barcoi.mb + barcoi.mA[0] * comps[1] + barcoi.mA[1] * comps[0]])
        #rellenamos los elementos de la matriz de tensiones correspondientes 
        #al primer eslabon (ocupa la posición 0) y los terminos independientes       
        M[0,0:4] =\
        1 + mbA[0] * (2 *mcA[0,0] - T1xa[0,0]\
        + (barcoi.ls/2. * comps[1])**2/barcoi.Ib),\
        - mbA[0] * (T1ya[0]+ barcoi.ls**2/4. * comps[0] * comps[1]/barcoi.Ib),\
        - mbA[0] * T1xa[0,0],\
        - mbA[0] * T1ya[0]
        
        M[1,0:4] =\
        - mbA[1] * (T1ya[0]+ barcoi.ls**2/4. * comps[0] * comps[1]/barcoi.Ib),\
        1 + mbA[1] * (2 * mcA[1,0] - T1xa[1,0]\
        + (barcoi.ls/2. * comps[0])**2/barcoi.Ib),\
        - mbA[1] * T1ya[0],\
        - mbA[1] * T1xa[1,0],
        
        B[0] = -barcoi.Fb[0]\
        + barcoi.mul * np.dot(barcoi.vr,comps) * comps[0]\
        + barcoi.mut * barcoi.ls * np.dot(barcoi.vr,compt) * compt[0]\
        + barcoi.mul2 * np.sign(np.dot(barcoi.vr,comps))\
        * np.dot(barcoi.vr,comps)**2 * comps[0]\
        + barcoi.mut2 * barcoi.ls * np.sign(np.dot(barcoi.vr,compt))\
        * np.dot(barcoi.vr,compt)**2 * compt[0]\
        - mbA[0] * factor1[0,0]\
        - mbA[0] * cad.L * cad.w[0]\
        * (cad.para[0,0]  * cad.w[0] - cad.A * cad.normal[0,0]/cad.I)\
        - mbA[0] * barcoi.ls/2 * (barcoi.wb**2 * comps[0]\
        + (barcoi.M  - barcoi.ls * barcoi.mua * barcoi.wb) * comps[1]/barcoi.Ib) 
        
        B[1] = -barcoi.Fb[1]\
        + barcoi.mul * np.dot(barcoi.vr,comps) * comps[1]\
        + barcoi.mut * barcoi.ls * np.dot(barcoi.vr,compt) * compt[1]\
        + barcoi.mul2 * np.sign(np.dot(barcoi.vr,comps))\
        * np.dot(barcoi.vr,comps)**2 * comps[1]\
        + barcoi.mut2 * barcoi.ls * np.sign(np.dot(barcoi.vr,compt))\
        * np.dot(barcoi.vr,compt)**2 * compt[1]\
        - mbA[1] * factor1[1,0]\
        - mbA[1] * cad.L * cad.w[0]\
        * (cad.para[1,0]  * cad.w[0] - cad.A * cad.normal[1,0]/cad.I)\
        - mbA[1] * barcoi.ls/2 * (barcoi.wb**2 * comps[1]\
        - (barcoi.M  - barcoi.ls * barcoi.mua * barcoi.wb) * comps[0]/barcoi.Ib)
    else:
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
            M[0,0] = 1
            M[1,1] = 1
            B[0] = cad.Fi[0]
            B[1] = cad.Fi[1]
 ########################BARCO EN POSICION INTERMEDIA##########################       
    if barcoc:
        
        '''si hay barco  central el valor de barcoc sera el del
        eslabon siguiente al barco: barcoc.link = i implica que el barco 
        esta enganchado entre el eslabon i-1 y el i.'''
        
        i = barcoc.link
        #copiamos las filas de la matriz M desde 2i-2 hasta el final de modo
        #desplacemos todas las filas dos posiciones hacia abajo y
        #hacia la derecha para introducir
        #las tensiones del barco en su sitio
        
        #M[2*i+2:,2*i+2:] = M[2*i:-2,2*i:-2]
        #idem con los terminos independientes
        #B[2*i+2:] = B[2*i:-2]
        #orientacion de barcoc y normal al barco
        comps = np.array([np.cos(barcoc.theta),np.sin(barcoc.theta)])
        compt = np.array([comps[1],-comps[0]])
        mbA =np.array(\
        [barcoc.mb + barcoc.mA[0] * comps[0] - barcoc.mA[1] * comps[1],
        barcoc.mb + barcoc.mA[1] * comps[1] + barcoc.mA[1] * comps[0]])        
        
        #algunos calculos de cantidades que se repiten
        f1 = mbA * barcoc.ls ** 2 / barcoc.Ib /4
        f2 = mbA * cad.L**2 / cad.I
        f1c2 = f1 * comps[0] ** 2            
        f1s2 = f1 * comps[1] ** 2
        f1cs = f1 * comps[0] * comps[1]
        
        f2xx = f2 * cad.normal[0,i-1] ** 2
        f2yy = f2 * cad.normal[1,i-1] ** 2
        f2xy = f2 * cad.normal[0,i-1] * cad.normal[1,i-1]
        
        f2xxa = f2 * cad.normal[0,i-2] ** 2
        f2yya = f2 * cad.normal[1,i-2] ** 2
        f2xya = f2 * cad.normal[0,i-2] * cad.normal[1,i-2]
        
        #coeficientes de la matriz de tensiones correpondientes a la
        #tension entre es barco y el eslabon siguiente
        M[2*i:2*i+2,2*i-2:2*i+4] = np.array(\
        [[-f1s2[0] - 1,\
        +f1cs[0],\
        1 + f1s2[0] + f2xx[0] + mbA[0] * mcA[0,i-1],\
        -f1cs[0] + f2xy[0],\
        f2xx[0] - mbA[0] * mcA[0,i-1],\
        + f2xy[0]],\
        [f1cs[1],\
        -f1c2[1] - 1,\
        -f1cs[1] + f2xy[1],\
        1 + f1c2[1] + f2yy[1] + mbA[1] * mcA[1,i-1],\
        + f2xy[1],\
        f2yy[1] - mbA[1] * mcA[1,i-1]]])
        
        #terminos independientes correspondientes         
        B[2*i] = - barcoc.Fb[0]\
        + barcoc.mul * np.dot(barcoc.vb,comps) * comps[0]\
        + barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[0]\
        + barcoc.mul2 * np.sign(np.dot(barcoc.vb,comps))\
        * np.dot(barcoc.vb,comps)**2 * comps[0]\
        + barcoc.mut2 * barcoc.ls * np.sign(np.dot(barcoc.vb,compt))\
        * np.dot(barcoc.vb,compt)**2 * compt[0]\
        - mbA[0] * cad.L * cad.w[i-1]**2 * cad.para[0,i-1]\
        - mbA[0] * barcoc.ls/2 * barcoc.wb**2 * comps[0]\
        - mbA[0] * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[1]\
        + mbA[0] * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
        * comps[1]\
        + mbA[0] * cad.L/cad.I * cad.A * cad.normal[0,i-1] * cad.w[i-1]\
        - mbA[0] * factor1[0,i-1]
        
        B[2*i+1] = -barcoc.Fb[1]\
        + barcoc.mul * np.dot(barcoc.vb,comps) * comps[1]\
        + barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[1]\
        + barcoc.mul2 * np.sign(np.dot(barcoc.vb,comps))\
        * np.dot(barcoc.vb,comps)**2 * comps[1]\
        + barcoc.mut2 * barcoc.ls * np.sign(np.dot(barcoc.vb,compt))\
        * np.dot(barcoc.vb,compt)**2 * compt[1]\
        - mbA[1] * cad.L * cad.w[i-1]**2 * cad.para[1,i-1]\
        - mbA[1] * barcoc.ls/2 * barcoc.wb**2 * comps[1]\
        + mbA[1] * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[0]\
        - mbA[1] * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
        * comps[0]\
        + mbA[1] * cad.L/cad.I * cad.A * cad.normal[1,i-1] * cad.w[i-1]\
        - mbA[1] * factor1[1,i-1]
        
        
        #Coeficientes de la ecuación de la tension entre el eslabon anterior
        #y el barco c
        M[2*i-2:2*i,2*i-4:2*i+2] = np.array(\
        [[-f2xxa[0] + mbA[0] * mcA[0,i-2],\
        -f2xya[0],\
        -mbA[0] * mcA[0,i-2] - f2xxa[0] - f1s2[0] - 1,\
        - f2xya[0] + f1cs[0],\
        f1s2[0] + 1,\
        - f1cs[0]],\
        [-f2xya[1],\
        -f2yya[1] + mbA[1] * mcA[1,i-2],\
        -f2xya[1] + f1cs[1],\
        -mbA[1] * mcA[1,i-2] - f2yya[1] - f1c2[1] - 1,\
        - f1cs[1],\
        f1c2[1] + 1]])
        
        #terminos independientes correspondientes         
        B[2*i-2] = - barcoc.Fb[0]\
        + barcoc.mul * np.dot(barcoc.vb,comps) * comps[0]\
        + barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[0]\
        + barcoc.mul2 * np.sign(np.dot(barcoc.vb,comps))\
        * np.dot(barcoc.vb,comps)**2 * comps[0]\
        + barcoc.mut2 * barcoc.ls * np.sign(np.dot(barcoc.vb,compt))\
        * np.dot(barcoc.vb,compt)**2 * compt[0]\
        + mbA[0] * cad.L * cad.w[i-2]**2 * cad.para[0,i-2]\
        - mbA[0] * barcoc.ls/2 * barcoc.wb**2 * comps[0]\
        - mbA[0] * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[1]\
        + mbA[0] * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
        * comps[1]\
        - mbA[0] * cad.L/cad.I * cad.A * cad.normal[0,i-2] * cad.w[i-2]\
        - mbA[0] * factor1[0,i-2]
        
        
        B[2*i-1] = -barcoc.Fb[1]\
        + barcoc.mul * np.dot(barcoc.vb,comps) * comps[1]\
        + barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[1]\
        + barcoc.mul2 * np.sign(np.dot(barcoc.vb,comps))\
        * np.dot(barcoc.vb,comps)**2 * comps[1]\
        + barcoc.mut2 * barcoc.ls * np.sign(np.dot(barcoc.vb,compt))\
        * np.dot(barcoc.vb,compt)**2 * compt[1]\
        + mbA[1] * cad.L * cad.w[i-2]**2 * cad.para[1,i-2]\
        - mbA[1] * barcoc.ls/2 *barcoc.wb**2 * comps[1]\
        + mbA[1] * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[0]\
        - mbA[1] * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
        * comps[0]\
        - mbA[1] * cad.L/cad.I * cad.A * cad.normal[1,i-2] * cad.w[i-2] \
        - mbA[1] * factor1[1,i-2]

        
###############################################################################        
        
    if barcod:
       comps = np.array([np.cos(barcod.theta),np.sin(barcod.theta)])
       compt = np.array([comps[1],-comps[0]])
       mbA =np.array(\
       [barcod.mb + barcod.mA[0] * comps[0] - barcod.mA[1] * comps[1],
       barcod.mb + barcod.mA[0] * comps[1] + barcod.mA[1] * comps[0]])
       M[-2,-4:] =\
       mbA[0] * T1xa[0,-1],\
       mbA[0] * T1ya[-1],\
       - 1 - mbA[0] * (2 * mcA[0,-1] - T1xa[0,-1]\
       + (barcod.ls/2. * comps[1])**2/barcod.Ib),\
       + mbA[0] * (T1ya[-1] + barcod.ls**2/4. * comps[0] * comps[1]/barcod.Ib)
       
       M[-1,-4:] =\
       mbA[1] * T1ya[-1],\
       mbA[1] * T1xa[1,-1],\
       +mbA[1] * (T1ya[-1] + barcod.ls**2/4. * comps[0] * comps[1]/barcod.Ib),\
       - 1 - mbA[1] * (2 * mcA[1,-1] - T1xa[1,-1]\
       + (barcod.ls/2. * comps[0])**2/barcod.Ib)
      
       B[-2] = - barcod.Fb[0]\
       + barcod.mul * np.dot(barcod.vr,comps) * comps[0]\
       + barcod.mut * barcod.ls * np.dot(barcod.vr,compt) * compt[0]\
       + barcod.mul2 * np.sign(np.dot(barcod.vr,comps))\
       * np.dot(barcod.vr,comps)**2 * comps[0]\
       + barcod.mut2 * barcod.ls * np.sign(np.dot(barcod.vr,compt))\
       * np.dot(barcod.vr,compt)**2 * compt[0]\
       - mbA[0] * factor1[0,-1]\
       + mbA[0] * cad.L * cad.w[-1]\
       * (cad.para[0,-1]  * cad.w[-1] - cad.A * cad.normal[0,-1]/cad.I)\
       - mbA[0] * barcod.ls/2 * (barcod.wb**2 * comps[0]\
       + (barcod.M  - barcod.ls * barcod.mua * barcod.wb) * comps[1]/barcod.Ib) 
       
       B[-1] = - barcod.Fb[1]\
       + barcod.mul * np.dot(barcod.vr,comps) * comps[1]\
       + barcod.mut * barcod.ls * np.dot(barcod.vr,compt) * compt[1]\
       + barcod.mul2 * np.sign(np.dot(barcod.vr,comps))\
       * np.dot(barcod.vr,comps)**2 * comps[1]\
       + barcod.mut2 * barcod.ls * np.sign(np.dot(barcod.vr,compt))\
       * np.dot(barcod.vr,compt)**2 * compt[1]\
       - mbA[1] * factor1[1,-1]\
       + mbA[1] * cad.L * cad.w[-1]\
       * (cad.para[1,-1]  * cad.w[-1] - cad.A * cad.normal[1,-1]/cad.I)\
       + mbA[1] * barcod.ls/2 * (-barcod.wb**2 * comps[1]\
       + (barcod.M  - barcod.ls * barcod.mua * barcod.wb) * comps[0]/barcod.Ib)    
    else:
       M[-1,-1] = 1
       M[-2,-2] = 1
       B[-2] = cad.Fd[0]
       B[-1] = cad.Fd[1]
      
           
    T = np.linalg.linalg.solve(M,B)
#    print('dentro', T)        
    return M,B,T       