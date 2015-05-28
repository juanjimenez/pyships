# -*- coding: utf-8 -*-
"""
Created on Mon Apr 06 22:43:48 2015
Este fichero contiene las matrices de tensiones de los eslabones de la cadena 
calculadas para distintas situaciones cadenas, barcos en las puntas extremos,
libres, reel en un extremo, en fin que es un archivo muuuuuy importante
@author: juan
"""
import numpy as np
import boom
import reel_solo
from scipy.sparse import csr_matrix, dia_matrix
from scipy.sparse.linalg import spsolve

#creamos una matriz pero esto depende de la situación
def Matriz_t_ini(cad , barcoint = 0, reel = 0):
    ''' inicio de la matriz de tensiones de la cadena, si barciint es distinto
    de cero, se está inicializando la matriz suponiendo que hay un barco en el
    centro, con lo cual hay que calcular una tensión más. cad es un objeto de
    tipo cadena (ver script boom)
    El resto de situaciones presuponen que no hay barco en el centro (de momento)
    Si hay bobina en un extremo se añade con reel = 1 se añaden dos ecuaciones
    más para las tensiones en la bobina y el tip (la punto del boom que se esta
    soltando en ese momento'''     
    if barcoint == 0:    

        #construimos una matriz sparse (csr) vacia de dimension el 2*numero de
        #eslabomes mas dos. Usamos csr (compressed sparse row format) porque es
        #el que mejor se adecua a la estructura de la matriz de coeficientes
        #de los formatos con lo que trabaja spsolve...
       
        M = csr_matrix((2 * cad.esl + 2 +2 * reel,2 * cad.esl + 2 + 2 * reel)) 
        # Costruimos el puntero de elementos no nulos de la matriz
        M.indptr[1:3] = [4,8]
        M.indptr[3:-2] = np.arange(14,14+6*(M.indptr.size-5),6)
        M.indptr[-2:] = [M.indptr[-3] + 4, M.indptr[-3] + 8]      
        
        #cosntruimos los indices de columna de los elementos no nulos de la 
        #matriz
        M.indices= np.array([0,1,2,3,0,1,2,3,0,1,2,3,4,5,0,1,2,3,4,5])
        for i in range(4,M.indptr.size-4,2):
            M.indices= np.concatenate((M.indices,i-2 + M.indices[8:20]))
        M.indices= np.concatenate((M.indices, i + M.indices[0:8]))
        M.data = np.zeros(M.indices.size)
        B = np.zeros(2 * cad.esl + 2 + 2 * reel)
        T = np.zeros(2 * cad.esl + 2 + 2 * reel)
       
        
    else:       
     
        M = csr_matrix((2 * cad.esl +4,2 * cad.esl +4))        
         # Costruimos el puntero de elementos no nulos de la matriz
        M.indptr[1:3] = [4,8]
        M.indptr[3:-2] = np.arange(14,14+6*(M.indptr.size-5),6)
        M.indptr[-2:] = [M.indptr[-3] + 4, M.indptr[-3] + 8]      
        
        #cosntruimos los indices de columna de los elementos no nulos de la 
        #matriz
        M.indices= np.array([0,1,2,3,0,1,2,3,0,1,2,3,4,5,0,1,2,3,4,5])
        for i in range(4,M.indptr.size-4,2):
            M.indices= np.concatenate((M.indices,i-2 + M.indices[8:20]))
        M.indices= np.concatenate((M.indices, i + M.indices[0:8]))
        M.data = np.zeros(M.indices.size)
        B = np.zeros(2 * cad.esl +4)
        T = np.zeros(2 * cad.esl +4)
        
    
    return M, B, T 
        
        
def Matriz_t(M, B, T, cad, barcoi = None, barcod = None, barcoc = None):
    ''' Calculo de la matriz de coeficientes para calcular las tensiones
    soportadas por los eslabones de la cadena. los extremos de la cadena
    pueden estan libres por defecto, además se puede enganchar un objeto
    barcoi al extremo izquierdo de la cadena, un objeto barcod al extremo
    derecho de la cadena y un objeto barcoc a un eslabon intermedio
    cualquiera. En principio, el eslabon concreto al que se engancha barcoc
    viene determinado por el atributo barcoc.link.'''
    
    normal_2 = cad.mL2 * cad.normal**2 #contiene los cuadrados de las
    #componentes del vector normal a los eslabones multiplicados por el 
    #factor m*L**2
    normal_cr = -cad.mL2 * cad.normal[0] * cad .normal[1] #producto de 
    #terminos cruzados permite obtenern T3ya T1ya de matlab 
    #factores para construir la matriz de coeficientes
    T1xa = cad.I -normal_2#definido como en el de matlab,tambien sirve
    #para T3xa
    if barcoc:
       #la segunda fila es T3yb y T1yb
        T2xa =np.concatenate((-2. * cad.I - (normal_2[:,:barcoc.link-2] +\
        normal_2[:,1:barcoc.link-1]),\
        [[0],[0]],[[0],[0]],\
        -2. * cad.I - (normal_2[:,barcoc.link-1:-1] +\
        normal_2[:,barcoc.link:])),axis = 1) #definido como
        #en matlab
        # al mantener las dos componentes, la de arriba es T2xa y la de abajo 
        #es T2yb
        #es tambien T2xb
        T2ya = np.concatenate((normal_cr[:barcoc.link-2]\
        + normal_cr[1:barcoc.link-1],\
        [0],[0],\
        normal_cr[barcoc.link-1:-1]\
        + normal_cr[barcoc.link:]), axis = 1)
        #definido como en matlab
        #por tanto, es también T2xb. 
        
      
        for i in range(0,barcoc.link-2): #modificamos las filas de
        #la matriz m 
#            M[2*i+2:2*i+4,2*i:2*i+6] = np.array(\
#            [[T1xa[0,i],normal_cr[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
#            normal_cr[i+1]],[normal_cr[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
#            normal_cr[i+1], T1xa[1,i+1]]])
            
            M.data[M.indptr[2*i+2]:M.indptr[2*i+4]] = \
            T1xa[0,i],normal_cr[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
            normal_cr[i+1],normal_cr[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
            normal_cr[i+1], T1xa[1,i+1]
           
        
        
        
        
        for i in range(barcoc.link,cad.normal.shape[1]): #modificamos las filas de
        #la matriz m 
            M.data[M.indptr[2*i+2]:M.indptr[2*i+4]] = \
            T1xa[0,i-1],normal_cr[i-1],T2xa[0,i],T2ya[i],T1xa[0,i],\
            normal_cr[i],normal_cr[i-1],T1xa[1,i-1],T2ya[i],T2xa[1,i],\
            normal_cr[i],T1xa[1,i]
        #calculo de los terminos independientes B
        #factores para construir los terminos independientes
        
        factor1 = cad.I * (cad.s * np.abs(np.sum(cad.v * cad.normal,0)) +\
        cad.q * abs(np.sum(cad.v * cad.para,0))) * cad.v / cad.vmod
        
        factor2 = cad.m * cad.L * (cad.I * cad.w**2 * cad.para - cad.A *\
        cad.w * cad.normal) 
        
        #cori = 2 * cad.m * cad.I * cad.w * cad.v
         
        
        B[range(2,2*barcoc.link-3,2)] =\
        + factor1[0,1:barcoc.link-1] - factor1[0,0:barcoc.link-2]\
        + factor2[0,1:barcoc.link-1] + factor2[0,0:barcoc.link-2]\
        #- cori[1,1:barcoc.link-1] + cori[1,0:barcoc.link-2]
        B[range(3,2*barcoc.link-2,2)] =\
        + factor1[1,1:barcoc.link-1] - factor1[1,0:barcoc.link-2]\
        + factor2[1,1:barcoc.link-1] + factor2[1,0:barcoc.link-2]\
        #- cori[0,1:barcoc.link-1] + cori[0,0:barcoc.link-2]
        
        B[range(2*barcoc.link+2,B.shape[0]-3,2)] =\
        + factor1[0,barcoc.link:] - factor1[0,barcoc.link-1:-1]\
        + factor2[0,barcoc.link:] + factor2[0,barcoc.link-1:-1]\
        #- cori[1,barcoc.link:] + cori[1,barcoc.link-1:-1]
        B[range(2*barcoc.link+3,B.shape[0]-2,2)] = factor1[1,barcoc.link:]\
        -factor1[1,barcoc.link-1:-1]\
        + factor2[1,barcoc.link:] +factor2[1,barcoc.link-1:-1]
        #- cori[0,barcoc.link:] + cori[0,barcoc.link-1:-1]
    else:
        #la segunda fila es T3yb y T1yb
        T2xa =-2. * cad.I - (normal_2[:,:-1] + normal_2[:,1:]) #definido como
        #en matlab
        # al mantener las dos componentes, la de arriba es T2xa y la de abajo 
        #es T2yb
        #es tambien T2xb
        T2ya = normal_cr[:-1] + normal_cr[1:] #definido como en matlab
        #por tanto, es también T2xb. 
        
      
        for i in range(0,cad.normal.shape[1]-1): #modificamos las filas de
        #la matriz m 
            M.data[M.indptr[2*i+2]:M.indptr[2*i+4]] = \
            [T1xa[0,i],normal_cr[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
            normal_cr[i+1],normal_cr[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
            normal_cr[i+1], T1xa[1,i+1]]
    
        #calculo de los terminos independientes B
        #factores para construir los terminos independientes
        
        factor1 = cad.I * (cad.s * np.abs(np.sum(cad.v * cad.normal,0)) +\
        cad.q * abs(np.sum(cad.v * cad.para,0))) * cad.v / cad.vmod
        
        factor2 = cad.m * cad.L * (cad.I * cad.w**2 * cad.para - cad.A *\
        cad.w * cad.normal)
        #cori = 2 * cad.m * cad.I * cad.w * np.array([-cad.v[1],cad.v[0]])        
        
        B[range(2,B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
        + factor2[0,1:] +factor2[0,0:-1]\
        #+ cori[0,1:] - cori[0,0:-1]
        B[range(3,B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
        + factor2[1,1:] +factor2[1,0:-1]
        #+ cori[1,1:] - cori[1,0:-1]
    
    #añadimos la contribución de los barcos si existen o el 'cierre' de las
    # matrices de coeficientes si no existen...(no me gusta la condicional)
    #pero de momento hasta que me aclare es lo que hay...
    if barcoi:
       
        comps = np.array([np.cos(barcoi.theta),np.sin(barcoi.theta)])
        compt = np.array([comps[1],-comps[0]])
        
        M.data[M.indptr[0]:M.indptr[1]] = barcoi.mb + cad.m\
        + barcoi.mb * (normal_2[0,0]/cad.I\
        + cad.m * (barcoi.ls/2. * comps[1])**2/barcoi.Ib),\
        - barcoi.mb * (normal_cr[0]/cad.I + cad. m\
        * barcoi.ls**2/4. * comps[0] * comps[1]/barcoi.Ib),\
        barcoi.mb * (normal_2[0,0]/cad.I  - 1.),\
        - barcoi.mb * normal_cr[0]/cad.I
        
        M.data[M.indptr[1]:M.indptr[2]] = M[0,1],\
        barcoi.mb + cad.m\
        + barcoi.mb * (normal_2[1,0]/cad.I
        + cad.m * (barcoi.ls/2. * comps[0])**2/barcoi.Ib),\
        M[0,3],\
        barcoi.mb * (normal_2[1,0]/cad.I  - 1.)
         
        B[0] = -cad.m * barcoi.Fm * comps[0]\
        + cad.m * barcoi.mul * np.dot(barcoi.vb,comps) * comps[0]\
        + cad.m * barcoi.mut * barcoi.ls * np.dot(barcoi.vb,compt) * compt[0]\
        - barcoi.mb * (np.abs(np.dot(cad.v[:,0],cad.normal[:,0]))\
        * cad.s + np.abs(np.dot(cad.v[:,0],cad.para[:,0])) * cad.q)\
        * cad.v[0,0] / cad.vmod[0]\
        - cad.m * barcoi.mb * cad.L * cad.w[0]\
        * (cad.para[0,0]  * cad.w[0] - cad.A * cad.normal[0,0]/cad.I)\
        - cad.m * barcoi.mb * barcoi.ls/2 * (barcoi.wb**2 * comps[0]
        + (barcoi.M  - barcoi.ls * barcoi.mua * barcoi.wb)\
        * comps[1]/barcoi.Ib) 
        
        B[1] = -cad.m * barcoi.Fm * comps[1]\
        + cad.m * barcoi.mul * np.dot(barcoi.vb,comps) * comps[1]\
        + cad.m * barcoi.mut * barcoi.ls * np.dot(barcoi.vb,compt) * compt[1]\
        - barcoi.mb * (np.abs(np.dot(cad.v[:,0],cad.normal[:,0]))\
        * cad.s + np.abs(np.dot(cad.v[:,0],cad.para[:,0])) * cad.q)\
        * cad.v[1,0] / cad.vmod[0]\
        - cad.m * barcoi.mb * cad.L * cad.w[0]\
        * (cad.para[1,0]  * cad.w[0] - cad.A * cad.normal[1,0]/cad.I)\
        - cad.m * barcoi.mb * barcoi.ls/2 * (barcoi.wb**2 * comps[1]\
        - (barcoi.M  - barcoi.ls * barcoi.mua * barcoi.wb)\
        * comps[0]/barcoi.Ib)
    else:
        M[0,0] = 1
        M[1,1] = 1
        B[0] = cad.Fi[0]
        B[1] = cad.Fi[1]
        
    if barcoc:
        
        '''si hay barco  central el valor de barcoc sera el del
        eslabon anterior al barco: barcoc.link = i implica que el barco 
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
        
        #algunos calculos de cantidades que se repiten
        f1 = cad.m * barcoc.mb * barcoc.ls ** 2 / barcoc.Ib /4
        f2 = barcoc.mb * cad.mL2 / cad.I
        f1c2 = f1 * comps[0] ** 2            
        f1s2 = f1 * comps[1] ** 2
        f1cs = f1 * comps[0] * comps[1]
        
        f2xx = f2 * cad.normal[0,i-1] ** 2
        f2yy = f2 * cad.normal[1,i-1] ** 2
        f2xy = f2 * cad.normal[0,i-1] * cad.normal[1,i-1]
        
        f2xxa = f2 * cad.normal[0,i-2] ** 2
        f2yya = f2 * cad.normal[1,i-2] ** 2
        f2xya = f2 * cad.normal[0,i-2] * cad.normal[1,i-2]
        
                
        M.data[M.indptr[2*i-2]:M.indptr[2*i]] = \
        -f2xxa + barcoc.mb,\
        -f2xya,\
        -barcoc.mb - f2xxa - f1s2 - cad.m,\
        - f2xya + f1cs,\
        f1s2 + cad.m,\
        -f1cs,\
        -f2xya,\
        -f2yya + barcoc.mb,\
        -f2xya + f1cs,\
        -barcoc.mb - f2yya - f1c2 - cad.m,\
        - f1cs,\
        f1c2 + cad.m
        #coeficientes de la matriz de tensiones correpondientes a la
        #tension entre es barco y el eslabon siguiente
        M[2*i:2*i+2,2*i-2:2*i+4] = np.array(\
        [[-(f1s2 + cad.m),\
        +f1cs,\
        cad.m + f1s2 + f2xx + barcoc.mb,\
        -f1cs + f2xy,\
        f2xx - barcoc.mb,\
        + f2xy],\
        [f1cs,\
        -(f1c2 + cad.m),\
        -f1cs + f2xy,\
        cad.m + f1c2 + f2yy + barcoc.mb,\
        + f2xy,\
        f2yy - barcoc.mb]])
        
        M.data[M.indptr[2*i]:M.indptr[2*i+2]] = \
        -(f1s2 + cad.m),\
        +f1cs,\
        cad.m + f1s2 + f2xx + barcoc.mb,\
        -f1cs + f2xy,\
        f2xx - barcoc.mb,\
        + f2xy,\
        f1cs,\
        -(f1c2 + cad.m),\
        -f1cs + f2xy,\
        cad.m + f1c2 + f2yy + barcoc.mb,\
        + f2xy,\
        f2yy - barcoc.mb
        
        #terminos independientes correspondientes         
        B[2*i] = -cad.m * barcoc.Fm * comps[0]\
        + cad.m * barcoc.mul * np.dot(barcoc.vb,comps) * comps[0]\
        + cad.m * barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[0]\
        - cad.m * barcoc.mb * cad.L * cad.w[i-1]**2 * cad.para[0,i-1]\
        - cad.m * barcoc.mb * barcoc.ls/2 * barcoc.wb**2 * comps[0]\
        - cad.m * barcoc.mb * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[1]\
        + cad.m * barcoc.mb * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
        * comps[1]\
        + cad.m * barcoc.mb* cad.L/cad.I * cad.A * cad.normal[0,i-1] * cad.w[i-1]\
        - barcoc.mb * (np.abs(np.dot(cad.v[:,i-1],cad.normal[:,i-1]))\
        * cad.s + np.abs(np.dot(cad.v[:,i-1],cad.para[:,i-1])) * cad.q)\
        * cad.v[0,i-1] / cad.vmod[i-1]
        
        B[2*i+1] = -cad.m * barcoc.Fm * comps[1]\
        + cad.m * barcoc.mul * np.dot(barcoc.vb,comps) * comps[1]\
        + cad.m * barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[1]\
        - cad.m * barcoc.mb * cad.L * cad.w[i-1]**2 * cad.para[1,i-1]\
        - cad.m * barcoc.mb * barcoc.ls/2 * barcoc.wb**2 * comps[1]\
        + cad.m * barcoc.mb * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[0]\
        - cad.m * barcoc.mb * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
        * comps[0]\
        + cad.m * barcoc.mb * cad.L/cad.I * cad.A * cad.normal[1,i-1] * cad.w[i-1]\
        - barcoc.mb * (np.abs(np.dot(cad.v[:,i-1],cad.normal[:,i-1]))\
        * cad.s + np.abs(np.dot(cad.v[:,i-1],cad.para[:,i-1])) * cad.q)\
        * cad.v[1,i-1] / cad.vmod[i-1]
        
        
        
        
        #terminos independientes correspondientes (elsabon anterior)       
        B[2*i-2] = -cad.m * barcoc.Fm * comps[0]\
        + cad.m * barcoc.mul * np.dot(barcoc.vb,comps) * comps[0]\
        + cad.m * barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[0]\
        + cad.m * barcoc.mb * cad.L * cad.w[i-2]**2 * cad.para[0,i-2]\
        - cad.m * barcoc.mb * barcoc.ls/2 * barcoc.wb**2 * comps[0]\
        - cad.m * barcoc.mb * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[1]\
        + cad.m * barcoc.mb * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
        * comps[1]\
        - cad.m * barcoc.mb * cad.L/cad.I * cad.A * cad.normal[0,i-2] * cad.w[i-2]\
        - barcoc.mb * (np.abs(np.dot(cad.v[:,i-2],cad.normal[:,i-2]))\
        * cad.s + np.abs(np.dot(cad.v[:,i-2],cad.para[:,i-2])) * cad.q)\
        * cad.v[0,i-2] / cad.vmod[i-2]
        
        
        B[2*i-1] = -cad.m * barcoc.Fm * comps[1]\
        + cad.m * barcoc.mul * np.dot(barcoc.vb,comps) * comps[1]\
        + cad.m * barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[1]\
        + cad.m * barcoc.mb * cad.L * cad.w[i-2]**2 * cad.para[1,i-2]\
        - cad.m * barcoc.mb * barcoc.ls/2 *barcoc.wb**2 * comps[1]\
        + cad.m * barcoc.mb * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[0]\
        - cad.m * barcoc.mb * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
        * comps[0]\
        - cad.m * barcoc.mb * cad.L/cad.I * cad.A * cad.normal[1,i-2] * cad.w[i-2] \
        - barcoc.mb * (np.abs(np.dot(cad.v[:,i-2],cad.normal[:,i-2]))\
        * cad.s + np.abs(np.dot(cad.v[:,i-2],cad.para[:,i-2])) * cad.q)\
        * cad.v[1,i-2] / cad.vmod[i-2]
        
        
    if barcod:
       comps = np.array([np.cos(barcod.theta),np.sin(barcod.theta)])
       compt = np.array([comps[1],-comps[0]])
       
#       M[-2,-4:] = - barcod.mb * (normal_2[0,-1]/cad.I  - 1.),\
#       + barcod.mb * normal_cr[-1]/cad.I,\
#       - barcod.mb - cad.m - barcod.mb * (normal_2[0,-1]/cad.I
#       + cad.m * (barcod.ls/2. * comps[1])**2/barcod.Ib),\
#       + barcod.mb * (normal_cr[-1]/cad.I + cad. m\
#       * barcod.ls**2/4. * comps[0] * comps[1]/barcod.Ib)
       
       M.data[M.indptr[-3]:M.indptr[-2]] = \
       - barcod.mb * (normal_2[0,-1]/cad.I  - 1.),\
       + barcod.mb * normal_cr[-1]/cad.I,\
       - barcod.mb - cad.m - barcod.mb * (normal_2[0,-1]/cad.I
       + cad.m * (barcod.ls/2. * comps[1])**2/barcod.Ib),\
       + barcod.mb * (normal_cr[-1]/cad.I + cad. m\
       * barcod.ls**2/4. * comps[0] * comps[1]/barcod.Ib)
       
#       M[-1,-4:] = M[-2,-3],\
#       - barcod.mb * (normal_2[1,-1]/cad.I  - 1.),\
#       M[-2,-1],\
#       - barcod.mb - cad.m - barcod.mb * (normal_2[1,-1]/cad.I\
#       + cad.m * (barcod.ls**2/4. * comps[0])**2/barcod.Ib)
       
       M.data[M.indptr[-2]:M.indptr[-1]] = M[-2,-3],\
       - barcod.mb * (normal_2[1,-1]/cad.I  - 1.),\
       M[-2,-1],\
       - barcod.mb - cad.m - barcod.mb * (normal_2[1,-1]/cad.I\
       + cad.m * (barcod.ls**2/4. * comps[0])**2/barcod.Ib)
      
       B[-2] = -cad.m * barcod.Fm * comps[0]\
       + cad.m * barcod.mul * np.dot(barcod.vb,comps) * comps[0]\
       + cad.m * barcod.mut * barcod.ls * np.dot(barcod.vb,compt) * compt[0]\
       - barcod.mb * (np.abs(np.dot(cad.v[:,-1],cad.normal[:,-1]))\
       * cad.s + np.abs(np.dot(cad.v[:,-1],cad.para[:,-1])) * cad.q)\
       * cad.v[0,-1] / cad.vmod[-1]\
       + cad.m * barcod.mb * cad.L * cad.w[-1]\
       * (cad.para[0,-1]  * cad.w[-1] - cad.A\
       * cad.normal[0,-1]/cad.I) - cad.m * barcod.mb * barcod.ls/2\
       * (barcod.wb**2 * comps[0] + (barcod.M  - barcod.ls * barcod.mua\
       * barcod.wb) * comps[1]/barcod.Ib) 
       
       B[-1] = -cad.m * barcod.Fm * comps[1]\
       + cad.m * barcod.mul * np.dot(barcod.vb,comps) * comps[1]\
       + cad.m * barcod.mut * barcod.ls\
       * np.dot(barcod.vb,compt) * compt[1]\
       - barcod.mb * (np.abs(np.dot(cad.v[:,-1],cad.normal[:,-1]))\
       * cad.s + np.abs(np.dot(cad.v[:,-1],cad.para[:,-1])) * cad.q)\
       * cad.v[1,-1] / cad.vmod[-1]\
       + cad.m * barcod.mb * cad.L * cad.w[-1]\
       * (cad.para[1,-1]  * cad.w[-1]\
       - cad.A * cad.normal[1,-1]/cad.I) + cad.m * barcod.mb\
       * barcod.ls/2 * (-barcod.wb**2 * comps[1] + (barcod.M  - barcod.ls\
       * barcod.mua * barcod.wb) * comps[0]/barcod.Ib)
    
    else:
       M[-1,-1] = 1
       M[-2,-2] = 1
       B[-2] = cad.Fd[0]
       B[-1] = cad.Fd[1]
      
           
    T = spsolve(M,B)
#    print('dentro', T)        
    return M,B,T

def Matriz_t_2b(M, B, T, cad, barcoi, barcod):
    ''' Calculo de la matriz de coeficientes para calcular las tensiones
    soportadas por los eslabones de la cadena. Esta funcion es específica 
    para el caso en que la cadena es arrastrada solo por dos barcos que estan
    localizados en los EXTREMOS de la cadena... para otras situaciones emplear
    Matriz_t o alguna otra de las variantes de Matriz_t'''
    
    normal_2 = cad.mL2 * cad.normal**2 #contiene los cuadrados de las
    #componentes del vector normal a los eslabones multiplicados por el 
    #factor m*L**2
    normal_cr = -cad.mL2 * cad.normal[0] * cad .normal[1] #producto de 
    #terminos cruzados permite obtenern T3ya T1ya de matlab 
    #factores para construir la matriz de coeficientes
    T1xa = cad.I -normal_2#definido como en el de matlab,tambien sirve
    #para T3xa
    
    
    #la segunda fila es T3yb y T1yb
    T2xa =-2. * cad.I - (normal_2[:,:-1] + normal_2[:,1:]) #definido como
    #en matlab
    # al mantener las dos componentes, la de arriba es T2xa y la de abajo 
    #es T2yb
    #es tambien T2xb
    T2ya = normal_cr[:-1] + normal_cr[1:] #definido como en matlab
    #por tanto, es también T2xb. 
    
      
    for i in range(0,cad.normal.shape[1]-1): #modificamos las filas de
    #la matriz m 
        M.data[M.indptr[2*i+2]:M.indptr[2*i+4]] = \
        [T1xa[0,i],normal_cr[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        normal_cr[i+1],normal_cr[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        normal_cr[i+1], T1xa[1,i+1]]
    
    #calculo de los terminos independientes B
    #factores para construir los terminos independientes
    
    factor1 = cad.I * (cad.s * np.abs(np.sum(cad.v * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.v * cad.para,0))) * cad.v / cad.vmod
    
    factor2 = cad.m * cad.L * (cad.I * cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal)
    #cori = 2 * cad.m * cad.I * cad.w * np.array([-cad.v[1],cad.v[0]])        
    
    B[range(2,B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
    + factor2[0,1:] +factor2[0,0:-1]\
    #+ cori[0,1:] - cori[0,0:-1]
    B[range(3,B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
    + factor2[1,1:] +factor2[1,0:-1]
       
    
       
    comps = np.array([np.cos(barcoi.theta),np.sin(barcoi.theta)])
    compt = np.array([comps[1],-comps[0]])
    
    M.data[M.indptr[0]:M.indptr[1]] = barcoi.mb + cad.m\
    + barcoi.mb * (normal_2[0,0]/cad.I\
    + cad.m * (barcoi.ls/2. * comps[1])**2/barcoi.Ib),\
    - barcoi.mb * (normal_cr[0]/cad.I + cad. m\
    * barcoi.ls**2/4. * comps[0] * comps[1]/barcoi.Ib),\
    barcoi.mb * (normal_2[0,0]/cad.I  - 1.),\
    - barcoi.mb * normal_cr[0]/cad.I
    
    M.data[M.indptr[1]:M.indptr[2]] = M[0,1],\
    barcoi.mb + cad.m\
    + barcoi.mb * (normal_2[1,0]/cad.I
    + cad.m * (barcoi.ls/2. * comps[0])**2/barcoi.Ib),\
    M[0,3],\
    barcoi.mb * (normal_2[1,0]/cad.I  - 1.)
     
    B[0] = -cad.m * barcoi.Fm * comps[0]\
    + cad.m * barcoi.mul * np.dot(barcoi.vb,comps) * comps[0]\
    + cad.m * barcoi.mut * barcoi.ls * np.dot(barcoi.vb,compt) * compt[0]\
    - barcoi.mb * (np.abs(np.dot(cad.v[:,0],cad.normal[:,0]))\
    * cad.s + np.abs(np.dot(cad.v[:,0],cad.para[:,0])) * cad.q)\
    * cad.v[0,0] / cad.vmod[0]\
    - cad.m * barcoi.mb * cad.L * cad.w[0]\
    * (cad.para[0,0]  * cad.w[0] - cad.A * cad.normal[0,0]/cad.I)\
    - cad.m * barcoi.mb * barcoi.ls/2 * (barcoi.wb**2 * comps[0]
    + (barcoi.M  - barcoi.ls * barcoi.mua * barcoi.wb)\
    * comps[1]/barcoi.Ib) 
    
    B[1] = -cad.m * barcoi.Fm * comps[1]\
    + cad.m * barcoi.mul * np.dot(barcoi.vb,comps) * comps[1]\
    + cad.m * barcoi.mut * barcoi.ls * np.dot(barcoi.vb,compt) * compt[1]\
    - barcoi.mb * (np.abs(np.dot(cad.v[:,0],cad.normal[:,0]))\
    * cad.s + np.abs(np.dot(cad.v[:,0],cad.para[:,0])) * cad.q)\
    * cad.v[1,0] / cad.vmod[0]\
    - cad.m * barcoi.mb * cad.L * cad.w[0]\
    * (cad.para[1,0]  * cad.w[0] - cad.A * cad.normal[1,0]/cad.I)\
    - cad.m * barcoi.mb * barcoi.ls/2 * (barcoi.wb**2 * comps[1]\
    - (barcoi.M  - barcoi.ls * barcoi.mua * barcoi.wb)\
    * comps[0]/barcoi.Ib)
        
    
    comps = np.array([np.cos(barcod.theta),np.sin(barcod.theta)])
    compt = np.array([comps[1],-comps[0]])
  
    M.data[M.indptr[-3]:M.indptr[-2]] = \
    - barcod.mb * (normal_2[0,-1]/cad.I  - 1.),\
    + barcod.mb * normal_cr[-1]/cad.I,\
    - barcod.mb - cad.m - barcod.mb * (normal_2[0,-1]/cad.I
    + cad.m * (barcod.ls/2. * comps[1])**2/barcod.Ib),\
    + barcod.mb * (normal_cr[-1]/cad.I + cad. m\
    * barcod.ls**2/4. * comps[0] * comps[1]/barcod.Ib)
  
    M.data[M.indptr[-2]:M.indptr[-1]] = M[-2,-3],\
    - barcod.mb * (normal_2[1,-1]/cad.I  - 1.),\
    M[-2,-1],\
    - barcod.mb - cad.m - barcod.mb * (normal_2[1,-1]/cad.I\
    + cad.m * (barcod.ls**2/4. * comps[0])**2/barcod.Ib)
  
    B[-2] = -cad.m * barcod.Fm * comps[0]\
    + cad.m * barcod.mul * np.dot(barcod.vb,comps) * comps[0]\
    + cad.m * barcod.mut * barcod.ls * np.dot(barcod.vb,compt) * compt[0]\
    - barcod.mb * (np.abs(np.dot(cad.v[:,-1],cad.normal[:,-1]))\
    * cad.s + np.abs(np.dot(cad.v[:,-1],cad.para[:,-1])) * cad.q)\
    * cad.v[0,-1] / cad.vmod[-1]\
    + cad.m * barcod.mb * cad.L * cad.w[-1]\
    * (cad.para[0,-1]  * cad.w[-1] - cad.A\
    * cad.normal[0,-1]/cad.I) - cad.m * barcod.mb * barcod.ls/2\
    * (barcod.wb**2 * comps[0] + (barcod.M  - barcod.ls * barcod.mua\
    * barcod.wb) * comps[1]/barcod.Ib) 
   
    B[-1] = -cad.m * barcod.Fm * comps[1]\
    + cad.m * barcod.mul * np.dot(barcod.vb,comps) * comps[1]\
    + cad.m * barcod.mut * barcod.ls\
    * np.dot(barcod.vb,compt) * compt[1]\
    - barcod.mb * (np.abs(np.dot(cad.v[:,-1],cad.normal[:,-1]))\
    * cad.s + np.abs(np.dot(cad.v[:,-1],cad.para[:,-1])) * cad.q)\
    * cad.v[1,-1] / cad.vmod[-1]\
    + cad.m * barcod.mb * cad.L * cad.w[-1]\
    * (cad.para[1,-1]  * cad.w[-1]\
    - cad.A * cad.normal[1,-1]/cad.I) + cad.m * barcod.mb\
    * barcod.ls/2 * (-barcod.wb**2 * comps[1] + (barcod.M  - barcod.ls\
    * barcod.mua * barcod.wb) * comps[0]/barcod.Ib)
               
    T = spsolve(M,B)
    #    print('dentro', T)        
    return M,B,T

def Matriz_t_3b(M, B, T, cad, barcoi, barcod, barcoc):
    ''' Calculo de la matriz de coeficientes para calcular las tensiones
    soportadas por los eslabones de la cadena. En este caso se supone que hay
    tres barcos arrastrando la cadena: dos en los extremos y uno en una posicion 
    intermedia definida en barcoc.link, para otras configuraciones distintas
    emplear Matrz_t o alguna de sus variantes'''
    
    normal_2 = cad.mL2 * cad.normal**2 #contiene los cuadrados de las
    #componentes del vector normal a los eslabones multiplicados por el 
    #factor m*L**2
    normal_cr = -cad.mL2 * cad.normal[0] * cad .normal[1] #producto de 
    #terminos cruzados permite obtenern T3ya T1ya de matlab 
    #factores para construir la matriz de coeficientes
    T1xa = cad.I -normal_2#definido como en el de matlab,tambien sirve
    #para T3xa
    
       #la segunda fila es T3yb y T1ybcad
    T2xa =np.concatenate((-2. * cad.I - (normal_2[:,:barcoc.link-2] +\
    normal_2[:,1:barcoc.link-1]),\
    [[0],[0]],[[0],[0]],\
    -2. * cad.I - (normal_2[:,barcoc.link-1:-1] +\
    normal_2[:,barcoc.link:])),axis = 1) #definido como
    #en matlab
    # al mantener las dos componentes, la de arriba es T2xa y la de abajo 
    #es T2yb
    #es tambien T2xb
    T2ya = np.concatenate((normal_cr[:barcoc.link-2]\
    + normal_cr[1:barcoc.link-1],\
    [0],[0],\
    normal_cr[barcoc.link-1:-1]\
    + normal_cr[barcoc.link:]), axis = 1)
    #definido como en matlab
    #por tanto, es también T2xb. 
    
  
    for i in range(0,barcoc.link-2): #modificamos las filas de
    #la matriz m 
        
        M.data[M.indptr[2*i+2]:M.indptr[2*i+4]] = \
        T1xa[0,i],normal_cr[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        normal_cr[i+1],normal_cr[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        normal_cr[i+1], T1xa[1,i+1]
       
    
    
    
    
    for i in range(barcoc.link,cad.normal.shape[1]): #modificamos las filas de
    #la matriz m 
        M.data[M.indptr[2*i+2]:M.indptr[2*i+4]] = \
        T1xa[0,i-1],normal_cr[i-1],T2xa[0,i],T2ya[i],T1xa[0,i],\
        normal_cr[i],normal_cr[i-1],T1xa[1,i-1],T2ya[i],T2xa[1,i],\
        normal_cr[i],T1xa[1,i]
    #calculo de los terminos independientes B
    #factores para construir los terminos independientes
    
    factor1 = cad.I * (cad.s * np.abs(np.sum(cad.v * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.v * cad.para,0))) * cad.v / cad.vmod
    
    factor2 = cad.m * cad.L * (cad.I * cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal) 
     
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
            
    #barco izquierda        
    comps = np.array([np.cos(barcoi.theta),np.sin(barcoi.theta)])
    compt = np.array([comps[1],-comps[0]])
    
    M.data[M.indptr[0]:M.indptr[1]] = barcoi.mb + cad.m\
    + barcoi.mb * (normal_2[0,0]/cad.I\
    + cad.m * (barcoi.ls/2. * comps[1])**2/barcoi.Ib),\
    - barcoi.mb * (normal_cr[0]/cad.I + cad. m\
    * barcoi.ls**2/4. * comps[0] * comps[1]/barcoi.Ib),\
    barcoi.mb * (normal_2[0,0]/cad.I  - 1.),\
    - barcoi.mb * normal_cr[0]/cad.I
    
    M.data[M.indptr[1]:M.indptr[2]] = M[0,1],\
    barcoi.mb + cad.m\
    + barcoi.mb * (normal_2[1,0]/cad.I
    + cad.m * (barcoi.ls/2. * comps[0])**2/barcoi.Ib),\
    M[0,3],\
    barcoi.mb * (normal_2[1,0]/cad.I  - 1.)
     
    B[0] = -cad.m * barcoi.Fm * comps[0]\
    + cad.m * barcoi.mul * np.dot(barcoi.vb,comps) * comps[0]\
    + cad.m * barcoi.mut * barcoi.ls * np.dot(barcoi.vb,compt) * compt[0]\
    - barcoi.mb * (np.abs(np.dot(cad.v[:,0],cad.normal[:,0]))\
    * cad.s + np.abs(np.dot(cad.v[:,0],cad.para[:,0])) * cad.q)\
    * cad.v[0,0] / cad.vmod[0]\
    - cad.m * barcoi.mb * cad.L * cad.w[0]\
    * (cad.para[0,0]  * cad.w[0] - cad.A * cad.normal[0,0]/cad.I)\
    - cad.m * barcoi.mb * barcoi.ls/2 * (barcoi.wb**2 * comps[0]
    + (barcoi.M  - barcoi.ls * barcoi.mua * barcoi.wb)\
    * comps[1]/barcoi.Ib) 
    
    B[1] = -cad.m * barcoi.Fm * comps[1]\
    + cad.m * barcoi.mul * np.dot(barcoi.vb,comps) * comps[1]\
    + cad.m * barcoi.mut * barcoi.ls * np.dot(barcoi.vb,compt) * compt[1]\
    - barcoi.mb * (np.abs(np.dot(cad.v[:,0],cad.normal[:,0]))\
    * cad.s + np.abs(np.dot(cad.v[:,0],cad.para[:,0])) * cad.q)\
    * cad.v[1,0] / cad.vmod[0]\
    - cad.m * barcoi.mb * cad.L * cad.w[0]\
    * (cad.para[1,0]  * cad.w[0] - cad.A * cad.normal[1,0]/cad.I)\
    - cad.m * barcoi.mb * barcoi.ls/2 * (barcoi.wb**2 * comps[1]\
    - (barcoi.M  - barcoi.ls * barcoi.mua * barcoi.wb)\
    * comps[0]/barcoi.Ib)
    
    #barco central    
    i = barcoc.link
    #copiamos las filas de la matriz M desde 2i-2 hasta el final de modo
    #desplacemos todas las filas dos posiciones hacia abajo y
    #hacia la derecha para introducir
    #las tensiones del barco en su sitio
    
    
    comps = np.array([np.cos(barcoc.theta),np.sin(barcoc.theta)])
    compt = np.array([comps[1],-comps[0]])
    
    #algunos calculos de cantidades que se repiten
    f1 = cad.m * barcoc.mb * barcoc.ls ** 2 / barcoc.Ib /4
    f2 = barcoc.mb * cad.mL2 / cad.I
    f1c2 = f1 * comps[0] ** 2            
    f1s2 = f1 * comps[1] ** 2
    f1cs = f1 * comps[0] * comps[1]
    
    f2xx = f2 * cad.normal[0,i-1] ** 2
    f2yy = f2 * cad.normal[1,i-1] ** 2
    f2xy = f2 * cad.normal[0,i-1] * cad.normal[1,i-1]
    
    f2xxa = f2 * cad.normal[0,i-2] ** 2
    f2yya = f2 * cad.normal[1,i-2] ** 2
    f2xya = f2 * cad.normal[0,i-2] * cad.normal[1,i-2]
    

    
    M.data[M.indptr[2*i-2]:M.indptr[2*i]] = \
    -f2xxa + barcoc.mb,\
    -f2xya,\
    -barcoc.mb - f2xxa - f1s2 - cad.m,\
    - f2xya + f1cs,\
    f1s2 + cad.m,\
    -f1cs,\
    -f2xya,\
    -f2yya + barcoc.mb,\
    -f2xya + f1cs,\
    -barcoc.mb - f2yya - f1c2 - cad.m,\
    - f1cs,\
    f1c2 + cad.m
    #coeficientes de la matriz de tensiones correpondientes a la
    #tension entre es barco y el eslabon siguiente
    M[2*i:2*i+2,2*i-2:2*i+4] = np.array(\
    [[-(f1s2 + cad.m),\
    +f1cs,\
    cad.m + f1s2 + f2xx + barcoc.mb,\
    -f1cs + f2xy,\
    f2xx - barcoc.mb,\
    + f2xy],\
    [f1cs,\
    -(f1c2 + cad.m),\
    -f1cs + f2xy,\
    cad.m + f1c2 + f2yy + barcoc.mb,\
    + f2xy,\
    f2yy - barcoc.mb]])
    
    M.data[M.indptr[2*i]:M.indptr[2*i+2]] = \
    -(f1s2 + cad.m),\
    +f1cs,\
    cad.m + f1s2 + f2xx + barcoc.mb,\
    -f1cs + f2xy,\
    f2xx - barcoc.mb,\
    + f2xy,\
    f1cs,\
    -(f1c2 + cad.m),\
    -f1cs + f2xy,\
    cad.m + f1c2 + f2yy + barcoc.mb,\
    + f2xy,\
    f2yy - barcoc.mb
    
    #terminos independientes correspondientes         
    B[2*i] = -cad.m * barcoc.Fm * comps[0]\
    + cad.m * barcoc.mul * np.dot(barcoc.vb,comps) * comps[0]\
    + cad.m * barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[0]\
    - cad.m * barcoc.mb * cad.L * cad.w[i-1]**2 * cad.para[0,i-1]\
    - cad.m * barcoc.mb * barcoc.ls/2 * barcoc.wb**2 * comps[0]\
    - cad.m * barcoc.mb * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[1]\
    + cad.m * barcoc.mb * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
    * comps[1]\
    + cad.m * barcoc.mb* cad.L/cad.I * cad.A * cad.normal[0,i-1] * cad.w[i-1]\
    - barcoc.mb * (np.abs(np.dot(cad.v[:,i-1],cad.normal[:,i-1]))\
    * cad.s + np.abs(np.dot(cad.v[:,i-1],cad.para[:,i-1])) * cad.q)\
    * cad.v[0,i-1] / cad.vmod[i-1]
    
    B[2*i+1] = -cad.m * barcoc.Fm * comps[1]\
    + cad.m * barcoc.mul * np.dot(barcoc.vb,comps) * comps[1]\
    + cad.m * barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[1]\
    - cad.m * barcoc.mb * cad.L * cad.w[i-1]**2 * cad.para[1,i-1]\
    - cad.m * barcoc.mb * barcoc.ls/2 * barcoc.wb**2 * comps[1]\
    + cad.m * barcoc.mb * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[0]\
    - cad.m * barcoc.mb * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
    * comps[0]\
    + cad.m * barcoc.mb * cad.L/cad.I * cad.A * cad.normal[1,i-1] * cad.w[i-1]\
    - barcoc.mb * (np.abs(np.dot(cad.v[:,i-1],cad.normal[:,i-1]))\
    * cad.s + np.abs(np.dot(cad.v[:,i-1],cad.para[:,i-1])) * cad.q)\
    * cad.v[1,i-1] / cad.vmod[i-1]
    
    
    
    
    #terminos independientes correspondientes (elsabon anterior)       
    B[2*i-2] = -cad.m * barcoc.Fm * comps[0]\
    + cad.m * barcoc.mul * np.dot(barcoc.vb,comps) * comps[0]\
    + cad.m * barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[0]\
    + cad.m * barcoc.mb * cad.L * cad.w[i-2]**2 * cad.para[0,i-2]\
    - cad.m * barcoc.mb * barcoc.ls/2 * barcoc.wb**2 * comps[0]\
    - cad.m * barcoc.mb * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[1]\
    + cad.m * barcoc.mb * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
    * comps[1]\
    - cad.m * barcoc.mb * cad.L/cad.I * cad.A * cad.normal[0,i-2] * cad.w[i-2]\
    - barcoc.mb * (np.abs(np.dot(cad.v[:,i-2],cad.normal[:,i-2]))\
    * cad.s + np.abs(np.dot(cad.v[:,i-2],cad.para[:,i-2])) * cad.q)\
    * cad.v[0,i-2] / cad.vmod[i-2]
    
    
    B[2*i-1] = -cad.m * barcoc.Fm * comps[1]\
    + cad.m * barcoc.mul * np.dot(barcoc.vb,comps) * comps[1]\
    + cad.m * barcoc.mut * barcoc.ls * np.dot(barcoc.vb,compt) * compt[1]\
    + cad.m * barcoc.mb * cad.L * cad.w[i-2]**2 * cad.para[1,i-2]\
    - cad.m * barcoc.mb * barcoc.ls/2 *barcoc.wb**2 * comps[1]\
    + cad.m * barcoc.mb * barcoc.ls/2/barcoc.Ib * barcoc.M * comps[0]\
    - cad.m * barcoc.mb * barcoc.ls**2/2/barcoc.Ib * barcoc.mua * barcoc.wb\
    * comps[0]\
    - cad.m * barcoc.mb * cad.L/cad.I * cad.A * cad.normal[1,i-2] * cad.w[i-2] \
    - barcoc.mb * (np.abs(np.dot(cad.v[:,i-2],cad.normal[:,i-2]))\
    * cad.s + np.abs(np.dot(cad.v[:,i-2],cad.para[:,i-2])) * cad.q)\
    * cad.v[1,i-2] / cad.vmod[i-2]
    
    
    #barco d
    comps = np.array([np.cos(barcod.theta),np.sin(barcod.theta)])
    compt = np.array([comps[1],-comps[0]])
        
    M.data[M.indptr[-3]:M.indptr[-2]] = \
    - barcod.mb * (normal_2[0,-1]/cad.I  - 1.),\
    + barcod.mb * normal_cr[-1]/cad.I,\
    - barcod.mb - cad.m - barcod.mb * (normal_2[0,-1]/cad.I
    + cad.m * (barcod.ls/2. * comps[1])**2/barcod.Ib),\
    + barcod.mb * (normal_cr[-1]/cad.I + cad. m\
    * barcod.ls**2/4. * comps[0] * comps[1]/barcod.Ib)
          
    M.data[M.indptr[-2]:M.indptr[-1]] = M[-2,-3],\
    - barcod.mb * (normal_2[1,-1]/cad.I  - 1.),\
    M[-2,-1],\
    - barcod.mb - cad.m - barcod.mb * (normal_2[1,-1]/cad.I\
    + cad.m * (barcod.ls**2/4. * comps[0])**2/barcod.Ib)
  
    B[-2] = -cad.m * barcod.Fm * comps[0]\
    + cad.m * barcod.mul * np.dot(barcod.vb,comps) * comps[0]\
    + cad.m * barcod.mut * barcod.ls * np.dot(barcod.vb,compt) * compt[0]\
    - barcod.mb * (np.abs(np.dot(cad.v[:,-1],cad.normal[:,-1]))\
    * cad.s + np.abs(np.dot(cad.v[:,-1],cad.para[:,-1])) * cad.q)\
    * cad.v[0,-1] / cad.vmod[-1]\
    + cad.m * barcod.mb * cad.L * cad.w[-1]\
    * (cad.para[0,-1]  * cad.w[-1] - cad.A\
    * cad.normal[0,-1]/cad.I) - cad.m * barcod.mb * barcod.ls/2\
    * (barcod.wb**2 * comps[0] + (barcod.M  - barcod.ls * barcod.mua\
    * barcod.wb) * comps[1]/barcod.Ib) 
   
    B[-1] = -cad.m * barcod.Fm * comps[1]\
    + cad.m * barcod.mul * np.dot(barcod.vb,comps) * comps[1]\
    + cad.m * barcod.mut * barcod.ls\
    * np.dot(barcod.vb,compt) * compt[1]\
    - barcod.mb * (np.abs(np.dot(cad.v[:,-1],cad.normal[:,-1]))\
    * cad.s + np.abs(np.dot(cad.v[:,-1],cad.para[:,-1])) * cad.q)\
    * cad.v[1,-1] / cad.vmod[-1]\
    + cad.m * barcod.mb * cad.L * cad.w[-1]\
    * (cad.para[1,-1]  * cad.w[-1]\
    - cad.A * cad.normal[1,-1]/cad.I) + cad.m * barcod.mb\
    * barcod.ls/2 * (-barcod.wb**2 * comps[1] + (barcod.M  - barcod.ls\
    * barcod.mua * barcod.wb) * comps[0]/barcod.Ib)
    
      
           
    T = spsolve(M,B)
#    print('dentro', T)        
    return M,B,T

def Matriz_t1(M, B, T, cad, barcoi):
    ''' Calculo de la matriz de coeficientes para calcular las tensiones
    soportadas por los eslabones de la cadena. Solo hay barco en un extremo
    que por defecto se toma a la izquierda de la posición inicial de la 
    cadena. El otro extremo queda libre o sometido a una fuerza exterior 
    asociada a la cadena.'''
    
    normal_2 = cad.mL2 * cad.normal**2 #contiene los cuadrados de las
    #componentes del vector normal a los eslabones multiplicados por el 
    #factor m*L**2
    normal_cr = -cad.mL2 * cad.normal[0] * cad .normal[1] #producto de 
    #terminos cruzados permite obtenern T3ya T1ya de matlab 
    #factores para construir la matriz de coeficientes
    T1xa = cad.I -normal_2#definido como en el de matlab,tambien sirve
    #para T3xa
    #la segunda fila es T3yb y T1yb
    T2xa =-2. * cad.I - (normal_2[:,:-1] + normal_2[:,1:]) #definido como
    #en matlab
    # al mantener las dos componentes, la de arriba es T2xa y la de abajo 
    #es T2yb
    #es tambien T2xb
    T2ya = normal_cr[:-1] + normal_cr[1:] #definido como en matlab
    #por tanto, es también T2xb. 
    
  
    for i in range(0,cad.normal.shape[1]-1): #modificamos las filas de
    #la matriz m 
        M.data[M.indptr[2*i+2]:M.indptr[2*i+4]] = \
        [T1xa[0,i],normal_cr[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        normal_cr[i+1],normal_cr[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        normal_cr[i+1], T1xa[1,i+1]]

    #calculo de los terminos independientes B
    #factores para construir los terminos independientes
    
    factor1 = cad.I * (cad.s * np.abs(np.sum(cad.v * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.v * cad.para,0))) * cad.v / cad.vmod
    
    factor2 = cad.m * cad.L * (cad.I * cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal)
    #cori = 2 * cad.m * cad.I * cad.w * np.array([-cad.v[1],cad.v[0]])        
    
    B[range(2,B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
    + factor2[0,1:] +factor2[0,0:-1]\
    #+ cori[0,1:] - cori[0,0:-1]
    B[range(3,B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
    + factor2[1,1:] +factor2[1,0:-1]
    #+ cori[1,1:] - cori[1,0:-1]
    # barcoi
    comps = np.array([np.cos(barcoi.theta),np.sin(barcoi.theta)])
    compt = np.array([comps[1],-comps[0]])
    
    M.data[M.indptr[0]:M.indptr[1]] = barcoi.mb + cad.m\
    + barcoi.mb * (normal_2[0,0]/cad.I\
    + cad.m * (barcoi.ls/2. * comps[1])**2/barcoi.Ib),\
    - barcoi.mb * (normal_cr[0]/cad.I + cad. m\
    * barcoi.ls**2/4. * comps[0] * comps[1]/barcoi.Ib),\
    barcoi.mb * (normal_2[0,0]/cad.I  - 1.),\
    - barcoi.mb * normal_cr[0]/cad.I
    
    M.data[M.indptr[1]:M.indptr[2]] = M[0,1],\
    barcoi.mb + cad.m\
    + barcoi.mb * (normal_2[1,0]/cad.I
    + cad.m * (barcoi.ls/2. * comps[0])**2/barcoi.Ib),\
    M[0,3],\
    barcoi.mb * (normal_2[1,0]/cad.I  - 1.)
     
    B[0] = -cad.m * barcoi.Fm * comps[0]\
    + cad.m * barcoi.mul * np.dot(barcoi.vb,comps) * comps[0]\
    + cad.m * barcoi.mut * barcoi.ls * np.dot(barcoi.vb,compt) * compt[0]\
    - barcoi.mb * (np.abs(np.dot(cad.v[:,0],cad.normal[:,0]))\
    * cad.s + np.abs(np.dot(cad.v[:,0],cad.para[:,0])) * cad.q)\
    * cad.v[0,0] / cad.vmod[0]\
    - cad.m * barcoi.mb * cad.L * cad.w[0]\
    * (cad.para[0,0]  * cad.w[0] - cad.A * cad.normal[0,0]/cad.I)\
    - cad.m * barcoi.mb * barcoi.ls/2 * (barcoi.wb**2 * comps[0]
    + (barcoi.M  - barcoi.ls * barcoi.mua * barcoi.wb)\
    * comps[1]/barcoi.Ib) 
    
    B[1] = -cad.m * barcoi.Fm * comps[1]\
    + cad.m * barcoi.mul * np.dot(barcoi.vb,comps) * comps[1]\
    + cad.m * barcoi.mut * barcoi.ls * np.dot(barcoi.vb,compt) * compt[1]\
    - barcoi.mb * (np.abs(np.dot(cad.v[:,0],cad.normal[:,0]))\
    * cad.s + np.abs(np.dot(cad.v[:,0],cad.para[:,0])) * cad.q)\
    * cad.v[1,0] / cad.vmod[0]\
    - cad.m * barcoi.mb * cad.L * cad.w[0]\
    * (cad.para[1,0]  * cad.w[0] - cad.A * cad.normal[1,0]/cad.I)\
    - cad.m * barcoi.mb * barcoi.ls/2 * (barcoi.wb**2 * comps[1]\
    - (barcoi.M  - barcoi.ls * barcoi.mua * barcoi.wb)\
    * comps[0]/barcoi.Ib)
    
#extremo libre        
    M[-1,-1] = 1
    M[-2,-2] = 1
    B[-2] = cad.Fd[0]
    B[-1] = cad.Fd[1]
  
           
    T = spsolve(M,B)
#    print('dentro', T)        
    return M,B,T

def Matriz_trb(M, B, T, cad,bob, barcod):
    ''' Calculo de la matriz de coeficientes para calcular las tensiones
    soportadas por los eslabones de la cadena. La matriz supone que en un 
    extremo de la cadena tenemos enganchada una bobina y en el otro un barco
    que tira. Se considera un barco a la derecha de la cadena
    y la bobina a la izquierda de la posicion inicial. '''
    
    normal_2 = cad.mL2 * cad.normal**2 #contiene los cuadrados de las
    #componentes del vector normal a los eslabones multiplicados por el 
    #factor m*L**2
    normal_cr = -cad.mL2 * cad.normal[0] * cad .normal[1] #producto de 
    #terminos cruzados permite obtenern T3ya T1ya de matlab 
    #factores para construir la matriz de coeficientes
    T1xa = cad.I -normal_2#definido como en el de matlab,tambien sirve
    #para T3xa
    
    #la segunda fila es T3yb y T1yb
    T2xa =-2. * cad.I - (normal_2[:,:-1] + normal_2[:,1:]) #definido como
    #en matlab
    # al mantener las dos componentes, la de arriba es T2xa y la de abajo 
    #es T2yb
    #es tambien T2xb
    T2ya = normal_cr[:-1] + normal_cr[1:] #definido como en matlab
    #por tanto, es también T2xb. 
        
      
    for i in range(0,cad.normal.shape[1]-1): #modificamos las filas de
    #la matriz m 
        M.data[M.indptr[2*i+4]:M.indptr[2*i+6]] = \
        [T1xa[0,i],normal_cr[i],T2xa[0,i],T2ya[i],T1xa[0,i+1],\
        normal_cr[i+1],normal_cr[i],T1xa[1,i],T2ya[i],T2xa[1,i],\
        normal_cr[i+1], T1xa[1,i+1]]
    
    #calculo de los terminos independientes B
    #factores para construir los terminos independientes
        
    factor1 = cad.I * (cad.s * np.abs(np.sum(cad.v * cad.normal,0)) +\
    cad.q * abs(np.sum(cad.v * cad.para,0))) * cad.v / cad.vmod
        
    factor2 = cad.m * cad.L * (cad.I * cad.w**2 * cad.para - cad.A *\
    cad.w * cad.normal)
              
        
    B[range(4,B.shape[0]-2,2)] = factor1[0,1:]-factor1[0,0:-1]\
    + factor2[0,1:] +factor2[0,0:-1]\
        
    B[range(3,B.shape[0]-2,2)] = factor1[1,1:]-factor1[1,0:-1]\
    + factor2[1,1:] +factor2[1,0:-1]
        
    
    
    #terminos del reel y la punta
    #primer par de ecuaciones (dinamica del reel)
    sb = np.sin(bob.theta)
    cb = np.cos(bob.theta)
    st = np.sin(bob.tht)
    ct = np.cos(bob.tht)
    pr1 = - (bob.empty_ra * sb - bob.cad.dh / np.pi \
    * (cb - bob.theta * sb / 2))        
    pr2 = (bob.empty_ra + bob.cad.dh * bob.theta / 2 / np.pi) / bob.I        
    pr3 = - bob.lt**2 * st / bob.It        
    pr4 = np.cos(bob.tht)**2 + 1
    pr5 = np.cos(bob.tht) * st
        
    pr6 = - (bob.empty_ra * cb - bob.cad.dh / np.pi \
    * (sb + bob.theta * cb / 2))            
    pr7 = - bob.lt**2 * ct / bob.It
    pr8 = st**2 + 1        
        
    M.data[M.indptr[0]:M.indptr[1]] = pr1 * - sb * pr2\
    - pr3 * st+ pr4,\
    pr1 * cb * pr2 + pr3 * ct + pr5,\
    pr3 * st - pr4,\
    pr3 * ct - pr5 
        
    M.data[M.indptr[1]:M.indptr[2]] = pr6 * - sb * pr2\
    - pr7 * st+ pr4,\
    pr6 * cb * pr2 + pr7 * ct + pr8,\
    pr7 * -st - pr5,\
    pr7 * ct - pr8
         
    B[0] =  (bob.empty_ra * cb + bob.cad.dh / np.pi \
    * (sb + bob.theta * cb / 2)) * bob.w ** 2\
    - pr1 * bob.Ar * bob.w / bob.I \
    + bob.lt * ct * bob.wt ** 2 \
    + bob.lt * st * bob.wt / bob.It \
    + 2 * (bob.vt[0] * ct + bob.vt[1] * st) * st * bob.wt
        
    B[1] = - (bob.empty_ra * cb + bob.cad.dh / np.pi \
    * (sb + bob.theta * cb / 2)) * bob.w ** 2\
    - pr6 * bob.Ar * bob.w / bob.I \
    + bob.lt * st * bob.w ** 2 \
    + 2 * (bob.vt[0] * ct + bob.vt[1] * st) * ct * bob.wt
    
    #segundo par de equaciones (dinamica del tip)
    pr1t = -cad.m * cad.I * bob.It + cad.m * bob.mt * cad.I * bob.lt ** 2\
    * np.array([st ** 2, ct ** 2])
    pr2t = - cad.m * bob.mt * cad.I * bob.lt ** 2 * ct * st
    pr3t = cad.m * cad.I * bob.It + cad.m * bob.mt * cad.I * bob.lt ** 2\
    * np.array([st ** 2, ct ** 2])\
    + cad.I * bob.mt * cad.L ** 2 * cad.normal[:,0] ** 2
    pr4t = - cad.m * bob.mt * cad.I * bob.lt ** 2 * ct * st \
    + bob.mt * cad.I * bob.lt ** 2 * cad.normal[0,0] * cad.normal[0,1]
    pr5t = bob.It * bob.mt * normal_2[0] - bob.It * cad.I * bob.mt
    pr6t = - bob.It * bob.mt * normal_cr[0]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    
    M.data[M.indptr[2]:M.indptr[3]] = pr1t[0], pr2t, pr3t[0], pr4t, pr5t[0],\
    pr6t
    
    M.data[M.indptr[3]:M.indptr[4]] = pr2t, pr1t[1], pr4t, pr3t[1], pr6t,\
    pr5t[1]  
    
    bp1 = - cad.m * cad.I * bob.It * bob.mt * bob.lt * bob.wt *2\
    * np.array([ct,st])
    bp2 = cad.m * cad.I * bob.mt * bob.lt * bob.At * bob.wt * np.array([-st,ct])
    bp3 = - cad.m * cad.I * bob.It * bob.mt * cad.L * cad.w[0] ** 2\
    * cad.para[:,0]
    bp4 = - cad.m * bob.It * bob.mt * cad.L * cad.A * cad.w[0]\
    * cad.normal[:,0]
    bp5 = - bob.mt * bob.It * factor1[:,0]
    
    B[2] =  bp1[0] + bp2[0] + bp3[0] + bp4[0] + bp5[0]
    B[3] = bp1[1] + bp2[1] + bp3[1] + bp4[1] + bp5[1]
    
        
   # barcod
    comps = np.array([np.cos(barcod.theta),np.sin(barcod.theta)])
    compt = np.array([comps[1],-comps[0]])              
    M.data[M.indptr[-3]:M.indptr[-2]] = \
    - barcod.mb * (normal_2[0,-1]/cad.I  - 1.),\
    + barcod.mb * normal_cr[-1]/cad.I,\
    - barcod.mb - cad.m - barcod.mb * (normal_2[0,-1]/cad.I
    + cad.m * (barcod.ls/2. * comps[1])**2/barcod.Ib),\
    + barcod.mb * (normal_cr[-1]/cad.I + cad. m\
    * barcod.ls**2/4. * comps[0] * comps[1]/barcod.Ib)
              
    M.data[M.indptr[-2]:M.indptr[-1]] = M[-2,-3],\
    - barcod.mb * (normal_2[1,-1]/cad.I  - 1.),\
    M[-2,-1],\
    - barcod.mb - cad.m - barcod.mb * (normal_2[1,-1]/cad.I\
    + cad.m * (barcod.ls**2/4. * comps[0])**2/barcod.Ib)
      
    B[-2] = -cad.m * barcod.Fm * comps[0]\
    + cad.m * barcod.mul * np.dot(barcod.vb,comps) * comps[0]\
    + cad.m * barcod.mut * barcod.ls * np.dot(barcod.vb,compt) * compt[0]\
    - barcod.mb * (np.abs(np.dot(cad.v[:,-1],cad.normal[:,-1]))\
    * cad.s + np.abs(np.dot(cad.v[:,-1],cad.para[:,-1])) * cad.q)\
    * cad.v[0,-1] / cad.vmod[-1]\
    + cad.m * barcod.mb * cad.L * cad.w[-1]\
    * (cad.para[0,-1]  * cad.w[-1] - cad.A\
    * cad.normal[0,-1]/cad.I) - cad.m * barcod.mb * barcod.ls/2\
    * (barcod.wb**2 * comps[0] + (barcod.M  - barcod.ls * barcod.mua\
    * barcod.wb) * comps[1]/barcod.Ib) 
       
    B[-1] = -cad.m * barcod.Fm * comps[1]\
    + cad.m * barcod.mul * np.dot(barcod.vb,comps) * comps[1]\
    + cad.m * barcod.mut * barcod.ls\
    * np.dot(barcod.vb,compt) * compt[1]\
    - barcod.mb * (np.abs(np.dot(cad.v[:,-1],cad.normal[:,-1]))\
    * cad.s + np.abs(np.dot(cad.v[:,-1],cad.para[:,-1])) * cad.q)\
    * cad.v[1,-1] / cad.vmod[-1]\
    + cad.m * barcod.mb * cad.L * cad.w[-1]\
    * (cad.para[1,-1]  * cad.w[-1]\
    - cad.A * cad.normal[1,-1]/cad.I) + cad.m * barcod.mb\
    * barcod.ls/2 * (-barcod.wb**2 * comps[1] + (barcod.M  - barcod.ls\
    * barcod.mua * barcod.wb) * comps[0]/barcod.Ib)
    
    
      
           
    T = spsolve(M,B)
#    print('dentro', T)        
    return M,B,T              