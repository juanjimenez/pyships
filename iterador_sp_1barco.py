# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:08:57 2014
Este programa esta pensado para iterar cadenas y barcos y ver si peta o que
pasa

@author: juan
"""
import barcosolo
import boom
#import barchain_matrix
import barchain_matrix_sp_all
import numpy as np
from matplotlib import pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches

delta = 0.01
#creamos los barcos
bi = barcosolo.barco()


#creamos la cadena de momento pequeñita porque es para probar
cd5 = boom.cadena(20)

#situamos la cadena ... toda en el eje x con los el primer y último eslabon.
#apuntando para arriba.

cd5.normal[:,0] = [1,0]
cd5.normal[:,0] = [1,0]
cd5.normal[:,-1] = [-1,0]
cd5.para[:,0] = [0,1]
cd5.para[:,-1] = [0,-1]
cd5.calcms()
cd5.A = 10
cd5.s = 20
#inicializamos los barquitos,poniendolos al principio y al final de la cadena
bi.pb[:] = [0,1.5]
bi.Fm = 500
bi.M = 0
nudoss = np.zeros([2,cd5.esl+1]) #posición del extremo derecho de cada eslabón

nudosa = np.zeros([2,cd5.esl+1]) #posisción del extremo izquierdo de cada eslabón
#inicializamos la matriz de coeficientes de las tensiones en los extremos
#de los eslabones de las cadenas.
errores = np.zeros([2,cd5.esl+1])
maxerr = np.zeros([2,cd5.esl+1])
#cd5.Matriz_t(bi,bd)
#movemos los elabones, estamos usando el paso de integracion por defecto
#delta = 5e-4
barras = np.zeros((2,cd5.cms.shape[1]+1))
tam = 100000#aqui definimos de momento el tamaño del experimento...
step =1000
# y aqui cada cuantas iteraciones guardamos datos y pintamos
birec = np.zeros((tam//step+1,16)) #el barco izquierdo guarda los valores
#de la tensión con el primer eslabon... (chapuzilla pero en fin) hay
#una tension más que eslabones

cadrec = np.zeros((tam//step+1,14,cd5.esl))
#cd5.Fd =[0., 10.]poder
M, B, T = barchain_matrix_sp_all.Matriz_t_ini(cd5,0)

for i in range(tam):
    #print i, i%step
    #M, B, T = barchain_matrix.Matriz_t_ini(cd5,0)
    M, B, T = barchain_matrix_sp_all.Matriz_t1(M, B, T, cd5,bi)
    cd5.T[:,0] = T[:-3:2]
    cd5.T[:,1] = T[1:-2:2]
    cd5.T[:,2] = T[2:-1:2]
    cd5.T[:,3] = T[3::2]
    cd5.movimiento(delta = delta)
   
    #y a continuación movemos los barcos
    #bi.controlador(pi/3,5,delta)
    bi.movimientotst(extf = T[0:2],dfcm = [-bi.ls/2.,0],delta = delta)
    #bd.controlador(pi/3,5,delta)    
    
#    cd5.Matriz_t(bi,bd)
    
    
      
#    #famoso parche para corregir errores de redondeo
#    popabi = bi.pb - bi.ls/2*np.array([np.cos(bi.theta),np.sin(bi.theta)]) 
#    primero = popabi - cd5.L*cd5.para[:,0]
#    ultimo = cd5.L*cd5.para[:,-1]
#    
#    
#
#    #aprovechamos para evaluar los errores cometidos, antes de recolocar los 
#    #eslabones. Esto solo tiene sentido en fase de depuración luego se quita.
#    nudoss[:,0] = popabi
#    nudosa[:,-1] = popabd
#    nudoss[:,1:] =   cd5.cms - cd5.L*cd5.para
#    nudosa[:,:-1] =  cd5.cms + cd5.L*cd5.para
#    errores = nudoss - nudosa
#    grandes = abs(errores) > abs(maxerr)
#    maxerr = grandes * errores + (1-grandes) * maxerr 
#        
#    #recoloco los eslabones...
#    cd5.calcms(primero,ultimo,-1*cd5.ord)
#    #recoloco los barcos... en realidad solo haría falta recolocar uno...    
#    bi.pb = cd5.cms[:,0]  + bi.ls/2*np.array([np.cos(bi.theta),np.sin(bi.theta)]) \
#    + cd5.L*cd5.para[:,0] 
    
    if i%step == 0:
        #print 'pinto'
        #pinto los errores maximos cometidos en los pasos que marca step
        # esto solo tiene sentido para depurar pero en una versión final
#        figure(1)
#        plot(arange(cd5.esl+1),errores[0,:],'o')
#        hold(True)
#        figure(2)        
#        plot(arange(cd5.esl+1),errores[1,:],'^')
#        hold(True)
#        pause(0.01)
         
        #guadamomos datos.
        birec[i//step] = bi.pb[0],bi.pb[1],bi.vb[0],bi.vb[1],bi.ab[0],bi.ab[1],\
        bi.theta,bi.wb,bi.alfab,bi.Fm,bi.M,bi.setl,bi.setw,bi.thewj,T[0],\
        T[1]

        
        

        cadrec[i//step] = cd5.cms[0],cd5.cms[1],cd5.v[0],cd5.v[1],cd5.a[0],\
        cd5.a[1],cd5.alfa,cd5.w,T[2:-1:2],T[2::2],cd5.normal[0],\
        cd5.normal[1],cd5.para[0],cd5.para[1]
#en la ultima fila guardamos propidades de los barcos y las cadenas y detalles
#del experimento       
        birec[-1,0:12] = bi.thewjmax,bi.Ac,bi.mb,bi.mul,bi.Ib,bi.mua,bi.ls,bi.mut,\
        bi.pmax,bi.thewjmax,bi.Ab,bi.Aw
        
        
        
        cadrec[-1,0,0:3] = cd5.s,cd5.q,cd5.A
        cadrec[-1,1,0:3] = cd5.L,cd5.m,cd5.I
        cadrec[-1,2,0:2] = delta,tam
np.savez('tugofwar2ship30links',bi =birec,cadena =cadrec)


def dibujar(bizq,cadena):
    '''dibuja a partir de datos recogidos en arrays de datos tipo barco y cadena
    ver el sistema de preparar matrices para guardar datos'''
    #figure()
    pl.hold(True)
    for i in range(bizq.shape[0]-1):    
        
        cms = np.array([cadena[i,0,:],cadena[i,1,:]])
        para = np.array([cadena[i,-2,:],cadena[i,-1,:]])
        pl.plot(cms[0,:],cms[1,:],'o')
        pl.hold(True)
        barrasi = cms + cadena[-1,1,0] * para
        barrasd = cms - cadena[-1,1,0] * para
        pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],'k')
        
        pl.plot(bizq[i,0],bizq[i,1],'+b')
        
        vertices = np.array([[-1.,-0.25],[-1.,0.25],[-0.25,0.35],[1,0],\
        [-0.25,-0.35],[-1.,-0.25]])
        rot = np.array([[np.cos(bizq[i,6]),- np.sin(bizq[i,6])],[np.sin(bizq[i,6]),\
        np.cos(bizq[i,6])]])       
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bizq[i,0],bizq[i,1]]       
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
       
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = 'blue')
        pl.gca().add_patch(patchi)
        
        
     
        
dibujar(birec,cadrec)
    
    