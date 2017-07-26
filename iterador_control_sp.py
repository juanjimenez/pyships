# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:08:57 2014
Este programa esta pensado para iterar cadenas y barcos incluyendo
control de rumbo y velocidad y ver si peta o que
pasa

03.06.2017. DEPRECATED: This file classes versions which are no longer
in use. See iterador_control_cbl.py for an up to date version of this script
@author: juan
"""
import barcosolo
import boom
import barchain_matrix_sp
import numpy as np
from matplotlib import pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches

delta = 0.01
tiempo= 0 

#creamos los barcos
bi = barcosolo.barco(1)
bd = barcosolo.barco(2)

#creamos la cadena de momento pequeñita porque es para probar
cd5 = boom.cadena(25)

#situamos la cadena ... toda en el eje x con los el primer y último eslabon.
#apuntando para arriba.

cd5.normal[:,0] = [1,0]
cd5.normal[:,0] = [1,0]

cd5.normal[:,-1] = [-1,0]
cd5.para[:,0] = [0,1]
cd5.para[:,-1] = [0,-1]
cd5.calcms()

cd5.m = 2.
cd5.A = 2.
cd5.s = 50.
#inicializamos los barquitos

bi.setl = 0
bi.setw = 0
bi.ls = 2.
bi.mb = 350.
bi.pmax = 5000.
bi.M = 0.
bi.rlink = [[bi.ide,bi.pb,bi.theta,bi.vb],[]]

#los ponemos al principio y al final de la cadena
bi.pb = cd5.cms[:,0]  + bi.ls/2.*np.array([np.cos(bi.theta),np.sin(bi.theta)]) \
+ cd5.L*cd5.para[:,0] 
    
#bi.Fm = 499.50
#bd.Fm = 499.50


bd.setl = 0
bd.setw = 0
bd.ls = 2.
bd.mb = 350.
bd.pmax = 5000.
bd.M = 0.
bd.rlink = [[bd.ide,bd.pb,bd.theta,bd.vb],[]]

bd.pb = cd5.cms[:,-1] + bd.ls/2.*np.array([np.cos(bd.theta),np.sin(bd.theta)]) \
- cd5.L*cd5.para[:,-1]
#probatinas de controladores a ver si ajustamos bien la cosa

bi.krali = [10.,0.,0.]
bi.krdis = [0.5,5.,0.01]
bi.krpid = [10.,10.,0.1]
bi.cfr = 500.
#
bi.kvali = [1500, 0, 0]
#bi.kvpid
#
bd.krali = [10.,0.,0.]
bd.krdis = [0.5,5.,0.01]
bd.krpid = [10.,10.,0.1]
bi.cfr = 500.
#
bd.kvali = [1500., 0., 0.]
#bd.kvpid

#link-amos las dos radios. Esto exige una solución más elegante
bi.rlink[1] = bd.rlink[0]
bd.rlink[1] = bi.rlink[0]

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
tam = 50 #50000 #aqui definimos de momento el tamaño del experimento...
step = 2 #1000
# y aqui cada cuantas iteraciones guardamos datos y pintamos
birec = np.zeros((tam//step+1,16)) #el barco izquierdo guarda los valores
#de la tensión con el primer eslabon... (chapuzilla pero en fin) hay
#una tension más que eslabones
bdrec = np.zeros((tam//step+1,14))
cadrec = np.zeros((tam//step+1,14,cd5.esl))
#cd5.Fd =[0., 10.]poder
M, B, T = barchain_matrix_sp.Matriz_t_ini(cd5,0)

for i in range(tam+1):
    #print i, i%step

    
    bi.propulsion(delta,extf = T[0:2],dfcm = [-bi.ls/2.,0])
    bd.propulsion(delta,extf = - T[-2:],dfcm = [-bi.ls/2.,0])
    
    M, B, T = barchain_matrix_sp.Matriz_t_ini(cd5,0)
    M, B, T = barchain_matrix_sp.Matriz_t(M, B, T, cd5,bi,bd)
         
    cd5.T[:,0] = T[:-3:2]
    cd5.T[:,1] = T[1:-2:2]
    cd5.T[:,2] = T[2:-1:2]
    cd5.T[:,3] = T[3::2]
    cd5.movimiento(delta = delta)

    
    #y a continuación movemos los barcos
    
 #   bi.movimientotst(delta = delta, extf = T[0:2],dfcm = [-bi.ls/2.,0])
    
    bi.movimiento(delta = delta )
 #   bd.rlink[1] = bi.rlink[0]
    bi.controlador(delta,pi/2.,1.,1000,5)
    
        
 #   bd.movimientotst(delta = delta, extf = - T[-2:],dfcm = [-bi.ls/2.,0])
        
    bd.movimiento(delta = delta)
    bd.controlador(delta,pi/2.,1.,1000,5)    
     
#    figure(2)
#    #plot(time,bi.vb[0],'.r')
#    #plot(time,bd.vb[0],'.g')
#    plot(time,bi.pb[0],'r+')
#    plot(time,bd.pb[0],'g+')
    
    tiempo= tiempo+ delta

#    
    
      
    #famoso parche para corregir errores de redondeo
#    popabi = bi.pb - bi.ls/2*np.array([np.cos(bi.theta),np.sin(bi.theta)]) 
#    primero = popabi - cd5.L*cd5.para[:,0]
#    popabd = bd.pb - bd.ls/2*np.array([np.cos(bd.theta),np.sin(bd.theta)]) 
#    ultimo = popabd + cd5.L*cd5.para[:,-1]
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
##        
##    #recoloco los eslabones...
#    cd5.calcms(primero,ultimo,-1*cd5.ord)
#    #recoloco los barcos... en realidad solo haría falta recolocar uno...    
#    bi.pb = cd5.cms[:,0]  + bi.ls/2*np.array([np.cos(bi.theta),np.sin(bi.theta)]) \
#    + cd5.L*cd5.para[:,0] 
#    bd.pb = cd5.cms[:,-1] + bd.ls/2*np.array([np.cos(bd.theta),np.sin(bd.theta)]) \
#    - cd5.L*cd5.para[:,-1]
##    print 'fuera izquierdo', bi.pb
##    print 'fuera derecho', bd.pb
    
    
    if (i%step == 0):
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

        bdrec[i//step] = bd.pb[0],bd.pb[1],bd.vb[0],bd.vb[1],bd.ab[0],bd.ab[1],\
        bd.theta,bd.wb,bd.alfab,bd.Fm,bd.M,bd.setl,bd.setw,bd.thewj

        cadrec[i//step] = cd5.cms[0],cd5.cms[1],cd5.v[0],cd5.v[1],cd5.a[0],\
        cd5.a[1],cd5.alfa,cd5.w,T[2:-1:2],T[2::2],cd5.normal[0],\
        cd5.normal[1],cd5.para[0],cd5.para[1]
#en la ultima fila guardamos propidades de los barcos y las cadenas y detalles
#del experimento       
        birec[-1,0:12] = bi.thewjmax,bi.Ac,bi.mb,bi.mul,bi.Ib,bi.mua,bi.ls,bi.mut,\
        bi.pmax,bi.pmax,bi.Ab,bi.Aw
        
        bdrec[-1,0:12] = bd.thewjmax,bd.Ac,bd.mb,bd.mul,bd.Ib,bd.mua,bd.ls,bd.mut,\
        bd.pmax,bd.pmax,bd.Ab,bd.Aw
        
        cadrec[-1,0,0:3] = cd5.s,cd5.q,cd5.A
        cadrec[-1,1,0:3] = cd5.L,cd5.m,cd5.I
        cadrec[-1,2,0:2] = delta,tam
np.savez('tugofwar2ship30links',bi =birec,bd =bdrec, cadena =cadrec)


def dibujar(bizq,bdcha,cadena):
    '''dibuja a partir de datos recogidos en arrays de datos tipo barco y cadena
    ver el sistema de preparar matrices para guardar datos'''
    figure(1)
    pl.hold(True)
    for i in range(bizq.shape[0]-1):    
        
        cms = np.array([cadena[i,0,:],cadena[i,1,:]])
        para = np.array([cadena[i,-2,:],cadena[i,-1,:]])
        pl.plot(cms[0,:],cms[1,:],'o')
        pl.hold(True)
        barrasi = cms + cadena[-1,1,0] * para
        barrasd = cms - cadena[-1,1,0] * para
        pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],'k')
        
        pl.plot(bizq[i,0],bizq[i,1],'+r')
        pl.plot(bdcha[i,0],bdcha[i,1],'+g')
#        pl.plot([bizq[i,0],bdcha[i,0]],[bizq[i,1],bdcha[i,1]])
#        normal = bizq[i,1] - bdcha[i,1], - bizq[i,0] + bdcha[i,0]
#        pl.plot([(bizq[i,0]+bdcha[i,0])/2,(bizq[i,0]+bdcha[i,0])/2+normal[0]],\
#        [(bizq[i,1]+bdcha[i,1])/2,(bizq[i,1]+bdcha[i,1])/2 + normal[1]])
#       revisa esto: pinta el barco separado al alargar o acortar la eslora        
#        vertices = np.array([[-1.,-0.25],[-1.,0.25],[-0.25,0.35],[4,0],\
#        [-0.25,-0.35],[-1.,-0.25]])
        
        vertices = np.array([[-bizq[-1,6]/2.,-0.25],[-bizq[-1,6]/2.,0.25],\
        [-0.25,0.35],[bizq[-1,6]/2.,0],[-0.25,-0.35],[-bizq[-1,6]/2.,-0.25]])
        
        rot = np.array([[np.cos(bizq[i,6]),- np.sin(bizq[i,6])],[np.sin(bizq[i,6]),\
        np.cos(bizq[i,6])]])       
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bizq[i,0],bizq[i,1]]       
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
       
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = 'red')
        pl.gca().add_patch(patchi)

        vertices = np.array([[-bdcha[-1,6]/2.,-0.25],[-bdcha[-1,6]/2.,0.25],\
        [-0.25,0.35],[bdcha[-1,6]/2.,0],[-0.25,-0.35],[-bdcha[-1,6]/2.,-0.25]])        
        rot = np.array([[np.cos(bdcha[i,6]),- np.sin(bdcha[i,6])],[np.sin(bdcha[i,6]),\
        np.cos(bdcha[i,6])]])       
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bdcha[i,0],bdcha[i,1]]
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
     
        pathd = Path(vertrot,codes)
        patchd = patches.PathPatch(pathd,facecolor = 'green')
        pl.gca().add_patch(patchd)
        
        
        pl.plot(bizq[i,0],bizq[i,1],'+r')
        pl.plot(bdcha[i,0],bdcha[i,1],'+g')
        
        pl.pause(0.01)        
        #hold(False) 
        
#        plot(cd5.cms[0,:],cd5.cms[1,:],'o')
#        barrasi = cd5.cms + cd5.L*cd5.para
#        barrasd= cd5.cms - cd5.L*cd5.para   
#        plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],'k')
#        plot(bd.pb[0],bd.pb[1],'+b')  
#        #plot([bd.pb[0]-np.cos(bd.theta),bd.pb[0]+np.cos(bd.theta)],\
#        #[bd.pb[1]-np.sin(bd.theta),bd.pb[1]+np.sin(bd.theta)])
#        plot(bi.pb[0],bi.pb[1],'+r')
#        #plot([bi.pb[0]-np.cos(bi.theta),bi.pb[0]+np.cos(bi.theta)],\
#        #[bi.pb[1]-np.sin(bi.theta),bi.pb[1]+np.sin(bi.theta)])
#        
#        vertices = array([[-1.,-0.25],[-1.,0.25],[-0.25,0.35],[1,0],\
#        [-0.25,-0.35],[-1.,-0.25]])
#        
#        rot = array([[cos(bi.theta),- sin(bi.theta)],[sin(bi.theta),\
#        cos(bi.theta)]])       
#        vertrot = array([dot(rot,i) for i in vertices]) + bi.pb
#        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
#        Path.CURVE3]
#        #xs, ys = zip(*vertrot)        
#        pathi = Path(vertrot,codes)
#        patchi = patches.PathPatch(pathi,facecolor = 'blue')
#        gca().add_patch(patchi)
#        
#        rot = array([[cos(bd.theta),- sin(bd.theta)],[sin(bd.theta),\
#        cos(bd.theta)]])       
#        vertrot = array([dot(rot,i) for i in vertices]) + bd.pb
#        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
#        Path.CURVE3]
#        #xs, ys = zip(*vertrot)        
#        pathd = Path(vertrot,codes)
#        patchd = patches.PathPatch(pathd,facecolor = 'red')
#        gca().add_patch(patchd)
#        

dibujar(birec,bdrec,cadrec)

figure(2)
a = []
for i in zip(bdrec[:-1,2],bdrec[:-1,3]):
    a.append(norm(i))
    
plot(a,'g')    

a = []
for i in zip(birec[:-1,2],birec[:-1,3]):
    a.append(norm(i))
    
plot(a,'r')

figure(3)
plot(ones(bdrec.shape[0])*0.0)
plot(bdrec[:-1,6],'g')
plot(birec[:-1,6],'r')    

figure(4)
a = []
for i in zip(birec[:-1,0],birec[:-1,1],bdrec[:-1,0],bdrec[:-1,1]):
    a.append(norm([i[2]-i[0],i[3]-i[1]]))
    
plot(a,'k')
   