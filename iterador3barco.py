# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:08:57 2014
Este programa esta pensado para iterar cadenas y barcos y ver si peta o que
pasa

@author: juan
"""
import barcosolo
import boom
import barchain_matrix
import numpy as np
from matplotlib import pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
delta = 0.01
#creamos los barcos
bi = barcosolo.barco()
bd = barcosolo.barco()
bc = barcosolo.barco()

#creamos la cadena de momento pequeñita porque es para probar
cd = boom.cadena(60)
#situamos la cadena ... toda en el eje x con los el primer y último eslabon.
#apuntando para arriba.
cd.normal[:,0] = [1,0]
cd.normal[:,0] = [1,0]
cd.normal[:,-1] = [-1,0]
cd.para[:,0] = [0,1]
cd.para[:,-1] = [0,-1]
cd.calcms()
cd.A = 0.1
#posicionamos el barco del centro
bc.link = 31

cd.calcms()

#inicializamos los barquitos,poniendolos al principio y al final de la cadena
bi.pb[:] = [0,1.5]
bi.Fm = 500
bi.M = 0
bc.pb[:] = [cd.cms[0,bc.link-1]-0.5,0.5]
bc.Fm = 500
bc.M = 0
bd.pb[:] = [cd.cms[0,-1],1.5]
bd.Fm = 500
bd.M = 0

#iniciamos la matriz de tensiones
M, B, T = barchain_matrix.Matriz_t_ini(cd,1)
#movemos los elabones, estamos usando el paso de integracion por defecto
#delta = 5e-4

tam = 10000 #aqui definimos de momento el tamaño del experimento...
step = 500# y aqui cada cuantas iteraciones guardamos datos y pintamos
birec = np.zeros((tam//step+1,16)) #el barco izquierdo guarda los
#de la tensión con el primer eslabon... (chapuzilla pero en fin) hay
#una tension más que eslabones
bcrec = np.zeros((tam//step+1,16))
bdrec = np.zeros((tam//step+1,14))
cadrec = np.zeros((tam//step+1,14,cd.esl))
rec = np.zeros((tam//step+1,14,cd.esl))





for i in range(tam):
    #M, B, T = barchain_matrix.Matriz_t_ini(cd,1)
    M, B, T = barchain_matrix.Matriz_t(M ,B , T, cd, bi, bd, bc)
    j = bc.link
    #pasamos las tensiones de la cadena a la cadena...
    cd.T[:,0] = np.concatenate([T[:2*j-3:2],T[2*j:-3:2]])
    cd.T[:,1] = np.concatenate([T[1:2*j-2:2],T[2*j+1:-2:2]])
    cd.T[:,2] = np.concatenate([T[2:2*j-1:2],T[2*j+2:-1:2]])
    cd.T[:,3] = np.concatenate([T[3:2*j:2],T[2*j+3::2]])
    cd.movimiento(delta = delta)
    
    #y a continuación movemos los barcos
    bi.movimientotst(extf = T[0:2],dfcm = [-bi.ls/2.,0],delta = delta)
    bc.movimientotst(extf = -T[2*j-2:2*j] +\
    T[2*j:2*j+2],dfcm = [-bc.ls/2.,0], delta = delta)
    bd.movimientotst(extf = - T[-2:],dfcm = [-bd.ls/2.,0],delta = delta)
#    cd1.Matriz_t(bi,bd)
     #famoso parche para corregir errores de redondeo
    popabi = bi.pb - bi.ls/2*np.array([np.cos(bi.theta),np.sin(bi.theta)]) 
    primero = popabi - cd.L*cd.para[:,0]
    popabd = bd.pb - bd.ls/2*np.array([np.cos(bd.theta),np.sin(bd.theta)]) 
    ultimo = popabd + cd.L*cd.para[:,-1]
    popabc =  bc.pb - bc.ls/2*np.array([np.cos(bc.theta),np.sin(bc.theta)])
    
    

    #aprovechamos para evaluar los errores cometidos, antes de recolocar los 
    #eslabones. Esto solo tiene sentido en fase de depuración luego se quita.
#    nudoss[:,0] = popabi
#    nudosa[:,-1] = popabd
#    nudoss[:,1:] =   cd5.cms - cd5.L*cd5.para
#    nudosa[:,:-1] =  cd5.cms + cd5.L*cd5.para
#    errores = nudoss - nudosa
#    grandes = abs(errores) > abs(maxerr)
#    maxerr = grandes * errores + (1-grandes) * maxerr 
        
    #recoloco los eslabones...
    cd.calcms(primero,ultimo,-1*cd.ord)
#    #recoloco los barcos... en realidad solo haría falta recolocar uno...    
    bi.pb = cd.cms[:,0]  + bi.ls/2*np.array([np.cos(bi.theta),np.sin(bi.theta)]) \
    + cd.L*cd.para[:,0] 
    bd.pb = cd.cms[:,-1] + bd.ls/2*np.array([np.cos(bd.theta),np.sin(bd.theta)]) \
    - cd.L*cd.para[:,-1]
    #El siguiente lo recoloco con el eslabon siguiente... Da igual hacerlo
    # con el anterior, puesto que los eslabones ya se han recolocado antes
#    
    bc.pb = cd.cms[:,bc.link-1] + bc.ls/2*np.array([np.cos(bc.theta),np.sin(bc.theta)]) \
    + cd.L*cd.para[:,bc.link-1]
    
#    hold(True) 
    if i%step == 0:
        #guadamomos datos.
        birec[i//step] = bi.pb[0],bi.pb[1],bi.vb[0],bi.vb[1],bi.ab[0],bi.ab[1],\
        bi.theta,bi.wb,bi.alfab,bi.Fm,bi.M,bi.setl,bi.setw,bi.thewj,T[0],\
        T[1]
        
        bcrec[i//step] = bc.pb[0],bc.pb[1],bc.vb[0],bc.vb[1],bc.ab[0],bc.ab[1],\
        bc.theta,bc.wb,bc.alfab,bc.Fm,bc.M,bc.setl,bc.setw,bc.thewj,T[0],\
        T[1]
        
        bdrec[i//step] = bd.pb[0],bd.pb[1],bd.vb[0],bd.vb[1],bd.ab[0],bd.ab[1],\
        bd.theta,bd.wb,bd.alfab,bd.Fm,bd.M,bd.setl,bd.setw,bd.thewj

        cadrec[i//step] = cd.cms[0],cd.cms[1],cd.v[0],cd.v[1],cd.a[0],\
        cd.a[1],cd.alfa,cd.w,T[2:-3:2],T[2:-2:2],cd.normal[0],\
        cd.normal[1],cd.para[0],cd.para[1]
        
        
#en la ultima fila guardamos propidades de los barcos y las cadenas y detalles
#del experimento       
        birec[-1,0:12] = bi.thewjmax,bi.Ac,bi.mb,bi.mul,bi.Ib,bi.mua,bi.ls,bi.mut,\
        bi.pmax,bi.pmaxw,bi.Ab,bi.Aw
        
        bcrec[-1,0:12] = bc.thewjmax,bc.Ac,bc.mb,bc.mul,bc.Ib,bc.mua,bc.ls,bc.mut,\
        bc.pmax,bc.pmaxw,bc.Ab,bc.Aw
        
        bdrec[-1,0:12] = bd.thewjmax,bd.Ac,bd.mb,bd.mul,bd.Ib,bd.mua,bd.ls,bd.mut,\
        bd.pmax,bd.pmaxw,bd.Ab,bd.Aw
        T[:2*j-1:2]
        cadrec[-1,0,0:8] = cd.s,cd.q,cd.A,cd.L,cd.m,cd.I,delta,tam
        
        np.savez('tugofwar2ship30links',bi =birec,bc = bcrec,bd =bdrec,\
        cadena =cadrec)


def dibujar(bizq,bctr,bdcha,cadena1):
    '''dibuja a partir de datos recogidos en arrays de datos tipo barco y cadena
    ver el sistema de preparar matrices para guardar datos'''
    pl.hold(True)
    for i in range(bizq.shape[0]-1):    
 
        cms = np.array([cadena1[i,0,:],cadena1[i,1,:]])
        para = np.array([cadena1[i,-2,:],cadena1[i,-1,:]])
        pl.plot(cms[0,:],cms[1,:],'o')
        barras1i = cms + cadena1[-1,0,3] * para
        barras1d = cms - cadena1[-1,0,3] * para
        pl.plot([barras1i[0,:],barras1d[0,:]],[barras1i[1,:],barras1d[1,:]],'k')
        pl.plot(bizq[i,0],bizq[i,1],'+b')
        pl.plot(bctr[i,0],bctr[i,1],'+k')
        pl.plot(bdcha[i,0],bdcha[i,1],'+r')
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
        
        rot = np.array([[np.cos(bctr[i,6]),- np.sin(bctr[i,6])],[np.sin(bctr[i,6]),\
        np.cos(bctr[i,6])]])       
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bctr[i,0],bctr[i,1]]       
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
       
        pathc = Path(vertrot,codes)
        patchc = patches.PathPatch(pathc,facecolor = 'green')
        pl.gca().add_patch(patchc)
        
        rot = np.array([[np.cos(bdcha[i,6]),- np.sin(bdcha[i,6])],[np.sin(bdcha[i,6]),\
        np.cos(bdcha[i,6])]])       
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bdcha[i,0],bdcha[i,1]]
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
     
        pathd = Path(vertrot,codes)
        patchd = patches.PathPatch(pathd,facecolor = 'red')
        pl.gca().add_patch(patchd)
        pl.pause(0.01)        
        
#        plot(cd1.cms[0,:],cd1.cms[1,:],'o')
#        barrasi = cd1.cms + cd1.L*cd1.para
#        barrasd= cd1.cms - cd1.L*cd1.para   
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
dibujar(birec,bcrec,bdrec,cadrec)
    
    