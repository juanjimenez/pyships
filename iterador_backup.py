# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:08:57 2014
Este programa esta pensado para iterar cadenas y barcos y ver si peta o que
pasa

@author: juan
"""
import barcosolo
import barchain_matrix
from matplotlib.path import Path
import matplotlib.patches as patches
delta = 0.1
#creamos los barcos
bi = barcosolo.barco()
bd = barcosolo.barco()
bc = barcosolo.barco()

#creamos la cadena de momento pequeñita porque es para probar
cd1 = barchain_matrix.cadena()
cd2 = barchain_matrix.cadena()
#situamos la cadena ... toda en el eje x con los el primer y último eslabon.
#apuntando para arriba.

cd1.normal[:,0] = [1,0]
cd1.normal[:,0] = [1,0]
cd1.normal[:,-1] = [-1,0]
cd1.para[:,0] = [0,1]
cd1.para[:,-1] = [0,-1]
cd1.calcms()

cd2.normal[:,0] = [1,0]
cd2.normal[:,0] = [1,0]
cd2.normal[:,-1] = [-1,0]
cd2.para[:,0] = [0,1]
cd2.para[:,-1] = [0,-1]
cd2.calcms(cd1.cms[:,-1])
#inicializamos los barquitos,poniendolos al principio y al final de la cadena
bi.pb[:] = [0,1.5]
bi.Fm = 100
bi.M = 0
bc.pb[:] = [cd1.cms[0,-1],1.5]
bc.Fm = 100
bd.M = 0
bd.pb[:] = [cd2.cms[0,-1],1.5]
bd.Fm = 100
bd.M = 0
#inicializamos la matriz de coeficientes de las tensiones en los extremos
#de los eslabones de las cadenas.

#cd1.Matriz_t(bi,bd)
#movemos los elabones, estamos usando el paso de integracion por defecto
#delta = 5e-4

tam = 1000 #aqui definimos de momento el tamaño del experimento...
step = 100# y aqui cada cuantas iteraciones guardamos datos y pintamos
birec = zeros((tam//step+1,16)) #el barco izquierdo guarda los
#de la tensión con el primer eslabon... (chapuzilla pero en fin) hay
#una tension más que eslabones
bcrec = zeros((tam//step+1,16))
bdrec = zeros((tam//step+1,14))
cad1rec = zeros((tam//step+1,14,cd1.esl))
rec = zeros((tam//step+1,14,cd1.esl))
for i in range(tam):

    cd1.Matriz_t(bi,bc)
    cd2.Matriz_t(bc,bd)
    cd1.movimiento(delta = delta)
    cd2.movimiento(delta = delta)   
    #y a continuación movemos los barcos
    bi.movimientotst(extf = cd1.T[0:2],dfcm = [-bi.ls/2.,0],delta = delta)
    bc.movimientotst(extf = -cd1.T[-2:] + cd2.T[0:2],dfcm = [-bc.ls/2.,0],\
    delta = delta)
    bd.movimientotst(extf = - cd2.T[-2:],dfcm = [-bd.ls/2.,0],delta = delta)
#    cd1.Matriz_t(bi,bd)
    
#    hold(True) 
    if i%step == 0:
        #guadamomos datos.
        birec[i//step] = bi.pb[0],bi.pb[1],bi.vb[0],bi.vb[1],bi.ab[0],bi.ab[1],\
        bi.theta,bi.wb,bi.alfab,bi.Fm,bi.M,bi.setl,bi.setw,bi.thewj,cd1.T[0],\
        cd1.T[1]
        
        bcrec[i//step] = bc.pb[0],bc.pb[1],bc.vb[0],bc.vb[1],bc.ab[0],bc.ab[1],\
        bc.theta,bc.wb,bc.alfab,bc.Fm,bc.M,bc.setl,bc.setw,bc.thewj,cd2.T[0],\
        cd2.T[1]
        
        bdrec[i//step] = bd.pb[0],bd.pb[1],bd.vb[0],bd.vb[1],bd.ab[0],bd.ab[1],\
        bd.theta,bd.wb,bd.alfab,bd.Fm,bd.M,bd.setl,bd.setw,bd.thewj

        cad1rec[i//step] = cd1.cms[0],cd1.cms[1],cd1.v[0],cd1.v[1],cd1.a[0],\
        cd1.a[1],cd1.alfa,cd1.w,cd1.T[2:-1:2],cd1.T[2::2],cd1.normal[0],\
        cd1.normal[1],cd1.para[0],cd1.para[1]
        
        cad2rec[i//step] = cd2.cms[0],cd2.cms[1],cd2.v[0],cd2.v[1],cd2.a[0],\
        cd2.a[1],cd2.alfa,cd2.w,cd2.T[2:-1:2],cd2.T[2::2],cd2.normal[0],\
        cd2.normal[1],cd2.para[0],cd2.para[1]
#en la ultima fila guardamos propidades de los barcos y las cadenas y detalles
#del experimento       
        birec[-1,0:12] = bi.thewjmax,bi.Ac,bi.mb,bi.mul,bi.Ib,bi.mua,bi.ls,bi.mut,\
        bi.pmax,bi.pmaxw,bi.Ab,bi.Aw
        
        bcrec[-1,0:12] = bc.thewjmax,bc.Ac,bc.mb,bc.mul,bc.Ib,bc.mua,bc.ls,bc.mut,\
        bc.pmax,bc.pmaxw,bc.Ab,bc.Aw
        
        bdrec[-1,0:12] = bd.thewjmax,bd.Ac,bd.mb,bd.mul,bd.Ib,bd.mua,bd.ls,bd.mut,\
        bd.pmax,bd.pmaxw,bd.Ab,bd.Aw
        
        cad1rec[-1,0,0:8] = cd1.s,cd1.q,cd1.A,cd1.L,cd1.m,cd1.I,delta,tam
        
        cad2rec[-1,0,0:8] = cd2.s,cd2.q,cd2.A,cd2.L,cd2.m,cd2.I,delta,tam
savez('tugofwar2ship30links',bi =birec,bc = bcrec,bd =bdrec,\
cadena1 =cad1rec, cadena2 = cad2rec)


def dibujar(bizq,bctr,bdcha,cadena1,cadena2):
    '''dibuja a partir de datos recogidos en arrays de datos tipo barco y cadena
    ver el sistema de preparar matrices para guardar datos'''
    hold(True)
    for i in range(bizq.shape[0]-1):    
 
        cms = array([cadena1[i,0,:],cadena1[i,1,:]])
        para = array([cadena1[i,-2,:],cadena1[i,-1,:]])
        plot(cms[0,:],cms[1,:],'o')
        barras1i = cms + cadena1[-1,0,3] * para
        barras1d = cms - cadena1[-1,0,3] * para
        plot([barras1i[0,:],barras1d[0,:]],[barras1i[1,:],barras1d[1,:]],'k')
        
        cms = array([cadena2[i,0,:],cadena2[i,1,:]])
        para = array([cadena2[i,-2,:],cadena2[i,-1,:]])
        plot(cms[0,:],cms[1,:],'o')
        barras2i = cms + cadena2[-1,0,3] * para
        barras2d = cms - cadena2[-1,0,3] * para
        plot([barras2i[0,:],barras2d[0,:]],[barras2i[1,:],barras2d[1,:]],'k')
        
        plot(bizq[i,0],bizq[i,1],'+b')
        plot(bctr[i,0],bctr[i,1],'+k')
        plot(bdcha[i,0],bdcha[i,1],'+r')
        vertices = array([[-1.,-0.25],[-1.,0.25],[-0.25,0.35],[1,0],\
        [-0.25,-0.35],[-1.,-0.25]])
        rot = array([[cos(bizq[i,6]),- sin(bizq[i,6])],[sin(bizq[i,6]),\
        cos(bizq[i,6])]])       
        vertrot = array([dot(rot,j) for j in vertices]) + [bizq[i,0],bizq[i,1]]       
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
       
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = 'blue')
        gca().add_patch(patchi)
        
        rot = array([[cos(bctr[i,6]),- sin(bctr[i,6])],[sin(bctr[i,6]),\
        cos(bctr[i,6])]])       
        vertrot = array([dot(rot,j) for j in vertices]) + [bctr[i,0],bctr[i,1]]       
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
       
        pathc = Path(vertrot,codes)
        patchc = patches.PathPatch(pathc,facecolor = 'green')
        gca().add_patch(patchc)
        
        rot = array([[cos(bdcha[i,6]),- sin(bdcha[i,6])],[sin(bdcha[i,6]),\
        cos(bdcha[i,6])]])       
        vertrot = array([dot(rot,j) for j in vertices]) + [bdcha[i,0],bdcha[i,1]]
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
     
        pathd = Path(vertrot,codes)
        patchd = patches.PathPatch(pathd,facecolor = 'red')
        gca().add_patch(patchd)        
        
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
dibujar(birec,bcrec,bdrec,cad1rec,cad2rec)
    
    