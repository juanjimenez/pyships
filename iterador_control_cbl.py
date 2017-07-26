# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:08:57 2014
Este programa esta pensado para iterar cadenas y barcos incluyendo
control de rumbo y velocidad y ver si peta o que
pasa
check passed 07.03.2017
@author: juan
"""
#import barcosolo_nv
import usv
import boom
import corrientes
import bezier_cvr as bz

import numpy as np
from numpy.linalg.linalg import norm 
from matplotlib import pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
from time import strftime
from os import mkdir, stat, path

delta = 0.00001
tiempo= 0

#fichero = strftime("%b%d_%Y_%H%M") para guardar resultados
fichero = "../"+strftime("%b%d_%Y").lower() + "/" +\
strftime("%b%d_%Y_%H%M").lower()
directorio = path.dirname(fichero)
try:
    stat(directorio)
except:
    mkdir(directorio)


#creamos los barcos
bi = usv.barco(1) #,'timon')
bd = usv.barco(2) #,'timon')

#creamos la cadena de momento pequeñitca porque es para probar
#se trata de una cadena de 200 eslabones con unos cables de 10cm je, je que
#tienen una constante elastica de 1000...
cd5 = boom.cadena(200,np.array([0.1,0.1,1000]))
#cd5.L = 0.5
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
cd5.q = 0.
cd5.s = 50. #0.
 

#cd5.s2 = 100
#cd5.q2 = 0.01
#cd5.I = cd5.m * cd5.L**2 / 12
#inicializamos los barquitos

#ajustamos los valores de los paramatros de los barcos...
 
bi.setw = 0
bi.ls = 2. #5.
bi.mb = 350.
bi.pmax = 5000.
#bi.mul2 = 100
#bi.mut2 = 1000
bi.M = 0.
bi.rlink = [[bi.ide,bi.pb,bi.theta,bi.vb],[]]

    
#bi.Fm = 1000.
#bd.Fm = 1000.


bd.setw = 0
bd.ls = 2.
bd.mb = 350.
bd.pmax = 5000.
#bd.mul2 = 100
#bd.mut2 = 1000
bd.M = 0.
bd.rlink = [[bd.ide,bd.pb,bd.theta,bd.vb],[]]



#los ponemos al principio y al final de la cadena (a~nadimos a la y la longitud
#del cable...en realidad multiplicada por 0.8 en dirección i para que
# no arranque tenso)

bi.pb = cd5.cms[:,0]  + bi.ls/2.*np.array([np.cos(bi.theta),np.sin(bi.theta)]) \
+ cd5.L*cd5.para[:,0] + [0, 0.8 * cd5.cbl[0]]

bd.pb = cd5.cms[:,-1] + bd.ls/2.*np.array([np.cos(bd.theta),np.sin(bd.theta)]) \
- cd5.L*cd5.para[:,-1] + [0,0.8 * cd5.cbl[1]]


##############################################################################
#probatinas de controladores a ver si ajustamos bien la cosa
###############################################################################




#controles de alineacion longitudinal
bi.krl = [2000.,100.,1000]
bd.krl = [2000.,100.,1000]
#controles de distancia transversal
bi.krd = [1.5,0.01,0.01] #[1,5.,0.1]
bd.krd = [1.5,0.01,0.01] #[1,5.,0.1]

#controles de rumbo standard  
bi.krpid = [5.,10.,0.1] #[5.,10.,0.1]
bd.krpid = [5.,10.,0.1] #[5.,10.,0.1]

#controles de velocidad standard
bi.kvpid = [2000,100,500]
bd.kvpid = [2000,100,500]

###############################################################################

#link-amos las dos radios. Esto exige una solución más elegante
bi.rlink[1] = bd.rlink[0]
bd.rlink[1] = bi.rlink[0]


#inicializamos limites de representación de corrientes.#######################
#lim = [bi.pb[0]-bi.ls,bd.pb[0]+ bd.ls,bi.pb[1]-bi.ls,bd.pb[1]+bd.ls]

###############################################################################
#tamaño experimento
###############################################################################
tam = 20000000 #aqui definimos de momento el tamaño del experimento...
step =200000#0000 #1000
# y aqui cada cuantas iteraciones guardamos datos y pintamos
birec = np.zeros((tam//step+2,16)) #el barco izquierdo guarda los valores
#de la tensión con el primer eslabon... (chapuzilla pero en fin) hay
#una tension más que eslabones
bdrec = np.zeros((tam//step+2,16))
cadrec = np.zeros((tam//step+2,14,cd5.esl))
#cd5.Fd =[0., 10.]poder

###############################################################################
#consignas...
###############################################################################
c_rumbo = np.pi/4
d_sway = 25
d_surge = 25 #-2
###############################################################################
#EMPIEZA BUCLE CALCULO
###############################################################################
for i in range(tam+1):
    #print i, i%step

    #bi.vr = - corrientes.gencor(bi.pb) + bi.vb
    #bd.vr = - corrientes.gencor(bd.pb) + bd.vb
    bi.propulsion(delta,extf = cd5.T[0:2],dfcm = [-bi.ls/2.,0])
    bd.propulsion(delta,extf = - cd5.T[-2:],dfcm = [-bd.ls/2.,0])
        
    
    
    cd5.mtr_s(bi,bd)
    #calculamos la velocidad de arrastre para la cadena
    #cd5.vr = - corrientes.gencor(cd5.cms) + cd5.v    
    cd5.movimiento(delta = delta)

    
    #y a continuación movemos los barcos
    
 #   bi.movimientotst(delta = delta, extf = T[0:2],dfcm = [-bi.ls/2.,0])
    
    bi.movimiento(delta = delta )
    bd.movimiento(delta = delta)
    
        
 #   bd.movimientotst(delta = delta, extf = - T[-2:],dfcm = [-bi.ls/2.,0])
        
    
    #print 'izq'
    bi.controlador(delta,c_rumbo,1.,1000,np.array([d_surge,-d_sway]))#(delta,0.,0.4,500,3)    
    #print 'dcha'    
    bd.controlador(delta,c_rumbo,1.,1000,np.array([-d_surge,d_sway]))#(delta,0.,0.4,500,3)
    

####actualizamos límites para pintar corrientes#################################
#    lim = [min([bi.pb[0]-bi.ls,bd.pb[0]-bd.ls,lim[0]]),\
#    max([bi.pb[0]+bi.ls,bd.pb[0]+bd.ls,lim[1]]),\
#    min([bi.pb[1]-bi.ls,bd.pb[1]-bd.ls,lim[2]]),\
#    max([bi.pb[1]+bi.ls,bd.pb[1]+bd.ls,lim[3]])]
##############################################################################
    
    
       
#    figure(2)
#    #plot(time,bi.vb[0],'.r')
#    #plot(time,bd.vb[0],'.g')
#    plot(time,bi.pb[0],'r+')
#    plot(time,bd.pb[0],'g+')
    
    tiempo= tiempo + delta

#    
    
#==============================================================================
#    #famoso parche para corregir errores de redondeo
#    popabi = bi.pb - bi.ls/2*np.array([np.cos(bi.theta),np.sin(bi.theta)]) 
#    primero = popabi - cd5.L*cd5.para[:,0]
#    popabd = bd.pb - bd.ls/2*np.array([np.cos(bd.theta),np.sin(bd.theta)]) 
#    ultimo = popabd + cd5.L*cd5.para[:,-1]
#    
#    
#
#    #aprovechamos para evaluar los errores cometidos, antes de recolocar los 
#    #eslabones. Esto solo tiene sentido en fase de depuración luego se quita.
##    nudoss[:,0] = popabi
##    nudosa[:,-1] = popabd
##    nudoss[:,1:] =   cd5.cms - cd5.L*cd5.para
##    nudosa[:,:-1] =  cd5.cms + cd5.L*cd5.para
##    errores = nudoss - nudosa
##    grandes = abs(errores) > abs(maxerr)
##    maxerr = grandes * errores + (1-grandes) * maxerr 
##        
##    #recoloco los eslabones...
#    cd5.calcms(primero,ultimo,-1*cd5.ord)
#    #recoloco los barcos... en realidad solo haría falta recolocar uno...    
#    bi.pb = cd5.cms[:,0]  + bi.ls/2*np.array([np.cos(bi.theta),np.sin(bi.theta)]) \
#    + cd5.L*cd5.para[:,0] 
#    bd.pb = cd5.cms[:,-1] + bd.ls/2*np.array([np.cos(bd.theta),np.sin(bd.theta)]) \
#    - cd5.L*cd5.para[:,-1]#       
#==============================================================================

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
#       print i
        #guadamomos datos.
       
        birec[i//step] = bi.pb[0],bi.pb[1],bi.vb[0],bi.vb[1],bi.ab[0],bi.ab[1],\
        bi.theta,bi.wb,bi.alfab,bi.Fm,bi.M,bi.setl,bi.setw,bi.thewj,cd5.T[0],\
        cd5.T[1]

        bdrec[i//step] = bd.pb[0],bd.pb[1],bd.vb[0],bd.vb[1],bd.ab[0],bd.ab[1],\
        bd.theta,bd.wb,bd.alfab,bd.Fm,bd.M,bd.setl,bd.setw,bd.thewj,cd5.T[-2],\
        cd5.T[-1]

        cadrec[i//step] = cd5.cms[0],cd5.cms[1],cd5.v[0],cd5.v[1],cd5.a[0],\
        cd5.a[1],cd5.alfa,cd5.w,cd5.T[2:-1:2],cd5.T[2::2],cd5.normal[0],\
        cd5.normal[1],cd5.para[0],cd5.para[1]
        
 ##############################################################################
 
#en la ultima fila guardamos propidades de los barcos y las cadenas y detalles
#del experimento
        
birec[-1,0:12] = bi.thewjmax,bi.Ac,bi.mb,bi.mul,bi.Ib,bi.mua,bi.ls,bi.mut,\
bi.pmax,bi.pmax,bi.Ab,bi.Aw

bdrec[-1,0:12] = bd.thewjmax,bd.Ac,bd.mb,bd.mul,bd.Ib,bd.mua,bd.ls,bd.mut,\
bd.pmax,bd.pmax,bd.Ab,bd.Aw

cadrec[-1,0,0:3] = cd5.s,cd5.q,cd5.A
cadrec[-1,1,0:3] = cd5.L,cd5.m,cd5.I
cadrec[-1,2,0:2] = delta,tam
cadrec[-1,3,0:3] = cd5.cbl
        
np.savez(fichero,bi =birec,bd =bdrec, cadena =cadrec)


def dibujar(bizq,bdcha,cadena):
    '''dibuja a partir de datos recogidos en arrays de datos tipo barco y cadena
    ver el sistema de preparar matrices para guardar datos'''
    pl.figure(1)
    pl.hold(True)
    for i in range(bizq.shape[0]-1):    
        
        cms = np.array([cadena[i,0,:],cadena[i,1,:]])
        para = np.array([cadena[i,-2,:],cadena[i,-1,:]])
        pl.plot(cms[0,:],cms[1,:],'o')
        pl.hold(True)
        barrasi = cms + cadena[-1,1,0] * para
        barrasd = cms - cadena[-1,1,0] * para
        pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],'k')
        
        pl.plot(bizq[i,0],bizq[i,1],'+r') #r
        pl.plot(bdcha[i,0],bdcha[i,1],'+g') #g
####################Dibujar cable de arrastre################################
        rot = np.array([[np.cos(bizq[i,6]),- np.sin(bizq[i,6])],[np.sin(bizq[i,6]),\
        np.cos(bizq[i,6])]])  
        popai =  np.dot(rot, np.array([-bizq[-1,6]/2.,0])) + [bizq[i,0],bizq[i,1]]  
        tipi = para[:,0] * cadena[-1,1,0] + cms[:,0]
        
        
        disti = norm(popai - tipi)
        di = disti/cadena[-1,3,0]
        print di
        if di > 1: di = 1
        r = bz.bezier4p([[tipi[0]],[tipi[1]]],[[popai[0]],[popai[1]]],1,1,1.5,\
        (1-di) * bizq[i,6]\
         +di * np.arctan2(popai[1]-tipi[1],popai[0] - tipi[0]),\
        (1-di) * np.arctan2(-para[0,0],-para[0,1])\
         +di * np.arctan2(popai[1]-tipi[1],popai[0] - tipi[0]),\
        100)
        bz.pintar_bezier(r[0],color = 'b')
        

        
###############################################################################        
        
#        pl.plot([bizq[i,0],bdcha[i,0]],[bizq[i,1],bdcha[i,1]])
#        normal = bizq[i,1] - bdcha[i,1], - bizq[i,0] + bdcha[i,0]
#        pl.plot([(bizq[i,0]+bdcha[i,0])/2,(bizq[i,0]+bdcha[i,0])/2+normal[0]],\
#        [(bizq[i,1]+bdcha[i,1])/2,(bizq[i,1]+bdcha[i,1])/2 + normal[1]])
#       revisa esto: pinta el barco separado al alargar o acortar la eslora        
#        vertices = np.array([[-1.,-0.25],[-1.,0.25],[-0.25,0.35],[4,0],\
#        [-0.25,-0.35],[-1.,-0.25]])
        
        vertices = np.array([[-bizq[-1,6]/2.,-0.25],[-bizq[-1,6]/2.,0.25],\
        [-0.25,0.35],[bizq[-1,6]/2.,0],[-0.25,-0.35],[-bizq[-1,6]/2.,-0.25]])
        
#        rot = np.array([[np.cos(bizq[i,6]),- np.sin(bizq[i,6])],[np.sin(bizq[i,6]),\
#        np.cos(bizq[i,6])]])       
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bizq[i,0],bizq[i,1]]       
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
       
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = 'red') #'red'
        pl.gca().add_patch(patchi)

#######################dibujar cable de arrastre derecha#######################
        rot = np.array([[np.cos(bdcha[i,6]),- np.sin(bdcha[i,6])],[np.sin(bdcha[i,6]),\
        np.cos(bdcha[i,6])]])  
        popad =  np.dot(rot, np.array([-bdcha[-1,6]/2.,0])) + [bdcha[i,0],bdcha[i,1]]  
        tipd = - para[:,-1] * cadena[-1,1,0] + cms[:,-1]
        
        
        distd = norm(popad - tipd)
        dd = distd/cadena[-1,3,0]
        print di
        if dd > 1: dd = 1
        r = bz.bezier4p([[tipd[0]],[tipd[1]]],[[popad[0]],[popad[1]]],1,1,1.5,\
        (1-dd) * bdcha[i,6]\
         +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
        (1-dd) * np.arctan2(-para[0,0],-para[0,1])\
         +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
        100)
        bz.pintar_bezier(r[0],color = 'b')
###############################################################################

        vertices = np.array([[-bdcha[-1,6]/2.,-0.25],[-bdcha[-1,6]/2.,0.25],\
        [-0.25,0.35],[bdcha[-1,6]/2.,0],[-0.25,-0.35],[-bdcha[-1,6]/2.,-0.25]])        
#        rot = np.array([[np.cos(bdcha[i,6]),- np.sin(bdcha[i,6])],[np.sin(bdcha[i,6]),\
#        np.cos(bdcha[i,6])]])       
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bdcha[i,0],bdcha[i,1]]
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
     
        pathd = Path(vertrot,codes)
        patchd = patches.PathPatch(pathd,facecolor = 'green') #'green'
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
######para pintar  corrientes (incluido el 13.04.2016)#########################
#corrientes.vercor(lim,[20,20],corrientes.camp)
###############################################################################
pl.figure(2)
a = []
for i in zip(bdrec[:-1,2],bdrec[:-1,3]):
    a.append(norm(i))
    
pl.plot(a,'g')
    

a = []
for i in zip(birec[:-1,2],birec[:-1,3]):
    a.append(norm(i))
    
pl.plot(a,'r')
pl.title('modulo de la velocidad')

a = []

for i in zip(bdrec[:-1,-2],bdrec[:-1,-1]):
    a.append(norm(i)) 

pl.figure(3)
pl.plot(a,'g')

a = []    
for i in zip(birec[:-1,-2],birec[:-1,-1]):
    a.append(norm(i)) 

pl.plot(a,'r')

pl.title('Tension en los extremos')
    
pl.figure(4)
pl.plot(np.ones(bdrec.shape[0])*0.0)
pl.plot(bdrec[:-1,6],'g')
pl.plot(birec[:-1,6],'r') 
pl.title('rumbo')   

pl.figure(5)
a = []
b = []
for i in zip(birec[:-1,0],birec[:-1,1],bdrec[:-1,0],bdrec[:-1,1]):
    a.append((i[2]-i[0])*np.cos(c_rumbo)+(i[3]-i[1])*np.sin(c_rumbo))
    b.append((i[2]-i[0])*np.sin(c_rumbo)-(i[3]-i[1])*np.cos(c_rumbo))
#distancia surge
pl.plot(a,'k')
#distancia en sway
pl.plot(b,'r')
pl.title('distancia entre b')   