# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:08:57 2014
Este programa esta pensado para iterar cadenas y barcos incluyendo
control de rumbo y velocidad y ver si peta o que
pasa

@author: juan
"""
import barcosolo_nv
import boom_nv
import barchain_matrix_nv_sp
import corrientes
import bezier_cvr


import numpy as np
from numpy.linalg.linalg import norm 
from matplotlib import pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
from time import strftime

delta = 0.001
tiempo= 0

fichero = strftime("%b%d_%Y_%H%M") 

#creamos los barcos
#bi = barcosolo_nv.barco(1) #,'timon')
bd = barcosolo_nv.barco(2) #,'timon')

#creamos la cadena de momento pequeñita porque es para probar
cd5 = boom_nv.cadena(14)
cd5.L = 0.5
#situamos la cadena ... toda en el eje x con los el primer y último eslabon.
#apuntando para arriba.

#cd5.normal[:,0] = [1,0]
cd5.normal[:,-1] = [-1,0]
#cd5.para[:,0] = [0,1]
cd5.para[:,-1] = [0,-1]
cd5.calcms()
cd5.Fi[0] = -1 #condicion de extremo izquierdo fijo
cd5.m = 2.
cd5.A = 2.
cd5.q = 0.
cd5.s = 0.
cd5.s2 = 100
cd5.q2 = 0.01
cd5.I = cd5.m * cd5.L**2 / 12
#inicializamos los barquitos

#ajustamos los valores de los paramatros de los barcos...
 
#bi.setl = 500 #valor de velocidad feedforward
#bi.setw = 0
#bi.ls = 5.
#bi.mb = 350.
#bi.pmax = 5000.
#bi.mul2 = 100
#bi.mut2 = 1000
#bi.M = 0.
#bi.rlink = [[bi.ide,bi.pb,bi.theta,bi.vb],[]]

    
#bi.Fm = 1000.
#bd.Fm = 1000.


bd.setl = 500#1000 #valor de velocidad feedforward
bd.setw = 0
bd.ls = 5.
bd.mb = 350.
bd.pmax = 5000.
bd.mul2 = 100
bd.mut2 = 1000
bd.M = 0.
#bd.rlink = [[bd.ide,bd.pb,bd.theta,bd.vb],[]]



#ponemos el barco al final de la cadena


bd.pb = cd5.cms[:,-1] + bd.ls/2.*np.array([np.cos(bd.theta),np.sin(bd.theta)]) \
- cd5.L*cd5.para[:,-1]
#probatinas de controladores a ver si ajustamos bien la cosa

#bi.krali = [10.,0.,0.]
#bi.krdis = [1,5.,0.1]
#bi.krpid = [5.,10.,0.1]
#bi.kvpid = [30,2,10]
#bi.cfr = 500.
#
#bi.kvali = [15, 0, 0]
#bi.kvpid
#
#bd.krali = [10.,0.,0.]
#bd.krdis = [1,5.,0.1]
bd.krpid = [5.,10.,0.1]
bd.kvpid = [30,2,10]
bd.cfr = 500.
#
#bd.kvali = [15., 0., 0.]
#bd.kvpid

#link-amos las dos radios. Esto exige una solución más elegante
#bi.rlink[1] = bd.rlink[0]
#bd.rlink[1] = bi.rlink[0]

#==============================================================================
# definimos un plan para rodear el buque....

#posicion de buque fondeado
buque = [10.,2.,[-6.,1.],-np.pi] 

#puntos de paso del barco que arastra la cadena
vertices = np.array([[-buque[0]/2.,-buque[1]/2.],\
[-buque[0]/2.,buque[1]/2],\
[buque[0]/4.,buque[1]/2],\
[buque[0]/2,0],\
[buque[0]/4.,-buque[1]/2]])        
rot = np.array([[np.cos(buque[3]),- np.sin(buque[3])],[np.sin(buque[3]),\
np.cos(buque[3])]])

vertrot = np.array([np.dot(rot,j) for j in vertices]) \
+ [buque[2][0],buque[2][1]]

xmax = np.max(vertrot[:,0])
xmin = np.min(vertrot[:,0])
ymax = np.max(vertrot[:,1])
ymin = np.min(vertrot[:,1])
puntos = np.array([[[bd.pb[0]],[bd.pb[1]]],\
[[bd.pb[0]+bd.ls/2.*np.cos(bd.theta)],[bd.pb[1]+bd.ls/2.*np.sin(bd.theta)]],\
[[xmax+1.],[ymax+1.]],\
[[xmin-1.],[ymax+1.]],[[xmin-1.],[ymin+0.5]]])
rumbos = np.array([bd.theta,bd.theta,-np.pi,-np.pi/2.,-np.pi])
plan = bd.planificador(puntos,rumbos,ndat = 20)

#==============================================================================

#inicializamos limites de representación de corrientes.#######################
#lim = [bi.pb[0]-bi.ls,bd.pb[0]+ bd.ls,bi.pb[1]-bi.ls,bd.pb[1]+bd.ls]

##############################################################################



#nudoss = np.zeros([2,cd5.esl+1]) #posición del extremo derecho de cada eslabón

#nudosa = np.zeros([2,cd5.esl+1]) #posisción del extremo izquierdo de cada eslabón
#inicializamos la matriz de coeficientes de las tensiones en los extremos
#de los eslabones de las cadenas.
#errores = np.zeros([2,cd5.esl+1])
#maxerr = np.zeros([2,cd5.esl+1])
#cd5.Matriz_t(bi,bd)
#movemos los elabones, estamos usando el paso de integracion por defecto
#delta = 5e-4
#barras = np.zeros((2,cd5.cms.shape[1]+1))


tam = 50000 #50000 #aqui definimos de momento el tamaño del experimento...
step = 1000 #1000
# y aqui cada cuantas iteraciones guardamos datos y pintamos
birec = np.zeros((tam//step+1,16)) #el barco izquierdo guarda los valores
#de la tensión con el primer eslabon... (chapuzilla pero en fin) hay
#una tension más que eslabones
bdrec = np.zeros((tam//step+1,16))
cadrec = np.zeros((tam//step+1,14,cd5.esl))
#cd5.Fd =[0., 10.]poder
M, B, T = barchain_matrix_nv_sp.Matriz_t_ini(cd5,0)

for i in range(tam+1):
    #print i, i%step

#    bi.vr = - corrientes.gencor(bi.pb) + bi.vb
    bd.vr = - corrientes.gencor(bd.pb) + bd.vb
#    bi.propulsion(delta,extf = T[0:2],dfcm = [-bi.ls/2.,0])
    bd.propulsion(delta,extf = - T[-2:],dfcm = [-bd.ls/2.,0])
        
    
    
    M, B, T = barchain_matrix_nv_sp.Matriz_t_ini(cd5,0)
    M, B, T = barchain_matrix_nv_sp.Matriz_t(M, B, T, cd5,barcod = bd)
         
#    T[0] = T[2]
#    T[1] = T[3]
    T = np.concatenate((np.array([T[0],T[1]]),T))  
    cd5.T[:,0] = T[:-3:2]
    cd5.T[:,1] = T[1:-2:2]
    cd5.T[:,2] = T[2:-1:2]
    cd5.T[:,3] = T[3::2]
    
    
    #calculamos la velocidad de arrastre para la cadena
    #cd5.vr = - corrientes.gencor(cd5.cms) + cd5.v    
    cd5.movimiento(delta = delta)

    
    #y a continuación movemos los barcos
    
 #   bi.movimientotst(delta = delta, extf = T[0:2],dfcm = [-bi.ls/2.,0])
    
#    bi.movimiento(delta = delta )
#    bi.controlador(delta,0.,0.4,500,3)
    
        
 #   bd.movimientotst(delta = delta, extf = - T[-2:],dfcm = [-bi.ls/2.,0])
        
    bd.movimiento(delta = delta)
    bd.controlador(delta,pi/2.,0.4,500,3)    

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
         
        #guadamomos datos.
        birec[i//step] = 0,0,0,0,0,0,\
        0,0,0,0,0,0,0,0,T[0],\
        T[1]

        bdrec[i//step] = bd.pb[0],bd.pb[1],bd.vb[0],bd.vb[1],bd.ab[0],bd.ab[1],\
        bd.theta,bd.wb,bd.alfab,bd.Fm,bd.M,bd.setl,bd.setw,bd.thewj,T[-2],T[-1]

        cadrec[i//step] = cd5.cms[0],cd5.cms[1],cd5.v[0],cd5.v[1],cd5.a[0],\
        cd5.a[1],cd5.alfa,cd5.w,T[2:-1:2],T[3::2],cd5.normal[0],\
        cd5.normal[1],cd5.para[0],cd5.para[1]
#en la ultima fila guardamos propidades de los barcos y las cadenas y detalles
#del experimento       
        birec[-1,0:12] = 0,0,0,0,0,0,0,0,\
        0,0,0,0
        
        bdrec[-1,0:12] = bd.thewjmax,bd.Ac,bd.mb,bd.mul,bd.Ib,bd.mua,bd.ls,bd.mut,\
        bd.pmax,bd.pmax,bd.Ab,bd.Aw
        
        cadrec[-1,0,0:3] = cd5.s,cd5.q,cd5.A
        cadrec[-1,1,0:3] = cd5.L,cd5.m,cd5.I
        cadrec[-1,2,0:2] = delta,tam
        
np.savez(fichero,bi =birec,bd =bdrec, cadena =cadrec)


def dibujar(bdcha,cadena,buque = [10.,2.,[0.,0.],0.]):
    '''dibuja a partir de datos recogidos en arrays de datos tipo barco y cadena
    ver el sistema de preparar matrices para guardar datos dibuja también
    un buque fondeado buque = [eslora, manga, posición, orientación,de     sup'''
    pl.figure(1)
    pl.hold(True)
    vertices = np.array([[-buque[0]/2.,-buque[1]/2.],\
    [-buque[0]/2.,buque[1]/2],\
    [buque[0]/4.,buque[1]/2],\
    [buque[0]/2,0],\
    [buque[0]/4.,-buque[1]/2],\
    [-buque[0]/2.,-buque[1]/2.]])        
    rot = np.array([[np.cos(buque[3]),- np.sin(buque[3])],[np.sin(buque[3]),\
    np.cos(buque[3])]])
    
    vertrot = np.array([np.dot(rot,j) for j in vertices]) \
    + [buque[2][0],buque[2][1]]
                       
    codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
    Path.CURVE3]
   
    pathi = Path(vertrot,codes)
    patchi = patches.PathPatch(pathi,facecolor = 'black') #'red'
    pl.gca().add_patch(patchi)

    for i in range(bdcha.shape[0]-1):    
        
        cms = np.array([cadena[i,0,:],cadena[i,1,:]])
        para = np.array([cadena[i,-2,:],cadena[i,-1,:]])
        pl.plot(cms[0,:],cms[1,:],'yo')
        pl.hold(True)
        barrasi = cms + cadena[-1,1,0] * para
        barrasd = cms - cadena[-1,1,0] * para
        pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],'y')
        
#        pl.plot(bizq[i,0],bizq[i,1],'+r') #r
#        pl.plot(bdcha[i,0],bdcha[i,1],'+g') #g
#        pl.plot([bizq[i,0],bdcha[i,0]],[bizq[i,1],bdcha[i,1]])
#        normal = bizq[i,1] - bdcha[i,1], - bizq[i,0] + bdcha[i,0]
#        pl.plot([(bizq[i,0]+bdcha[i,0])/2,(bizq[i,0]+bdcha[i,0])/2+normal[0]],\
#        [(bizq[i,1]+bdcha[i,1])/2,(bizq[i,1]+bdcha[i,1])/2 + normal[1]])
#       revisa esto: pinta el barco separado al alargar o acortar la eslora        
#        vertices = np.array([[-1.,-0.25],[-1.,0.25],[-0.25,0.35],[4,0],\
#        [-0.25,-0.35],[-1.,-0.25]])
        
        
#       pintamos el barco que supuestamente queremos envolver
        
        

        vertices = np.array([[-bdcha[-1,6]/2.,-0.25],[-bdcha[-1,6]/2.,0.25],\
        [-0.25,0.35],[bdcha[-1,6]/2.,0],[-0.25,-0.35],[-bdcha[-1,6]/2.,-0.25]])        
        rot = np.array([[np.cos(bdcha[i,6]),- np.sin(bdcha[i,6])],[np.sin(bdcha[i,6]),\
        np.cos(bdcha[i,6])]])       
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bdcha[i,0],bdcha[i,1]]
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
     
        pathd = Path(vertrot,codes)
        patchd = patches.PathPatch(pathd,facecolor = 'green') #'green'
        pl.gca().add_patch(patchd)
        
        
#        pl.plot(bizq[i,0],bizq[i,1],'+r')
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

dibujar(bdrec,cadrec,buque)
for i in plan:
    bezier_cvr.pintar_bezier(i[0])
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
for i in zip(birec[:-1,0],birec[:-1,1],bdrec[:-1,0],bdrec[:-1,1]):
    a.append(norm([i[2]-i[0],i[3]-i[1]]))
    
pl.plot(a,'k')
pl.title('distancia entre b')   