# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:08:57 2014

Review: 20.04.2017. 
Simulation for a long boom to enclose a 100m lenght ship
The planing has been defined using Dubins' trayectory 
@author: juan
"""
import usv
import boom
import corrientes
import dubing
import bezier_cvr

import numpy as np
from numpy.linalg.linalg import norm 
from matplotlib import pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
from time import strftime
from os import mkdir, stat, path
#==============================================================================
#    #famoso parche para corregir errores de redondeo(convertirlo en funcion) 
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

delta = 0.00001 #0.000001
t= 0

#fichero = strftime("%b%d_%Y_%H%M") para guardar resultados
fichero = "../"+strftime("%b%d_%Y").lower() + "/" +\
strftime("%b%d_%Y_%H%M").lower()
directorio = path.dirname(fichero)
try:
    stat(directorio)
except:
    mkdir(directorio)

# A usv instance 
#bi = barcosolo_nv.barco(1) #,'timon')
bd = usv.barco(2) #,'timon')

# A boom instance
cd5 = boom.cadena(200,np.array([0,10,1000]),np.array([1,0]))
cd5.L = 0.5
#the boom is collocated folded 

# the links are collocated parallel to the dock
up = 1
for i in range(cd5.esl):
    cd5.normal[:,i] = [0,up]
    cd5.para[:,i] = [up,0]
    up = - up
#cd5.normal[:,0] = [1,0]
cd5.normal[:,-1] = [-1,0]
#cd5.para[:,0] = [0,1]
cd5.para[:,-1] = [0,-1]

cd5.calcms()
cd5.Fi[0] = -1 #the left boom  tip is attached to the dock
cd5.m = 2.
cd5.A = 2.
cd5.q = 0.
cd5.s = 0.
cd5.s2 = 100
cd5.q2 = 0.01
cd5.I = cd5.m * cd5.L**2 / 12

#usv inialisation

#usv parameters fitting.
 
#bd.Fm = 1000.

bd.setl = 1000
bd.setw = 0
bd.ls = 5.
bd.mb = 350.
bd.pmax = 5000.
bd.mul2 = 100
bd.mut2 = 1000
bd.M = 0.
bd.theta = np.pi/2
#The usv is locate at the end of the boom

bd.pb = cd5.cms[:,-1] + bd.ls/2.*np.array([np.cos(bd.theta),np.sin(bd.theta)]) \
- cd5.L*cd5.para[:,-1] + [0,0.5 * cd5.cbl[1]]
#controller parameters fitting

#bd.krali = [10.,0.,0.]
#bd.krdis = [1,5.,0.1]
bd.krpid = [5.,10.,0.1]
bd.kvpid = [1000,450,250]
bd.cfr = 500.
#
#bd.kvali = [15., 0., 0.]
#bd.kvpid


#==============================================================================
# Dubins-like plan to enclose the anchored ship

#anchored ship position 
buque = np.array([100.,20.,-55.,10.,-np.pi]) 

#some contours points of the ship to draw it
vertices = np.array([[-buque[0]/2.,-buque[1]/2.],\
[-buque[0]/2.,buque[1]/2],\
[buque[0]/4.,buque[1]/2],\
[buque[0]/2,0],\
[buque[0]/4.,-buque[1]/2],\
[-buque[0]/2.,-buque[1]/2]])        
rot = np.array([[np.cos(buque[4]),- np.sin(buque[4])],[np.sin(buque[4]),\
np.cos(buque[4])]])

vertrot = np.array([np.dot(rot,j) for j in vertices]) \
+ [buque[2],buque[3]]

#Ship drawing                      
codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
Path.CURVE3]
   
pathi = Path(vertrot,codes)
patchi = patches.PathPatch(pathi,facecolor = 'black') #'red'
pl.gca().add_patch(patchi)

#primary way point for the usv that should enclose the ship

#we define an open dubins trajectory an select as the first point the current
#(starting point) of the usv

wp = [np.array([bd.pb[0],bd.pb[1],0,1]),\
np.array([6 + bd.pb[0],5 + bd.pb[1],6,1]),\
np.array([150,30,6,-1]),\
np.array([-100,30,6,-1]),\
np.array([-130,6,6,1]),\
np.array([-140,0,0,1])]

#secondary waypoints,

wp2 = dubing.secnwayp(wp)

#pintamos una linea de costa:
pl.plot([10,wp[-1][0]-10],[0,0],'k')

#pintamos la trayectoria a seguir

dubing.pintadubing(wp2)
pl.pause(1)

#==============================================================================

#inicializamos limites de representación de corrientes.#######################
# lim = np.array([]) #lim hay que definirlo como un array vacio si no hay
# corrientes
lim = [min([i[0] for i in wp]),max([i[0] for i in wp]) + 10,\
0,max(i[1] for i in wp) + 10]
#y definimos el tipo de corriente
c = 0 
##############################################################################

#tam = 1 #50000 #aqui definimos de momento el tamaño del experimento...
 #1000
# y aqui cada cuantas iteraciones guardamos datos y pintamos

tiempo = 1500. #1500 #total simulation time
#size of the arrays to save experimente results....
trecord = delta * 200000. #it should be a multiple of the integration step...
steps = int(tiempo//trecord) #total number of data to be saved
if steps <1:
    print "warning no results will be saved"
step = 0  #current step
bdrec = np.zeros((steps+1,16))
cadrec = np.zeros((steps+1,14,cd5.esl))

maxvueltas = 1 #numero de veces que se repite el circuito marcado por los wp.
#si se trata de un circuito con principio y fin lo lógico es que el valor de 
#maxvueltas sea 1
nvueltas = 0
#####################Main simulation loop#####################################
vset = 0.5 #usv speed setpoint
ant = np.array([0,1]) #i take the initial heading reference ad hoc... because
#i know were the two first waypoint are located...
while (t <= tiempo) & (nvueltas < maxvueltas):
 for i in wp2:
    print 'hola\n', i
    #always try to arrive to the first waypoint...
    vp = i[1] - bd.pb #vector from the ship actual position to the waypoint
    di = np.linalg.norm(vp) #intitial distance between the ship position 
                            #and the waypoint
    d = di #actual distance to be updated during the simulation.
    k = i[1] - ant
    k = k/np.linalg.norm(k) * vset  #speed vector (reference)
    #first step: reaching the first secondary waypoint before starting
    #the curve
    
    while d >= 3.0: #le damos 3 metros como margen de convergencia pero esto
                    #es una guarrada lamentable        
        bd.vr = - corrientes.gencor(bd.pb,fun = corrientes.campos[c])\
        + bd.vb
        bd.propulsion(delta,extf = - cd5.T[-2:],dfcm = [-bd.ls/2.,0])   
        cd5.mtr_s(bd)
                    
        #calculamos la velocidad de arrastre para la cadena
        cd5.vr = - corrientes.gencor(cd5.cms,fun = corrientes.campos[c])\
        + cd5.v
        #movemos la cadena...
        cd5.movimiento(delta = delta)    
        #y a continuación movemos los barcos          
        bd.movimiento(delta = delta)            
        bd.controlador(delta,\
        np.arctan2(vp[1],vp[0])*d/di + np.arctan2(k[1],k[0])* (1.-d/di),\
        norm(k),500)
       
        vp = i[1] - bd.pb
        d = norm(vp)
        di = d
        t += delta
        #print t, bd.pb, vp
        #print t
        if t >= step * trecord:
            #filling the matrix for saving results
            bdrec[step] = bd.pb[0],bd.pb[1],bd.vb[0],bd.vb[1],\
            bd.ab[0],bd.ab[1],\
            bd.theta,bd.wb,bd.alfab,\
            bd.Fm,bd.M,bd.setl,bd.setw,bd.thewj,cd5.T[-2],cd5.T[-1]
                
            cadrec[step] = cd5.cms[0],cd5.cms[1],cd5.v[0],cd5.v[1],cd5.a[0],\
            cd5.a[1],cd5.alfa,cd5.w,cd5.T[2:-1:2],cd5.T[3::2],cd5.normal[0],\
            cd5.normal[1],cd5.para[0],cd5.para[1]
#          
            step += 1            
        if t >= tiempo:
            break
        
    if t >= tiempo:
        break
    
    if i[0][2] == 0: #start and end way points has not to be turning around.
        ant = i[1]
        continue
    d = np.linalg.norm(i[2] - bd.pb)    
    vr = bd.pb-i[0][0:2] #vectorfrom the centre of dubin circle to the ship
    
    vt = np.array([-vr[1],vr[0]]) *(-i[0][3]) #normal vector in the moving sense
    ctr = i[0][2]-np.linalg.norm(vr)# error in the radious
    vrd = vr*ctr #adding size and sign
    
    vadd = vt + vrd
    acon = np.arctan2(vadd[1],vadd[0])
    
    print 'hola 2\n'
    while d >= 1.:
        
        bd.vr = - corrientes.gencor(bd.pb,fun = corrientes.campos[c])\
        + bd.vb
        bd.propulsion(delta,extf = - cd5.T[-2:],dfcm = [-bd.ls/2.,0])   
        cd5.mtr_s(bd)
                    
        #calculamos la velocidad de arrastre para la cadena
        cd5.vr = - corrientes.gencor(cd5.cms,fun = corrientes.campos[c])\
        + cd5.v
        #movemos la cadena...
        cd5.movimiento(delta = delta)    
        #y a continuación movemos los barcos          
        bd.movimiento(delta = delta)            
        bd.controlador(delta,\
        acon,\
        norm(k),500)
        
        
        
        d = norm(i[2] - bd.pb)
        
        vr = bd.pb-i[0][0:2] #vectorfrom the centre of dubin circle to the ship
    
        vt = np.array([-vr[1],vr[0]]) *(-i[0][3]) #normal vector in the moving sense
        ctr = i[0][2]-norm(vr)# error in the radious
        vrd = vr*ctr #adding size and sign
    
        vadd = vt + vrd
        acon = np.arctan2(vadd[1],vadd[0])
        t += delta
        if t >= step * trecord:
            #filling the matrix for saving results
            bdrec[step] = bd.pb[0],bd.pb[1],bd.vb[0],bd.vb[1],\
            bd.ab[0],bd.ab[1],\
            bd.theta,bd.wb,bd.alfab,\
            bd.Fm,bd.M,bd.setl,bd.setw,bd.thewj,cd5.T[-2],cd5.T[-1]
                
            cadrec[step] = cd5.cms[0],cd5.cms[1],cd5.v[0],cd5.v[1],cd5.a[0],\
            cd5.a[1],cd5.alfa,cd5.w,cd5.T[2:-1:2],cd5.T[3::2],cd5.normal[0],\
            cd5.normal[1],cd5.para[0],cd5.para[1]
#          
            step += 1
        if t>= tiempo:
            break
    ant = i[1]
    if t >= tiempo:
        break
 nvueltas += 1
#en la ultima fila guardamos propidades de los barcos y las cadenas y detalles
#del experimento       

bdrec = bdrec[0:step]
cadrec = cadrec[0:step]
        
bdrec[-1,0:13] = bd.thewjmax,bd.Ac,bd.mb,bd.mul,bd.Ib,bd.mua,bd.ls,bd.mut,\
bd.pmax,bd.pmax,bd.Ab,bd.Aw,c
        
cadrec[-1,0,0:3] = cd5.s,cd5.q,cd5.A
cadrec[-1,1,0:3] = cd5.L,cd5.m,cd5.I
cadrec[-1,2,0:3] = delta,step,tiempo
cadrec[-1,3,0:3] = cd5.cbl
      
np.savez(fichero, bd =bdrec, cadena =cadrec, bq = buque,\
pln = wp2,limites = lim)


def dibujar(bdcha,cadena,buque = [10.,2.,0.,0.,0.],pasos = 3):
    '''dibuja a partir de datos recogidos en arrays de datos tipo barco y cadena
    ver el sistema de preparar matrices para guardar datos dibuja también
    un buque fondeado buque = [eslora, manga, posicion_x,
    posicion_y, orientación,de     sup'''
    pl.figure(1)
    pl.hold(True)
    #       pintamos el barco que supuestamente queremos envolver
    vertices = np.array([[-buque[0]/2.,-buque[1]/2.],\
    [-buque[0]/2.,buque[1]/2],\
    [buque[0]/4.,buque[1]/2],\
    [buque[0]/2,0],\
    [buque[0]/4.,-buque[1]/2],\
    [-buque[0]/2.,-buque[1]/2]])        
    rot = np.array([[np.cos(buque[4]),- np.sin(buque[4])],[np.sin(buque[4]),\
    np.cos(buque[4])]])

    vertrot = np.array([np.dot(rot,j) for j in vertices]) \
    + [buque[2],buque[3]]
                       
    codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
    Path.CURVE3]
   
    pathi = Path(vertrot,codes)
    patchi = patches.PathPatch(pathi,facecolor = 'black') #'red'
    pl.gca().add_patch(patchi)

    for i in range(0,bdcha.shape[0]-1,pasos):    
        
        cms = np.array([cadena[i,0,:],cadena[i,1,:]])
        para = np.array([cadena[i,-2,:],cadena[i,-1,:]])
        pl.plot(cms[0,:],cms[1,:],'o')
        pl.hold(True)
        barrasi = cms + cadena[-1,1,0] * para
        barrasd = cms - cadena[-1,1,0] * para
        pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]])

                        
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
        
 #######################dibujar cable de arrastre derecha#######################
        rot = np.array([[np.cos(bdcha[i,6]),- np.sin(bdcha[i,6])],[np.sin(bdcha[i,6]),\
        np.cos(bdcha[i,6])]])  
        popad =  np.dot(rot, np.array([-bdcha[-1,6]/2.,0])) + [bdcha[i,0],bdcha[i,1]]  
        tipd = - para[:,-1] * cadena[-1,1,0] + cms[:,-1]
        
        
        distd = norm(popad - tipd)
        dd = distd/cadena[-1,3,1]
        print dd
        if dd > 1: dd = 1
        r = bezier_cvr.bezier4p([[tipd[0]],[tipd[1]]],[[popad[0]],[popad[1]]],1,1,1.5,\
        (1-dd) * bdcha[i,6]\
         +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
        (1-dd) * np.arctan2(-para[0,0],-para[0,1])\
         +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
        100)
        bezier_cvr.pintar_bezier(r[0],color = 'b')
###############################################################################       

        pl.plot(bdcha[i,0],bdcha[i,1],'+g')
        
#        pl.pause(0.01)        
        #hold(False) 
        
#        

dibujar(bdrec,cadrec,buque,pasos=1)

######para pintar  corrientes (incluido el 13.04.2016)#########################
corrientes.vercor(lim,[20,20],fun = corrientes.campos[c])
###############################################################################
#pl.figure(2)
#a =np.array([norm(i) for i in zip(bdrec[:-1,2],bdrec[:-1,3])])
##for i in zip(bdrec[:-1,2],bdrec[:-1,3]):
##    a.append(norm(i))    
##offset = np.concatenate(([0],np.cumsum([i[2][-1] for i in plan[:-1]])))
#
#tmp = np.array(np.arange(0,tiempo,2))
#pl.plot(tmp,a,'g')
##for i in plan:
##    pl.plot(i[2] + t0,norm(i[1],axis = 0))
##    t0 = i[2][-1]    
#
#pl.title('modulo de la velocidad')
#
#a = []
#
#for i in zip(cadrec[:,8,0],cadrec[:,9,0]):
#    a.append(norm(i)) 
#
#pl.figure(3)
#pl.plot(a,'r')
#
#a = []    
#for i in zip(cadrec[:,8,-1],cadrec[:,9,-1]):
#    a.append(norm(i)) 
#
#pl.plot(a,'g')
#
#pl.title('Tension en los extremos')
#    
#pl.figure(4)
#pl.plot(np.ones(bdrec.shape[0])*0.0)
#pl.plot(bdrec[:-1,6],'g')
#pl.title('rumbo')   
#
##situacion final
#pl.figure(5)
#pl.axis('equal')
#vertices = np.array([[-buque[0]/2.,-buque[1]/2.],\
#[-buque[0]/2.,buque[1]/2],\
#[buque[0]/4.,buque[1]/2],\
#[buque[0]/2,0],\
#[buque[0]/4.,-buque[1]/2],\
#[-buque[0]/2.,-buque[1]/2]])        
#rot = np.array([[np.cos(buque[4]),- np.sin(buque[4])],[np.sin(buque[4]),\
#np.cos(buque[4])]])
#
#vertrot = np.array([np.dot(rot,s) for s in vertices]) \
#+ [buque[2],buque[3]]
#codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
#Path.CURVE3]
#   
#pathb = Path(vertrot,codes)
#patchb = patches.PathPatch(pathb,facecolor = 'black') #'red'
#pl.gca().add_patch(patchb)
#pl.gca().add_patch(patchb)
#pl.plot([10,wp[-1][0]-10],[0,0],'k')
#cd5.dibujar('k')
#a = []
#for i in zip(birec[:-1,0],birec[:-1,1],bdrec[:-1,0],bdrec[:-1,1]):
#    a.append(norm([i[2]-i[0],i[3]-i[1]]))
#    
#pl.plot(a,'k')
#pl.title('distancia entre b')   