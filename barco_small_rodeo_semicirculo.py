# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:08:57 2014

Review: 20.04.2017. 
smal ship and small boom 
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
#wp = [np.array([0.,0.,0,-1]),\
#np.array([-21.177710264244599,3.8918224322146973,14.,-1]),\
#np.array([-21.177710264244599-14,0,0,-1])]
import read_log
wp = [np.array([0.,0.,0,-1]),\
np.array([-18.654246346914539,1.5,14.,-1]),\
np.array([-18.654246346914539-14,0,0,-1])]
#dubin secondary waypoints:
wp2 = dubing.secnwayp(wp)
dubing.pintadubing(wp2)
pl.axis('equal')
pl.pause(1)

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
cd5 = boom.cadena(50,np.array([0,4,1000]),np.array([1,0])) #to moore the tip,np.array([1,0]))
cd5.L = 0.415
#the boom is collocated folded 
# the links are collocated parallel to the dock
up = 1
for i in range(cd5.esl):
    cd5.normal[:,i] = [0,up]
    cd5.para[:,i] = [up,0]
    up = - up
    
#except the last one which is located tilted  45º with the dock to attach it to
#the usv    
#cd5.normal[:,0] = [1,0]
#cd5.normal[:,-1] = [-1,0]
#cd5.para[:,0] = [0,1]
#cd5.para[:,-1] = [0,-1]

cd5.calcms(np.array([0,-1.2]),np.array([0.-1.2]))

#left tip is free Fi[0] = 0 #the left boom  tip is attached to the dock
# Fi = -np.infty
cd5.m = 0.05
cd5.A = 0.01
cd5.q = 0.
cd5.s = 0.
cd5.s2 = 0.05
cd5.q2 = 0.001
cd5.I = cd5.m * cd5.L**2 / 12

#usv inialisation

#usv parameters fitting.
 
#bd.Fm = 1000.

bd.ls = 1.
bd.Ac = 0.4 
bd.mb = 3.4
bd.Ib = 0.85
bd.mua = 0. #caveat in the equations damping is bd.mua*bd.ls*bd.wb  
bd.mua2 = 1.#9.7577 #caveat in the equations damping is bd.mua*bd.ls*bd.wb**2 
bd.mut2 = 4. #2.226
bd.mut = 0.1#0.1
bd.mul2= 1.219
bd.mul = 0.01
bd.krpid = [2,0.5,0]
bd.kvpid = [0.,0.,0.]#np.array([5,0.,5])
bd.krd = [0.2,0,0]
bd.krl = [0.,0.,0.]
#########################valores cte pid para barco libre#####################
#bd.krpid = [1.5,0.5,0]
#bd.kvpid = [0.,0.,0.]#np.array([5,0.,5])
#bd.krd = [0.35,0,0]
#bd.krl = [0.,0.,0.]
##############################################################################
bd.thewjmax = 20. * np.pi / 180.
bd.pmax = 10.
#bd.theta = np.pi*5/6 #con este funcionaba ok

bd.theta = read_log.heading[0]
#bd.thewj = 20 * pi /180


#The usv is locate at the end of the boom

#bd.pb = cd5.cms[:,-1] + bd.ls/2.*np.array([np.cos(bd.theta),np.sin(bd.theta)]) \
#- cd5.L*cd5.para[:,-1] + [0,0.5 * cd5.cbl[1]]
#bd.pb = np.array([50.,0.])
bd.pb = np.array(read_log.datos_m).T[0]#np.array([0.,0])
#vset = 1 #0.5 #usv speed setpoint
ffv = 0.5 #fijamos la potencia transmitida al motor en lugar de controlar la velocidad
#===============Con control de distancia como si tuviera dos barcos=========#
#Es una trampa, el segundo barco no existe pero usamos el control 
#de distancia entre dos barcos que arrastran juntos una red. aquí la distacia
#que se controla es la del barco a la curva de Dubins. 18.05.2017

bd.rlink = [[bd.ide,bd.pb,bd.theta,bd.vb],[0,np.zeros(2),0.,np.zeros(2)]]

#controller parameters fitting

#bd.krali = [10.,0.,0.]
#bd.krdis = [1,5.,0.1]
#bd.krpid = [5.,10.,0.1]
#bd.kvpid = [1000,450,250]
#bd.cfr = 500.
#
#bd.kvali = [15., 0., 0.]
#bd.kvpid


#==============================================================================
# Simple Dubins-like plan 

#anchored ship position 
buque = np.array([])
#buque = np.array([8.,2.,-5.,1.,-np.pi]) 
#
##some contours points of the ship to draw it
#vertices = np.array([[-buque[0]/2.,-buque[1]/2.],\
#[-buque[0]/2.,buque[1]/2],\
#[buque[0]/4.,buque[1]/2],\
#[buque[0]/2,0],\
#[buque[0]/4.,-buque[1]/2],\
#[-buque[0]/2.,-buque[1]/2]])        
#rot = np.array([[np.cos(buque[4]),- np.sin(buque[4])],[np.sin(buque[4]),\
#np.cos(buque[4])]])
#
#vertrot = np.array([np.dot(rot,j) for j in vertices]) \
#+ [buque[2],buque[3]]
#
##Ship drawing                      
#codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
#Path.CURVE3]
#   
#pathi = Path(vertrot,codes)
#patchi = patches.PathPatch(pathi,facecolor = 'black') #'red'
#pl.gca().add_patch(patchi)




#==============================================================================

#inicializamos limites de representación de corrientes.#######################
# lim = np.array([]) #lim hay que definirlo como un array vacio si no hay
# corrientes
lim = [min([i[0] for i in wp]),max([i[0] for i in wp]),\
0,max(i[1] for i in wp)+15]
#y definimos el tipo de orriente
c = 0#2#1
#y las pintamos
corrientes.vercor(lim,[20,20],fun = corrientes.campos[c])
##############################################################################

delta = 0.000001 #0.0001
t= 0.
tiempo =120 #ojo 80 solo para corriente a favor
#tiempo= 130. #600. #total simulation time-> it cannot be an INTEGERRRRR
#size of the arrays to save experimente results....
trecord = delta * 50000. #it should be a multiple of the integration step...or not
tpintar = trecord * 20.
steps = int(tiempo//trecord) #total number of data to be saved
if steps <1:
    print "warning no results will be saved"
step = 0  #current step
step_pintar = 0
bdrec = np.zeros((steps+1,16))
cadrec = np.zeros((steps+1,14,cd5.esl))

maxvueltas = 1 #numero de veces que se repite el circuito marcado por los wp.
#si se trata de un circuito con principio y fin lo lógico es que el valor de 
#maxvueltas sea 1
nvueltas = 0
#####################Main simulation loop#####################################

ant = np.array([0.,0.])#np.array([0,1]) #i take the initial heading reference ad hoc... because
#i know were the two first waypoint are located...

while (t <= tiempo) & (nvueltas < maxvueltas):
 
 for i in wp2: 
    print 'hola\n', i
    
    p = bd.pb - i[1] #vector from the ship actual position to the waypoint
    
    d = np.linalg.norm(p) #distance between the ship position 
                            #and the waypoint
    dmin = d
    #first step: reaching the first secondary point before starting the curve
    #actual distance to be updated during the simulation.
    k = i[1] - ant #vector from the previous waypoint to the actual desired wp
    k = k/(np.linalg.norm(k) + (np.linalg.norm(k) == 0)) # deseared heading direction (unit vector)
    
    dop = np.dot(p,k) #the rationale behind this condition is: each waypoint
    #has to be overtaken in the heading direction, 
    #so once a waypoint  is outran there is not point on coming back to it
    dp = dop * k
    dl = p + dp  #this should be a vector  normal to the path...
    
    s = p + dl #just try to approach the path, keeping a lookout (to the next wp)    
    ###############################distance to path control terms ###########
    bd.rlink[1][1:3] = [i[1] + dp,np.arctan2(dp[1],dp[0])]
    #########################################################################
    while d > 3.:        
        bd.vr = - corrientes.gencor(bd.pb,fun = corrientes.campos[c])\
        + bd.vb
        bd.propulsion(delta,extf = - cd5.T[-2:],dfcm = [-bd.ls/2.,0])   
#        bd.propulsion(delta)
        cd5.mtr_s(bd)
                    
        #calculamos la velocidad de arrastre para la cadena
        cd5.vr = - corrientes.gencor(cd5.cms,fun = corrientes.campos[c])\
        + cd5.v
        #movemos la cadena...
        cd5.movimiento(delta = delta)    
        #y a continuación movemos los barcos          
        bd.movimiento(delta = delta)            
        
        p = bd.pb - i[1]
        d = np.linalg.norm(p)
        if d < dmin:
         dmin = d
        dop = np.dot(p,k)
        dp =dop * k
        dl = p + dp
        s = p + dl        
        bd.rlink[1][1:3] = [i[1]+ dp,np.arctan2(dp[1],dp[0])]
        bd.controlador(delta, np.arctan2(k[1],k[0]),0,ffv,np.zeros(2))
       
        #print bd.thewj * 180/np.pi
        
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
            
            step += 1 
####para pintar en continuo no tiene sentido en simulaciones grandes #########
        if t >= step_pintar * tpintar:    
            #print bd.thewj * 180 /np.pi, bd.cstb(bd.vb)
            
            pl.plot(bd.rlink[1][1][0],bd.rlink[1][1][1],'ob')
            bd.dibujar()
            cd5.dibujar()
            
            popad =  bd.csbt(np.array([-bd.ls/2.,0])) + bd.pb  
            tipd = - cd5.para[:,-1] * cd5.L + cd5.cms[:,-1]
            
            
            distd = norm(popad - tipd)
            dd = distd/cd5.cbl[1]
            
            if dd > 1: dd = 1
            r = bezier_cvr.bezier4p([[tipd[0]],[tipd[1]]],[[popad[0]],\
            [popad[1]]],1,1,1.5,\
            (1-dd) * bd.theta\
             +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
            (1-dd) * np.arctan2(-cd5.para[0,0],-cd5.para[0,1])\
             +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
            100)
            bezier_cvr.pintar_bezier(r[0],color = 'b')
            pl.pause(1)
            print bd.thewj * 180/np.pi
            
            step_pintar += 1 
###############################################################################
                        
           
        if t >= tiempo:
            break
#Nuevo invento peligroso: si me alejo del waypoint, me voy al siguiente..######        
        if d > dmin:
            break
###############################################################################        
    if t >= tiempo:
        break
    
    if i[0][2] == 0: #start and end waypoints have not to be turned around.
        ant = i[1]
        continue
    
 ####################We are here fiddling around #########################   
        
    vr = bd.pb-i[0][0:2] #vectorfrom the centre of dubin circle to the ship
    ctr = i[0][2]-np.linalg.norm(vr)# error in the radious    
    vr = vr/(np.linalg.norm(vr)  + (np.linalg.norm(vr) == 0))  #normal vector from the center of circle to ship
    vt = np.array([-vr[1],vr[0]]) *(-i[0][3]) #normal vector in the moving sense
    
    
    vrd = vr*ctr 
###########################distance to path control terms###################
    bd.rlink[1][1:3] = [bd.pb + vrd,np.arctan2(vt[1],vt[0])]

############################################################################    
    vadd = vt + vrd
#    acon = np.arctan2(vadd[1],vadd[0])
    acon = np.arctan2(vt[1],vt[0])
    d = norm(i[2] - bd.pb)
    dmin = d 
    print 'hola 2\n'
    while d >= 3.:
        
        bd.vr = - corrientes.gencor(bd.pb,fun = corrientes.campos[c])\
        + bd.vb
        bd.propulsion(delta,extf = - cd5.T[-2:],dfcm = [-bd.ls/2.,0]) 
#        bd.propulsion(delta)   
        cd5.mtr_s(bd)
                    
        #calculamos la velocidad de arrastre para la cadena
        cd5.vr = - corrientes.gencor(cd5.cms,fun = corrientes.campos[c])\
        + cd5.v
        #movemos la cadena...
        cd5.movimiento(delta = delta)    
        #y a continuación movemos los barcos          
        bd.movimiento(delta = delta)
            
        bd.controlador(delta,acon,0.,ffv,np.zeros(2))
        
        
        
        d = norm(i[2] - bd.pb)
        if d < dmin:
            dmin = d
        vr = bd.pb-i[0][0:2] #vectorfrom the centre of dubin circle to the ship
        ctr = i[0][2]-norm(vr)# error in the radious
        vr = vr/(np.linalg.norm(vr) + (np.linalg.norm(vr) == 0))
        vt = np.array([-vr[1],vr[0]]) *(-i[0][3]) #normal vector in the moving sense
        
        vrd = vr*ctr #adding size and sign
###########################distance to path control terms#####################        
        bd.rlink[1][1:3] = [bd.pb + vrd,np.arctan2(vt[1],vt[0])]

##############################################################################    
#        vadd = vt + vrd
#        acon = np.arctan2(vadd[1],vadd[0])
        acon = np.arctan2(vt[1],vt[0])
###############################################################################        
#adding the heading correction to keep the bow pointing inside the circle        
            #real heading calculation 
        den = i[0][2]**2*bd.Ac*bd.mut2 - bd.mua2
        if den <0 :
            print 'infeasible trajectory taking a wider radious.'
            den = 2
        the = np.arctan(np.sqrt(bd.mua2/den)) * -i[0][3]
#        rot = np.array([[np.cos(the),-np.sin(the)],[np.sin(the),np.cos(the)]])
        acon = acon + the
###############################################################################        
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
            
            step += 1 
####para pintar en continuo no tiene sentido en simulaciones grandes #########
        if t >= step_pintar * tpintar:    
            #print bd.thewj * 180 /np.pi, bd.cstb(bd.vb)
            
            pl.plot(bd.rlink[1][1][0],bd.rlink[1][1][1],'ob')
            bd.dibujar()
            cd5.dibujar()
            
            popad =  bd.csbt(np.array([-bd.ls/2.,0])) + bd.pb  
            tipd = - cd5.para[:,-1] * cd5.L + cd5.cms[:,-1]
            
            
            distd = norm(popad - tipd)
            dd = distd/cd5.cbl[1]
            
            if dd > 1: dd = 1
            r = bezier_cvr.bezier4p([[tipd[0]],[tipd[1]]],[[popad[0]],\
            [popad[1]]],1,1,1.5,\
            (1-dd) * bd.theta\
             +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
            (1-dd) * np.arctan2(-cd5.para[0,0],-cd5.para[0,1])\
             +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
            100)
            bezier_cvr.pintar_bezier(r[0],color = 'b')
            pl.pause(1)
            print bd.thewj * 180/np.pi
            
            step_pintar += 1 
###############################################################################          
           
        if t>= tiempo:
            break
###############nuevo invento peligro si me alejo del waypoin corto############        
        if d > dmin:
            break
###############################################################################        
    ant = i[2]#i[1]
    if t >= tiempo:
        break
 nvueltas += 1   
#en la ultima fila guardamos propidades de los barcos y las cadenas y detalles
#del experimento       
        
bdrec[-1,0:13] = bd.thewjmax,bd.Ac,bd.mb,bd.mul,bd.Ib,bd.mua,bd.ls,bd.mut,\
bd.pmax,bd.pmax,bd.Ab,bd.Aw,c
        
cadrec[-1,0,0:3] = cd5.s,cd5.q,cd5.A
cadrec[-1,1,0:3] = cd5.L,cd5.m,cd5.I
cadrec[-1,2,0:3] = delta,step,tiempo
cadrec[-1,3,0:3] = cd5.cbl
      
np.savez(fichero, bd =bdrec, cadena =cadrec, bq = buque,\
pln = wp2,limites = lim)


def dibujar(bdcha,cadena,buque = [10.,2.,0.,0.,0.]):
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

    #forma del barco remolcador
    vertices = np.array([[-bdcha[-1,6]/2.,-0.2*bdcha[-1,6]],\
    [-bdcha[-1,6]/2.,0.2*bdcha[-1,6]],\
    [-0.2*bdcha[-1,6],0.3*bdcha[-1,6]],\
    [bdcha[-1,6]/2.,0],[-0.2*bdcha[-1,6],-0.3*bdcha[-1,6]],\
    [-bdcha[-1,6]/2.,-0.2*bdcha[-1,6]]]) 

    for i in range(bdcha.shape[0]-1):    
        
        cms = np.array([cadena[i,0,:],cadena[i,1,:]])
        para = np.array([cadena[i,-2,:],cadena[i,-1,:]])
        pl.plot(cms[0,:],cms[1,:],'o')
        pl.hold(True)
        barrasi = cms + cadena[-1,1,0] * para
        barrasd = cms - cadena[-1,1,0] * para
        pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]])

                        
#        vertices = np.array([[-bdcha[-1,6]/2.,-0.25],[-bdcha[-1,6]/2.,0.25],\
#        [-0.25,0.35],[bdcha[-1,6]/2.,0],[-0.25,-0.35],[-bdcha[-1,6]/2.,-0.25]])        
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
        #print dd
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

##dibujar(bdrec,cadrec,buque)
#dibujar(bdrec,cadrec)
#######para pintar  corrientes (incluido el 13.04.2016)#########################
#corrientes.vercor(lim,[20,20],fun = corrientes.campos[c])
#################################################################################
#pl.figure(2)
#a =np.array([norm(i) for i in zip(bdrec[:-1,2],bdrec[:-1,3])])
##for i in zip(bdrec[:-1,2],bdrec[:-1,3]):
##    a.append(norm(i))    
##offset = np.concatenate(([0],np.cumsum([i[2][-1] for i in plan[:-1]])))
#
#tmp = np.array(np.arange(0,tiempo,trecord))
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

#situacion final
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