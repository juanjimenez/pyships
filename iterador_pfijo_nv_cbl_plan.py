# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:08:57 2014
Este programa esta pensado para iterar cadenas y barcos incluyendo
control de rumbo y velocidad y ver si peta o que
pasa

Review: 06.03.2017. (Simulation for a short boom to enclose a 9m lenght ship)

@author: juan
"""
import usv
import boom
import corrientes
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

delta = 0.00001
tiempo= 0
steps = 20 #20 #pasos en que se divide cada tramo planificado para calcular los
            #datos de la liebre

#fichero = strftime("%b%d_%Y_%H%M") para guardar resultados
fichero = "../"+strftime("%b%d_%Y").lower() + "/" +\
strftime("%b%d_%Y_%H%M").lower()
directorio = path.dirname(fichero)
try:
    stat(directorio)
except:
    mkdir(directorio)

#creamos los barcos
#bi = barcosolo_nv.barco(1) #,'timon')
bd = usv.barco(2) #,'timon')

#creamos la cadena de momento pequeñita porque es para probar
cd5 = boom.cadena(10,np.array([0,1,1000]))
cd5.L = 0.5
#situamos la cadena ... plegada en el punto cero,excepto el ulitmo eslabon

#apuntando para arriba.
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

#bd.pb = cd5.cms[:,-1] + bd.ls/2.*np.array([np.cos(bd.theta),np.sin(bd.theta)])\
#- cd5.L*cd5.para[:,-1]

bd.pb = cd5.cms[:,-1] + bd.ls/2.*np.array([np.cos(bd.theta),np.sin(bd.theta)]) \
- cd5.L*cd5.para[:,-1] + [0,0.5 * cd5.cbl[1]]
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
bd.kvpid = [1000,450,250]
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
buque = np.array([10.,2.,-7.,1.,-np.pi]) 

#puntos contorno del barco para pintarlo
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

#pintamos el buque                      
codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
Path.CURVE3]
   
pathi = Path(vertrot,codes)
patchi = patches.PathPatch(pathi,facecolor = 'black') #'red'
pl.gca().add_patch(patchi)

#puntos de paso del barco que arrastra la cadena
xmax = np.max(vertrot[:,0])
xmin = np.min(vertrot[:,0])
ymax = np.max(vertrot[:,1])
ymin = np.min(vertrot[:,1])
puntos = np.array([[[bd.pb[0]],[bd.pb[1]+0.5]],\
[[xmax+6],[ymax+6.]],\
[[xmin+2],[ymax+4.]],[[xmin-10.],[ymin+1]]])
rumbos = np.array([np.pi/2,np.pi/2,-np.pi*5./8.,-np.pi])
veloc = np.array([0.4,0.4,0.4,0.4])
plan = bd.planificador(puntos[0:3],rumbos,veloc = veloc[0:3],ndat = steps)
#pintamos una linea de costa:
pl.plot([10,xmin-10],[0,0],'k')

#pintamos la trayectoria a seguir

for i in plan:
    bezier_cvr.pintar_bezier(i[0])
pl.pause(1)

#==============================================================================

#inicializamos limites de representación de corrientes.#######################
# lim = np.array([]) #lim hay que definirlo como un array vacio si no hay corrientes
lim = [-25,10,0,15]
#y definimos el tipo de corriente
c = 2
##############################################################################

#tam = 1 #50000 #aqui definimos de momento el tamaño del experimento...
 #1000
# y aqui cada cuantas iteraciones guardamos datos y pintamos


bdrec = np.zeros((len(plan)*steps+1,16))
cadrec = np.zeros((len(plan)*steps+1,14,cd5.esl))


j = 0 #indice de las matrices que guardan los resultados
paso = 5
for i in plan:
    r = i[0]
    v = i[1]
    tp = i[2]
    t = 0.0    
    for s,k,l in zip(r.transpose(),v.transpose(),tp):                
        vp = np.array([s[0] - bd.pb[0],s[1] - bd.pb[1]])
        di = norm(vp)
        d = di
        count = 0 
        while t < l and np.dot(-vp,k)<= -0.0 and norm(vp) > 0.3:
            bd.vr = - corrientes.gencor(bd.pb,corrientes.campos[c]) + bd.vb
            bd.propulsion(delta,extf = - cd5.T[-2:],dfcm = [-bd.ls/2.,0])   

            cd5.mtr_s(bd)
                    
            #calculamos la velocidad de arrastre para la cadena
            cd5.vr = - corrientes.gencor(cd5.cms,corrientes.campos[c]) + cd5.v    
            cd5.movimiento(delta = delta)    
            #y a continuación movemos los barcos          
            bd.movimiento(delta = delta)            
            bd.controlador(delta,\
            np.arctan2(vp[1],vp[0])*d/di + np.arctan2(k[1],k[0])* (1.-d/di),\
            norm(k),500)
            
            if count%paso == 0:            
#             pl.plot(b.pb[0],b.pb[1],'.')
             vp = np.array([s[0] - bd.pb[0],s[1] - bd.pb[1]])
            t += delta
            count += 1                            
            
        bdrec[j] = bd.pb[0],bd.pb[1],bd.vb[0],bd.vb[1],\
        bd.ab[0],bd.ab[1],\
        bd.theta,bd.wb,bd.alfab,\
        bd.Fm,bd.M,bd.setl,bd.setw,bd.thewj,cd5.T[-2],cd5.T[-1]
                
        cadrec[j] = cd5.cms[0],cd5.cms[1],cd5.v[0],cd5.v[1],cd5.a[0],\
        cd5.a[1],cd5.alfa,cd5.w,cd5.T[2:-1:2],cd5.T[3::2],cd5.normal[0],\
        cd5.normal[1],cd5.para[0],cd5.para[1]
        j = j + 1
    tiempo += t
#en la ultima fila guardamos propidades de los barcos y las cadenas y detalles
#del experimento       
        
bdrec[-1,0:13] = bd.thewjmax,bd.Ac,bd.mb,bd.mul,bd.Ib,bd.mua,bd.ls,bd.mut,\
bd.pmax,bd.pmax,bd.Ab,bd.Aw,c
        
cadrec[-1,0,0:3] = cd5.s,cd5.q,cd5.A
cadrec[-1,1,0:3] = cd5.L,cd5.m,cd5.I
cadrec[-1,2,0:3] = delta,steps,tiempo
cadrec[-1,3,0:3] = cd5.cbl
      
np.savez(fichero, bd =bdrec, cadena =cadrec, bq = buque,\
pln = plan,limites = lim)


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

    for i in range(bdcha.shape[0]-1):    
        
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

dibujar(bdrec,cadrec,buque)

######para pintar  corrientes (incluido el 13.04.2016)#########################
corrientes.vercor(lim,[20,20],corrientes.campos[c])
###############################################################################
pl.figure(2)
a =np.array([norm(i) for i in zip(bdrec[:-1,2],bdrec[:-1,3])])
#for i in zip(bdrec[:-1,2],bdrec[:-1,3]):
#    a.append(norm(i))    
offset = np.concatenate(([0],np.cumsum([i[2][-1] for i in plan[:-1]])))
tmp = np.concatenate([i[2]+k for i,k in zip(plan,offset)])
vel = np. concatenate([norm(i[1],axis = 0) for i in plan])
pl.plot(tmp,vel)
pl.plot(tmp,a,'g')
#for i in plan:
#    pl.plot(i[2] + t0,norm(i[1],axis = 0))
#    t0 = i[2][-1]    

pl.title('modulo de la velocidad')

a = []

for i in zip(cadrec[:,8,0],cadrec[:,9,0]):
    a.append(norm(i)) 

pl.figure(3)
pl.plot(a,'r')

a = []    
for i in zip(cadrec[:,8,-1],cadrec[:,9,-1]):
    a.append(norm(i)) 

pl.plot(a,'g')

pl.title('Tension en los extremos')
    
pl.figure(4)
pl.plot(np.ones(bdrec.shape[0])*0.0)
pl.plot(bdrec[:-1,6],'g')
pl.title('rumbo')   

#situacion final
pl.figure(5)
pl.axis('equal')
vertices = np.array([[-buque[0]/2.,-buque[1]/2.],\
[-buque[0]/2.,buque[1]/2],\
[buque[0]/4.,buque[1]/2],\
[buque[0]/2,0],\
[buque[0]/4.,-buque[1]/2],\
[-buque[0]/2.,-buque[1]/2]])        
rot = np.array([[np.cos(buque[4]),- np.sin(buque[4])],[np.sin(buque[4]),\
np.cos(buque[4])]])

vertrot = np.array([np.dot(rot,s) for s in vertices]) \
+ [buque[2],buque[3]]
codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
Path.CURVE3]
   
pathb = Path(vertrot,codes)
patchb = patches.PathPatch(pathb,facecolor = 'black') #'red'
pl.gca().add_patch(patchb)
pl.gca().add_patch(patchb)
pl.plot([10,xmin-10],[0,0],'k')
cd5.dibujar('k')
#a = []
#for i in zip(birec[:-1,0],birec[:-1,1],bdrec[:-1,0],bdrec[:-1,1]):
#    a.append(norm([i[2]-i[0],i[3]-i[1]]))
#    
#pl.plot(a,'k')
#pl.title('distancia entre b')   