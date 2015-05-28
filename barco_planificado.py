# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 12:12:40 2015

@author: juan
Ejemplillo planificador bezier
"""
import barcosolo
import matplotlib.pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
from bezier_cvr import pintar_bezier
b = barcosolo.barco() #como no especificamos el gobierno es con waterjet

puntos = array([[[0.],[0.]],[[20.],[20.]],[[10.],[-20.]],[[0.],[0.]]]) #und so weiter...
rumbos = array([pi/2, 0., -pi, pi/2])


plan = b.planificador(puntos,rumbos,ndat = 50)


delta = 0.01 #paso de integracion
paso = 5

for i in plan:
    r = i[0]
    v = i[1]
    tp = i[2]
    t = 0.0
    pintar_bezier(r)
    for s,k,l in zip(r.transpose(),v.transpose(),tp):
#        print s, '\n'
#        print k, '\n'
        print l, '\n'
        pl. plot(s[0],s[1],'o')        
        vp = array([s[0] - b.pb[0],s[1] - b.pb[1]])
        di = norm(vp)
        d = di
        count = 0
        while t < l and dot(-vp,k)<= -0.0 and norm(vp) > 0.3:
            
            b.controlador(arctan2(vp[1],vp[0])*d/di\
            + arctan2(k[1],k[0])* (1.-d/di),\
            norm(k),delta)
            b.movimiento(delta)            
            if count%paso == 0:            
             pl.plot(b.pb[0],b.pb[1],'.')
             vp = array([s[0] - b.pb[0],s[1] - b.pb[1]])
             if count%(10*paso) == 0:            
               vertices = np.array([[-1.,-0.25],[-1.,0.25],[-0.25,0.35],[1,0],\
               [-0.25,-0.35],[-1.,-0.25]])
               rot = np.array([[np.cos(b.theta),- np.sin(b.theta)],[np.sin(b.theta),\
               np.cos(b.theta)]])  
               vertrot = np.array([np.dot(rot,j) for j in vertices]) +\
               [b.pb[0],b.pb[1]]       
               codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
               Path.CURVE3]
               pathi = Path(vertrot,codes)
               patchi = patches.PathPatch(pathi,facecolor = 'blue')
               pl.gca().add_patch(patchi)
             pause(0.1)
            t += delta
            count += 1
        
    
        
 




   
   

    
     