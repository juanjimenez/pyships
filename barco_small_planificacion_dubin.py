# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 13:50:24 2017
Example o control with Dubing curves
@author: juan
"""
import dubing
import usv
import numpy as np
import matplotlib.pyplot as pl

#close('all')
#dubin primary way points:
wp = [np.array([0,0,10,1]),np.array([80,0,10,-1])] #,np.array([0,40,10,1]),\
#np.array([40,40,10,-1])]

#dubin secondary waypoints:
wp2 = dubing.secnwayp(wp)
dubing.pintadubing(wp2)
#A usb to follow the waypoints:
b =usv.barco(1)
b.ls = 1
b.Ac = 0.4 
b.mb = 3.4
b.Ib = 0.85
b.mua = 0 #caveat in the equations damping is b.mua*b.ls*b.wb  
b.mua2 = 9.7577 #caveat in the equations damping is b.mua*b.ls*b.wb**2 
b.mut2 = 2.226
b.mut = 0.1
b.mul2= 1.219
b.mul = 0.01
b.krpid = [5.,0.5,0]
b.kvpid = np.array([5,0.,5])
b.pb = np.array([40.,-20.])
b.thewjmax = 20 * np.pi / 180
b.pmax = 10
b.theta = 0.
vset = 1.8 #we fixed a constant speed to be reached and hold during  the trip 
#I should send the ship to the first way point, but... I wonder if this is 
#so important if it has to follow a close path...
delta = 0.01 #integration time step
ant = wp2[-1][2] #for close trajectories the initial heading is from the last
#secondary waypoint to the first one...
step = 0
tiempo = 600
t = 0
while (t <= tiempo):
    
 for i in wp2:
    print 'hola\n'
    #always try to arrive to the first waypoint...
    vp = i[1] - b.pb #vector from the ship actual position to the waypoint
    di = np.linalg.norm(vp) #intitial distance between the ship position 
                            #and the waypoint
    d = di #actual distance to be updated during the simulation.
    k = i[1] - ant
    k = k/np.linalg.norm(k) * vset  #speed vector (reference)
    #first step: reaching the first secondary waypoint before starting
    #the curve
    n_u = np.array([-k[1],k[0]]) #normal vector pointing towards the line
    n_u = n_u/np.linalg.norm(n_u) #unit normal vector
    dl = abs(np.dot(vp,n_u))*n_u # straight distance from ship to trajectory
    s = vp #+ dl #the angle of thsi vector is uses a set point for straight paths
    pl.figure(1)
    pl.plot(i[1][0],i[1][1],'or')
    while d >= 1.:
        b.controlador(delta,\
        np.arctan2(s[1],s[0]),\
        np.linalg.norm(k))
        b.propulsion(delta)
        b.movimiento(delta)
        vp = i[1] - b.pb
        d = np.linalg.norm(vp)
        di = d
        dl = np.dot(vp,n_u)*n_u # straight distance from ship to trajectory
        s = vp + 5.* dl #the angle of thsi vector is the setpoint for straight paths
        t += delta
        if t >= step:
           pl.figure(1)
           pl.axis('equal')
           b.dibujar()
#           arrow(b.pb[0],b.pb[1],s[0],s[1])
#           arrow(b.pb[0],b.pb[1],vp[0],vp[1],color = 'r')
#           arrow(b.pb[0],b.pb[1],dl[0],dl[1],color = 'b')
           pl.figure(2)
           pl.plot(t,np.linalg.norm(b.vb),'o')
           step += 1
           pl.pause(0.01)

        if t >= tiempo:
             break
#   #################
#    t = tiempo + 10.
#   #################
#    break
    pl.figure(1)        
    pl.plot(i[1][0],i[1][1],'og')
    pl.plot(i[2][0],i[2][1],'or')
    
   

    #second step going around the primary waypoint#############################
    d = np.linalg.norm(i[2] - b.pb)    
    vr = b.pb-i[0][0:2] #vectorfrom the centre of dubin circle to the ship
    vt = np.array([-vr[1],vr[0]]) *(-i[0][3]) #normal vector in the moving sense
    
    
    #real heading calculation 
    den = i[0][2]**2*b.Ac*b.mut2 - b.mua2
    if den <0 :
        print 'infeasible trajectory taking a wider radious.'
        den = 2
    the = np.arctan(np.sqrt(b.mua2/den)) * -i[0][3]
    rot = np.array([[np.cos(the),-np.sin(the)],[np.sin(the),np.cos(the)]])
    vt = np.dot(rot,vt)
    
    ctr = i[0][2]-np.linalg.norm(vr)# error in the radious
    
    
    vrd = vr*ctr/i[0][2] #adding size and sign
    
    vadd = vt + vrd 
    acon = np.arctan2(vadd[1],vadd[0])
    #cotroluto a lo bruto
    
#    normal = np.arctan2(vr[1],vr[0]) # vr angle    
#    hdep = np.arctan2(vr[0],-vr[1])
#    hdep = hdep- np.pi * np.sign(hdep) * (i[0][3] > 0) #heading set point
#    ctr = i[0][2]-np.linalg.norm(vr)# error in the radious
#    actr = abs(ctr)
#    normal = normal - np.pi * np.sign(normal) * (ctr < 0)
#    
#    #para el promedio hay que tener todos los angulos en el primer circulo.
#    normal = normal - 2 * np.pi * (normal < 0)  
#    hdep = hdep + 2 * np.pi * (hdep < 0)    
    
    #acon = normal   * actr/(1 + actr) + hdep/(1 + actr)
    
    #acon = acon - 2 * np.pi * (abs(acon) > np.pi) * np.sign(acon)
    while d >= 1.:
        
        b.controlador(delta,\
        acon,\
        np.linalg.norm(k))
        b.propulsion(delta)
        b.movimiento(delta)
        d = np.linalg.norm(i[2] - b.pb)
#        vr = b.pb-i[0][0:2]
        
         
#        normal = np.arctan2(vr[1],vr[0]) # vr angle  
#        print 'normal inicial', normal
#        hdep = np.arctan2(vr[0],-vr[1])
#        print 'hdep inicial', hdep
#        hdep = hdep - np.pi * np.sign(hdep) * (i[0][3] > 0) #heading set point
#        ctr = i[0][2]-np.linalg.norm(vr)# error in the radious
#        actr = abs(ctr)
#        normal = normal - np.pi * np.sign(normal) * (ctr < 0)
#        print 'ctr ',ctr, ' normal ', normal    
        #para el promedio hay que tener todos los angulos en el primer circulo.
        #normal = normal + 2 * np.pi * (normal < 0)  
        
        #hdep = hdep + 2 * np.pi * (hdep < 0)    
        #print 'normal y hdep p circulo', normal, hdep   
        
        
#        acon = normal   * actr/(1 + actr) + hdep/(1 + actr)
#           
#        acon = acon - 2 * np.pi * (abs(acon) > np.pi) * np.sign(acon)
        
        vr = b.pb-i[0][0:2] #vectorfrom the centre of dubin circle to the ship
    
        vt = np.array([-vr[1],vr[0]]) *(-i[0][3]) #normal vector in the moving sense
        #real heading calculation        
        vt = np.dot(rot,vt)/i[0][2]
        
        ctr = i[0][2]-np.linalg.norm(vr)# error in the radious
        vrd = vr*ctr/i[0][2] #adding size and sign
         
        
       
        vadd = vt + vrd 
        
        acon = np.arctan2(vadd[1],vadd[0])
        #print acon, b.setw,  b.thewj, b.thewjmax      
        #cotroluto a lo bruto
        
        t +=delta
        if t >= step:
            pl.figure(1)
            b.dibujar()
            pl.arrow(b.pb[0],b.pb[1],vt[0],vt[1])
            pl.arrow(b.pb[0],b.pb[1],vadd[0],vadd[1],color='r')
            pl.arrow(b.pb[0],b.pb[1],vrd[0],vrd[1],color='g')
            #pl.arrow(i[0][0],i[0][1],vr[0],vr[1],color='b')
            pl.figure(2)
            pl.plot(t,np.linalg.norm(b.vb),'o')
            step += 1
            pl.pause(0.01)
        
        if t >= tiempo:
            break
    pl.figure(1)        
    pl.plot(i[2][0],i[2][1],'og')        
    ant = i[1]
    if t >= tiempo:
        break
    