# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 12:32:53 2017
This files contains all functions employed to describe a trajectory following
a set of waypoints. The trayectory is eventually saved as a list of sucessive
waypoints the ship should cover to fufill the trayectory. There is a little
bit of confussion ir our notation; there are two kinds of waypoints:

Primary waypoints:
Those defined to build the trayectory. It is supposed the ship are 
going to shirt the way point describing a circle of radius r around it.
The direction of rotation around the way point should be also
supplied. So, a primary wypoint is a tupla wp [x,y,r,s] x and y coordinates of
the waypoint, circle radius, s=+1 rotate clockwise s= -1
rotate counter-clockwise. 

Secodary waypoints:

They are calculated as the end points of the tangent segment between
the two circles defined by two consecutive primary waypoints.
See below for details...

@author: juan
"""

import numpy as np
import matplotlib.pyplot as pl

def tancircles(c1,c2):
    """ 
    c1 and c2  are numpy arrays, defining two circles, with elements [x,y,r]:
    x and y are the coordinates of the circle center and r is the radius
    The program return an array: with the coordenates of the intersection 
    points
    between each circle and the common tangent. The order is always 
    [x1,y1,x2,y2] So, the two first values are de coordinates of 
    the intersection with c1 and the last two are for the
     intersection with c2
     There are four diferent cases:
     Both circles are apart each other, then, there are four common tangents.
     The circles intersect in a single point, there are three common tangent
     one of them pass both circles for their common point.
     the circles intersect in two  points there are only two tangents.
     One circle is inside the other and intersec in a point ther is a single
     tangent which pass both circles for their intersection point
     
     One circle inside the other there are not tangents at all
     return an array with the tangents obtained, 
     
    """ 
    
    
    #first we calculate the distance between the circle center
    
    difc2c1=  c2 - c1
    d = np.linalg.linalg.norm(difc2c1[:-1])
    
        
    #and a unit vector in the direction frm c1 center to c2 center, the last 
    #entry is (r2-r1)/d is also employed hereafter
    uc1c2 = difc2c1/d
    
    #Para que existan targentes externas es suficiente que la diferencia de los
    #radios sea menor que la distancia entre los centros...
    
    if d >= np.abs(difc2c1[-1]):
        nplus = -uc1c2[-1] *uc1c2[:-1] +\
        np.array([-1,1])*np.sqrt(1-uc1c2[-1]**2)*uc1c2[-2::-1]
        
        tangents= [[c1[:-1] +\
        c1[-1]*nplus,c2[:-1] + c2[-1]*nplus]]
        
        nless = -uc1c2[-1] *uc1c2[:-1] -\
        np.array([-1,1])*np.sqrt(1-uc1c2[-1]**2)*uc1c2[-2::-1]
        
        tangents.append([c1[:-1] +\
        c1[-1]*nless,c2[:-1] + c2[-1]*nless])
    else:
        tangents = []
        return tangents
        
    #si la distancia es ademas mayor tambien que la suma tenemos el segundo
    #juego de tangentes...
    #Empezamos por sustituir la diferencia de radios por la suma de radios...
    difc2c1[-1] = c1[-1]+c2[-1]
    uc1c2[-1] =  difc2c1[-1]/d    
    if d >= np.abs(difc2c1[-1]):
        
        nplus = uc1c2[-1] *uc1c2[:-1] +\
        np.array([-1,1])*np.sqrt(1-uc1c2[-1]**2)*uc1c2[-2::-1]
        
        tangents.append([c1[:-1] +\
        c1[-1]*nplus,c2[:-1] - c2[-1]*nplus])
        
        nless = uc1c2[-1] *uc1c2[:-1] -\
        np.array([-1,1])*np.sqrt(1-uc1c2[-1]**2)*uc1c2[-2::-1] 
        
        tangents.append([c1[:-1] +\
        c1[-1]*nless,c2[:-1] - c2[-1]*nless])
        
        
    return tangents
        
def pintacirculo(c,ang=[0,2*np.pi],lin='--k',s=[]):
    '''
    Sirve para pintar arcos de circulos entre los angulos ang[0] y ang[1]
    de centro x0,y0 y radio r c=[x0,y0,r]. Por defecto pinta el circulo
    completo en linea discontinua y en negro
        
    '''
    
    theta = np.linspace(ang[0],ang[1],720)
    x = c[0]+c[2]*np.cos(theta)
    y = c[1]+c[2]*np.sin(theta)
    
        
    pl.plot(x,y,lin)
    if s:
     pl.arrow(x[359],y[359],x[361]-x[359],y[361]-y[359],head_width = 1)        
    
def pintactang(c1,c2):
    ''' llama a las anteriores y dibuja los circulos, los puntos de corte
    con las tangentes y los segmentos tangentes...
    Su función es ver que aspecto gereral tienen dos circulos y sus tangentes
    pero no se usa para nada mas
    '''
    tangentes = tancircles(c1,c2) 
    pintacirculo(c1)
    pintacirculo(c2)

    if tangentes:
        for i in tangentes:
            for j in i:
                pl.plot(j[0],j[1],'o')
            pl.plot([i[0][0],i[1][0]],[i[0][1],i[1][1]])
    
    
    
def secnwayp(l):
    ''' 
    l is a list of primary way points the ship should follow in order. 
    The function
    expects each waypoint as an array with elements [x,y,r,s] x and y are the
    coordinates of the waypoint, r is the radius of the circle to turn around
    the primary waypuint, s=+1 rotate clockwise s= -1 rotate counter-clockwise
    when turning around the primary waypoint.
    
    The funtion check if  the first and last primary waypoints have radii
    grater than 0. If not the trajetory described by the list of primary
    waypoints is open: it begins in the first waypoint and ends in the last one
    otherwise, the trajectory is considered close: the last waypoint is
    connected with the first one.

    Except for the case of the first and last primary waypoints in an open
    trajectory, each primary waypoint (pwp) can be associated with two secudary 
    waypoints (swp1 and swp2). swp1 is the common point between the circle 
    around the pwp and the common tangent with the previous pwp. swp2 is the
    common point between the circle of the pwp and the common tangent with the 
    next pwp. In the case first and last pwps of an open trayectory, there are
    a single swp associated and it is taken as swp = pwp.

    The function returns a list: lsec. Each entry of lsec has the following
    structure [pwp,spw1,spw2] pwp is e¡the original primary way point
    taken from the list l, spw1 is an array of coordenates x,y of spw1 
    spw2 is an array of coordenates x,y of spw2. 
    In the case of an open trayectory the first an last entries of lsec only
    two elements: [pwp,spw]         
    '''

    lsec  = [[i] for i in l]
    for n, i in enumerate(zip(l[:-1],l[1:])):
        
        tangentes = tancircles(i[0][:-1],i[1][:-1])
        print tangentes, '\n'
        if (i[0][3] == 1) and (i[1][3] == 1): #clockwise-clockwise
            try:               
                lsec[n].append(tangentes[0][0])
                lsec[n+1].append(tangentes[0][1])
                
            except:
                print 'infeasible trajectory between: ', i
                return lsec
        elif (i[0][3] == -1) and (i[1][3] == -1): #conterclckws-conterclkws
            try:       
                lsec[n].append(tangentes[1][0])
                lsec[n+1].append(tangentes[1][1])

            except:
                print 'infeasible trajectory between: ', i
                return lsec
        elif (i[0][3] == 1) and (i[1][3] == -1): #clckws-conterclkws
            try:       
                lsec[n].append(tangentes[2][0])
                lsec[n+1].append(tangentes[2][1])

            except:
                print 'infeasible trajectory between: ', i
                return lsec
        elif (i[0][3] == -1) and (i[1][3] == 1): #conterclckws-clkws
            try:       
                lsec[n].append(tangentes[3][0])
                lsec[n+1].append(tangentes[3][1])

            except:
                print 'infeasible trajectory between: ', i
                return lsec
    #only applied to close trajectories...            
    if (lsec[0][0][2] > 0.) and (lsec[-1][0][2]>0.):
        tangentes = tancircles(lsec[-1][0][:-1],lsec[0][0][:-1])
        if (lsec[-1][0][3] == 1) and (lsec[0][0][3] == 1):
        #clockwise-clockwise
            try:               
                lsec[-1].append(tangentes[0][0])
                lsec[0].insert(1,tangentes[0][1])
                
            except:
                print 'infeasible trajectory between: ', i
                return lsec
        elif (lsec[-1][0][3] == -1) and (lsec[0][0][3] == -1):
        #conterclckws-conterclkws
            try:       
                lsec[-1].append(tangentes[1][0])
                lsec[0].insert(1,tangentes[1][1])

            except:
                print 'infeasible trajectory between: ', i
                return lsec
        elif (lsec[-1][0][3] == 1) and (lsec[0][0][3] == -1):
        #clckws-conterclkws
            try:       
                lsec[-1].append(tangentes[2][0])
                lsec[0].insert(1,tangentes[2][1])

            except:
                print 'infeasible trajectory between: ', i
                return lsec
        elif (lsec[-1][0][3] == -1) and (lsec[0][0][3] == 1):
        #conterclckws-clkws
            try:       
                lsec[-1].append(tangentes[3][0])
                lsec[0].insert(1,tangentes[3][1])

            except:
                print 'infeasible trajectory between: ', i
                return lsec
    return lsec

def pintadubing(lsec,color = 'b'):
    '''
    This function draws a complete dubing trajectory, departing for a list
    of primary and secundary waypoint such as those generated by secnwayp()
    waypoints and trayectory are draw in blue
    circles are draw in dashed black 
    '''
    
    currutaca = lsec[0][1]
    for i in lsec:
        pintacirculo(i[0]) #Circles involved in Dubbing trajectory
        for j in i:
            pl.plot(j[0],j[1],'o'+color) #waypoint 
         
        pl.plot([currutaca[0],i[1][0]],[currutaca[1],i[1][1]],color)
        
        print i
        #draw the arc covered is radius > 0...            
        if i[0][2]>0:
            print i, 'paso'
            currutaca = i[2] 
            ang = [np.arctan2(i[1][1]-i[0][1],i[1][0]-i[0][0])\
            ,np.arctan2(i[2][1]-i[0][1],i[2][0]-i[0][0])]
            print ang
            ang[1] = ang[1] - i[0][3]*((ang[1]>ang[0])*(i[0][3]>0.)+\
            (ang[1]<ang[0])*(i[0][3]<0.))*2*np.pi
        
            print ang, '\n'        
            pintacirculo(i[0],ang,color,1)
        
    if (lsec[0][0][2] > 0.) and (lsec[-1][0][2]>0.):
        pl.plot([lsec[-1][2][0],lsec[0][1][0]],\
        [lsec[-1][2][1],lsec[0][1][1]],color)

