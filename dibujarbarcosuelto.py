# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 14:19:58 2015

@author: juan
"""
from matplotlib.path import Path
import matplotlib.patches as patches
import barcosolo
b1 = barcosolo.barco()
b1.theta = pi/6
hold(True)
arrow(0,0,1,1,color = 'r')
arrow(1.5,-2,0,0.5,color = 'r')
plot([0,2],[0,0],'k')
#plot([0,0],[-2,2],'k')
plot([-2.1,2],[-2.1,-2.1],'k')
plot([-2,-2],[-2.2,2],'k')
arrow(-2,-2.1,2,2.1)
arrow(0,0,0.5,1,color = 'g')
plot([-2*cos(b1.theta),2*cos(b1.theta)],[-2*sin(b1.theta),2*sin(b1.theta)],'b')
plot([2*sin(b1.theta),-2*sin(b1.theta)],[-2*cos(b1.theta),2*cos(b1.theta)],'b')
plot(b1.pb[0],b1.pb[1],'+')

angulillo =arange(0,b1.theta,b1.theta/100)
plot(cos(angulillo),sin(angulillo),'b')

vertices = np.array([[-1.,-0.25],[-1.,0.25],[-0.25,0.35],[1,0],\
[-0.25,-0.35],[-1.,-0.25]])
rot = array([[np.cos(b1.theta),- sin(b1.theta)],[sin(b1.theta),\
cos(b1.theta)]])  
vertrot = array([np.dot(rot,j) for j in vertices]) +\
[b1.pb[0],b1.pb[1]]       
codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
Path.CURVE3]
pathi = Path(vertrot,codes)
patchi = patches.PathPatch(pathi,facecolor = 'yellow')


gca().add_patch(patchi)


#plot([b1.pb[0]-b1.ls/2.0 * cos(b1.theta),b1.pb[0]-b1.ls/2.0 * cos(b1.theta)-cos(-pi/16)],\
#[b1.pb[0]-b1.ls/2.0*sin(b1.theta),b1.pb[0]-b1.ls/2.0*sin(b1.theta)-sin(-pi/16)],'r')
#
#plot(b1.pb[0]-b1.ls/2.0 * cos(b1.theta)-cos(-pi/16)/2.0,\
#b1.pb[0]-b1.ls/2.0*sin(b1.theta)-sin(-pi/16)/2.0,'ro')
#
#plot([b1.pb[0]-b1.ls/2.0 * cos(b1.theta),b1.pb[0]-b1.ls/2.0 * cos(b1.theta)-cos(pi/3)],\
#[b1.pb[0]-b1.ls/2.0*sin(b1.theta),b1.pb[0]-b1.ls/2.0*sin(b1.theta)-sin(pi/3)],'b')
#
#plot(b1.pb[0]-b1.ls/2.0 * cos(b1.theta)-cos(pi/3)/2.0,\
#b1.pb[0]-b1.ls/2.0*sin(b1.theta)-sin(pi/3)/2.0,'bo')
#
#plot([b1.pb[0]-b1.ls/2.0 * cos(b1.theta)-cos(-pi/16),\
#b1.pb[0]-b1.ls/2.0 * cos(b1.theta)-cos(-pi/16)- cos(-pi/2)],\
#[b1.pb[0]-b1.ls/2.0*sin(b1.theta)-sin(-pi/16),\
#b1.pb[0]-b1.ls/2.0*sin(b1.theta)-sin(-pi/16)-sin(-pi/2)],'r')
#
#plot(b1.pb[0]-b1.ls/2.0 * cos(b1.theta)-cos(-pi/16)- cos(-pi/2)/2.0,\
#b1.pb[0]-b1.ls/2.0*sin(b1.theta)-sin(-pi/16)-sin(-pi/2)/2.0,'ro')
#
#plot([b1.pb[0]-b1.ls/2.0 * cos(b1.theta)-cos(pi/3),\
#b1.pb[0]-b1.ls/2.0 * cos(b1.theta)-cos(pi/3)-cos(3*pi/4)],\
#[b1.pb[0]-b1.ls/2.0*sin(b1.theta)-sin(pi/3),\
#b1.pb[0]-b1.ls/2.0*sin(b1.theta)-sin(pi/3)-sin(3*pi/4)],'b')
#
#plot(b1.pb[0]-b1.ls/2.0 * cos(b1.theta)-cos(pi/3)-cos(3*pi/4)/2.0,\
#b1.pb[0]-b1.ls/2.0*sin(b1.theta)-sin(pi/3)-sin(3*pi/4)/2.0,'bo')

axis('equal')

#Ejes tierra y ejes barco
text(-2.12,1.9,r'$y$')
text(1.9,-2.22,r'$x$')
text(1.9*cos(b1.theta),1.8*sin(b1.theta),r'$x_b$',rotation=b1.theta*180/pi)
text(-2.2*sin(b1.theta),1.9*cos(b1.theta),r'$y_b$',rotation=b1.theta*180/pi)
text(0.1,0.6,r'$\vec{v}_b$',rotation=180/pi*arctan(2))
text(-1.2,-1.3,r'$\vec{r}_{cm}$',rotation=45)
text(1.05*cos(pi/15),1.05*sin(pi/15),r'$\theta_b$')
axis('off')

#anotaciones tensiones etc..
#text(-1.5,-0.35,r'$i-1, \vec{T}_{i-1,b}$')
#text(-1.1,-1.0,r'$i, \vec{T}_{b,i}$') 
#text(-1.8,0.25,r'$i-2$')
#text(-1.0,-1.7,r'$i+1$')
text(0.5,0.8,r'$\vec{F}_m$',rotation=45)
text(1.52,-1.6,r'$\vec{M}_b$')
text(1.52, -2,r'$z$')