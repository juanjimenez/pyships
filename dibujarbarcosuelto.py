# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 14:19:58 2015
Dibuja un esquema del USV y de las fuerzas, velocidad etc...
Vale como esquema para las publicaciones
reviewed 8.03.2017
@author: juan
"""
from matplotlib.path import Path
import matplotlib.patches as patches
import usv
b1 = usv.barco(1)
b1.ls = 1
b1.theta = pi/8
b1.pb=[1.,1.]
hold(True)
b1.dibujar(color='grey')
arrow(0,0,1,1,length_includes_head=True,color='k') #,color = 'r')

#angulo orientacion del barco
plot([b1.pb[0],b1.pb[0]+0.7],[b1.pb[1],b1.pb[1]],'--k')
#
angulillo =arange(0,b1.theta/2,b1.theta/100)
plot(0.5*b1.pb[0]+cos(angulillo),b1.pb[1]+sin(angulillo),'k')
text(b1.pb[0]+0.52*cos(b1.theta/2),b1.pb[1]+0.42*sin(b1.theta/2),r'$\theta_b$')
#

##Ejes tierra 
arrow(-0.2,0,2.2,0,length_includes_head=True,color='k')
arrow(0,-0.2,0,2.2,length_includes_head=True,color = 'k')
text(-0.12,1.9,r'$y$')
text(1.9,-0.12,r'$x$')
#ejes barco
arrow(b1.pb[0]-0.07*cos(b1.theta),b1.pb[1]-0.07*sin(b1.theta),\
0.85*cos(b1.theta),0.85*sin(b1.theta),\
length_includes_head=True,color='k')
arrow(b1.pb[0]+0.07*sin(b1.theta),b1.pb[1]-0.07*cos(b1.theta),\
-0.85*sin(b1.theta),0.85*cos(b1.theta),\
length_includes_head=True,color = 'k')
text(b1.pb[0]+0.7*cos(b1.theta),b1.pb[1]+0.5*sin(b1.theta),\
r'$x_b$',rotation=b1.theta*180/pi)
text(b1.pb[0]-1*sin(b1.theta),b1.pb[1]+0.65*cos(b1.theta),\
r'$y_b$',rotation=b1.theta*180/pi)

#Fuera ejercida y velocidad instantanea
arrow(b1.pb[0],b1.pb[1],0.25,0.5)
text(b1.pb[0]+0.05,b1.pb[1]+0.35,r'$\vec{v}_b$',rotation=180/pi*arctan(2))
arrow(b1.pb[0],b1.pb[1],0.1,0.9)
text(b1.pb[0]-0.09,b1.pb[1]+0.45,r'$\vec{F}_m$',rotation=180/pi*arctan(9))

text(0.3*b1.pb[0],0.5*b1.pb[1],r'$\vec{r}_{cm}$',rotation=45)


axis([-0.25,2.25,-0.25,2.25])
grid('on')
gca().axes.get_xaxis().set_ticklabels([])
gca().axes.get_yaxis().set_ticklabels([])
hold(False)