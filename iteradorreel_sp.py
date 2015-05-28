# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 18:33:51 2015
Este programa itera el movimiento de un boom desde un reel se supone 
que el reel esta situado a la 'derecha' del boom y el barco que tira de el a la
izquierda... (en proyecto todavia)
@author: juan
"""

import numpy as np
import numpy as np
from numpy.linalg.linalg import norm 
import matplotlib.pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
#import bezier_cvr
import boom
import reel_solo
import barcosolo
import barchain_matrix_sp_all

#supongo que ya he importado bastante
delta = 0.01 #paso de integracion, tiempo

#creamos un reel (recordad que en principio suponemos en reel colocado en (0,0)
#tomamos el valor por defecto de la candena contenida en el interior que 
#contiene 10 eslabones... de los cuales dos enteros están fuera del reel
#Esto deja el tip del reel a cero... lo cual es una fuente de marrones...
rl = reel_solo.reel()

#creo una cadena que contiene solo los dos eslabones sueltos 
#del reel
cad = boom.cadena(2) 

#el reel tiene el extremo del boom apuntando hacia abajo, asique engachamos los
#dos elabones por la derecha al extremo del boom enrollado en el reel
#posicion del cm del primer eslabon
cad.cms[:,0] = [cad.L,- rl.r]
#orientacion del primer eslabon

cad.calcms(primero = cad.cms[:,0])

#por ultimo construimos un barquito a la derecha y lo colocamos en su sitio

bd = barcosolo.barco()
bd.theta = 0
bd.pb[:] = [3,-rl.r]

#voy a pintar todo un momento a ver que coño pasa...
rl.dibujar()
cad.dibujar()
bd.dibujar()

#iniciamos la matriz de tensiones
M, B, T = barchain_matrix_sp_all.Matriz_t_ini(cad,reel = 1)

#damos un paso de integracion, just to see what happens
M, B, T = barchain_matrix_sp_all.Matriz_trb(M,B,T,cad,rl,bd)
 



