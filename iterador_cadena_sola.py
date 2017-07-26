# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:08:57 2014
Itera una cadena sometida a fuerzas en sus extremos... solo sirve para hacer
una demo de un escenario elemental de arrastre de la cadena

reviewed 08.03.2017
@author: juan
"""
import usv
import boom
#import barchain_matrix
import numpy as np
from matplotlib import pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
rcParams.update({'font.size': 26})
delta = 0.0001 #0.001
#creamos la cadena de momento pequeñita porque es para probar
cd5 = boom.cadena(6)#(5)
cd5.Fd = [250,1200]
cd5.Fi = [-35,-850]
#situamos la cadena ... toda en el eje x con los el primer y último eslabon.
#apuntando para arriba.
cd5.m = 6 #25.
#cd5.s = 0.
#cd5.q = 0.
#cd5.A = 0.
cd5.normal[:,0] = [1,0]
cd5.normal[:,0] = [1,0]
cd5.normal[:,-1] = [-1,0]
cd5.para[:,0] = [0,1]
cd5.para[:,-1] = [0,-1]
cd5.calcms()


tam = 2500#1000 #aqui definimos de momento el tamaño del experimento...
step =2500#165 #10
# y aqui cada cuantas iteraciones guardamos datos y pintamos
#M, B, T = barchain_matrix.Matriz_t_ini(cd5,0)
tiempo = 0
for i in range(tam):
    #print i, i%step
#    M, B, T = barchain_matrix.Matriz_t_ini(cd5,0)
#    M, B, T = barchain_matrix.Matriz_t(M, B, T, cd5)
#    cd5.T[:,0] = T[:-3:2]
#    cd5.T[:,1] = T[1:-2:2]
#    cd5.T[:,2] = T[2:-1:2]
#    cd5.T[:,3] = T[3::2]
    cd5.mtr_s()    
    cd5.movimiento(delta = delta)
    tiempo = tiempo + delta   
    if (i+1)%step == 0:
        print tiempo
        cd5.dibujar('-.k')