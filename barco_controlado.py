# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 20:53:38 2015
Esto es una prueba para mover un barco con control de velocidad y rumbo
reviewed and updated: 8.03.2017
@author: juan
"""
# Valores del experimento
import usv
import numpy as np
rumbo =pi/4 #-pi/2 #0.0
velocidad = 1.0
delta = 0.01
tf = 5
paso = 60
#creamos un barquito
barquito = usv.barco(3)

#barquito.pmax = 500
barquito.kvpid = [500,50,50]
barquito.krpid = [4,2,0.0]
#lo hacemos moverse a ver que pasa

usv.trazar(barquito,delta,rumbo,velocidad,tf,paso,1000.)
