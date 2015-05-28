# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 20:53:38 2015
Esto es una prueba para mover un barco con control de velocidad y rumbo
@author: juan
"""
# Valores del experimento
import barcosolo
import numpy as np
rumbo = 0.0
velocidad = 2.0
delta = 0.01
tf = 3600
paso = 100
#creamos un barquito
barquito = barcosolo.barco()

#barquito.pmax = 500
#barquito.kvpid = [1000, 10, 10]
#lo hacemos moverse a ver que pasa

barcosolo.trazar(barquito,delta,rumbo,velocidad,tf,paso)
