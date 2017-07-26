# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 20:53:38 2015
This script performs a test on USV speed set point values.
It obtains the relationships among values of setl (the variable of class 
barco that stablish the USB speed setpoint), USV stationary speed and time
taken by the USV to reach the speed
The USV controller is cut out and the set-point is feed directly to the engine

Reviewed 06.03.2017 Probably there are some paramaters bad scaled the result
are somewhat strange
@author: juan
"""
# Valores del experimento
import usv
import numpy as np
rumbo = 0.0 #pi/2 #-pi/2 #0.0
velocidad = 4.0
delta = 0.01

#creamos un barquito
bd =usv.barco(3)

bd.ls = 1
bd.Ac = 0.4 
bd.mb = 3.4
bd.Ib = 0.85
bd.mua = 0 #caveat in the equations damping is bd.mua*bd.ls*bd.wb  
bd.mua2 = 9.7577 #caveat in the equations damping is bd.mua*bd.ls*bd.wb**2 
bd.mut2 = 2.226
bd.mut = 0.1
bd.mul2= 1.219
bd.mul = 0.01
bd.krpid = [5.,0.5,0]
bd.kvpid = np.array([5,0.,5])
bd.krd = [3.,0.,0.]
bd.krl = [0.,0.,0.]
bd.thewjmax = 20 * np.pi / 180
bd.pmax = 10
bd.theta = np.pi/2


bd.controlador = lambda a,b,c,d:[]




