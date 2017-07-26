# -*- coding: utf-8 -*-
"""
Created on Sat Jul 15 12:15:08 2017
ESte archivo se creo para obtener los graficos del articulo de las pelotas
@author: Juanfran Jiménez C
"""
import graficos
import numpy as np
from numpy.linalg.linalg import norm 
from matplotlib import pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches

import dubing
import bezier_cvr
import corrientes
#import read_log
def graf(nombref,pasos=1,ppo=0,fin=inf):
 
# wp = [np.array([0.,0.,0,-1]),\
# np.array([-18.654246346914539,1.5,14.,-1]),\
# np.array([-18.654246346914539-14,0,0,-1])]
##dubin secondary waypoints:
# wp2 = dubing.secnwayp(wp)
# dubing.pintadubing(wp2)   
 data = np.load(nombref)        
 bdrec = data['bd']
 cadrec = data['cadena']
 buque = data['bq']
 plan = data['pln']
 dubing.pintadubing(plan,'k')
 pl.figure(5)
 dubing.pintadubing(plan)
 pl.figure(1)
 pl.plot([plan[0][0][0],plan[2][0][0]],[-1,-1],'k',linewidth=4)
 try:
  lim =  (data['limites']).tolist()
 except:
  lim = []
 if lim:
  corrientes.vercor(lim,[20,20],\
  fun = corrientes.campos[int(bdrec[-1,12])]) 
 graficos.dibujarpf(bdrec,cadrec,buque,pasos,ppo,fin)
 return data
 