# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 12:11:00 2017
Para usar unto a autopsia para sacar los post-mortem de los ensayos
@author: juan
"""
import numpy as np
from numpy.linalg.linalg import norm 
from matplotlib import pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches

import dubing
import bezier_cvr
import corrientes

def dibujarpf(bdcha,cadena,pasos=1,ppo=0, fin = np.inf,multi = False):
    '''dibuja a partir de datos recogidos en arrays de datos tipo barco y cadena
    ver el sistema de preparar matrices para guardar datos dibuja también
    un buque fondeado buque = [eslora, manga, posición, orientación,de     sup'''
    
    
    
    if fin > bdcha.shape[0]-1:
        fin = bdcha.shape[0]-1
        
    for i in range(ppo,fin,pasos):
        if multi:
            figure()
        
        cms = np.array([cadena[i,0,:],cadena[i,1,:]])
        para = np.array([cadena[i,-2,:],cadena[i,-1,:]])
        pl.plot(cms[0,:],cms[1,:],'bo')
        pl.hold(True)
        barrasi = cms + cadena[-1,1,0] * para
        barrasd = cms - cadena[-1,1,0] * para
        pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],'k')
                        

        
        

        vertices = np.array([[-bdcha[-1,6]/2.,-0.25*bdcha[-1,6]/2],\
        [-bdcha[-1,6]/2.,0.25*bdcha[-1,6]/2],\
        [-0.25*bdcha[-1,6]/2,0.35*bdcha[-1,6]/2],[bdcha[-1,6]/2.,0],\
        [-0.25*bdcha[-1,6]/2,-0.35*bdcha[-1,6]/2],[-bdcha[-1,6]/2.,\
        -0.25*bdcha[-1,6]/2]])        
        rot = np.array([[np.cos(bdcha[i,6]),- np.sin(bdcha[i,6])],[np.sin(bdcha[i,6]),\
        np.cos(bdcha[i,6])]])       
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bdcha[i,0],bdcha[i,1]]
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
     
        pathd = Path(vertrot,codes)
        patchd = patches.PathPatch(pathd,facecolor = 'green') #'green'
        pl.gca().add_patch(patchd)
         #######################dibujar cable de arrastre derecha#######################
        rot = np.array([[np.cos(bdcha[i,6]),- np.sin(bdcha[i,6])],[np.sin(bdcha[i,6]),\
        np.cos(bdcha[i,6])]])  
        popad =  np.dot(rot, np.array([-bdcha[-1,6]/2.,0])) + [bdcha[i,0],bdcha[i,1]]  
        tipd = - para[:,-1] * cadena[-1,1,0] + cms[:,-1]
        
        
        distd = norm(popad - tipd)
        dd = distd/cadena[-1,3,1]
        #print dd
        if dd > 1: dd = 1
        r = bezier_cvr.bezier4p([[tipd[0]],[tipd[1]]],[[popad[0]],[popad[1]]],1,1,1.5,\
        (1-dd) * bdcha[i,6]\
         +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
        (1-dd) * np.arctan2(-para[0,0],-para[0,1])\
         +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
        100)
        bezier_cvr.pintar_bezier(r[0],color = 'b')
        ###############################################################################       

