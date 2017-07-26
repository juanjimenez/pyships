# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 16:28:49 2016
Para dibujar resultados de experimentos de barcos guardados en archivos
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

def dibujar(bizq,bdcha,cadena):
    '''dibuja a partir de datos recogidos en arrays de datos tipo barco y cadena
    ver el sistema de preparar matrices para guardar datos'''
    pl.figure(1)
    pl.hold(True)
    for i in range(bizq.shape[0]-1):    
        
        cms = np.array([cadena[i,0,:],cadena[i,1,:]])
        para = np.array([cadena[i,-2,:],cadena[i,-1,:]])
        pl.plot(cms[0,:],cms[1,:],'o')
        pl.hold(True)
        barrasi = cms + cadena[-1,1,0] * para
        barrasd = cms - cadena[-1,1,0] * para
        pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],'k')
        
        pl.plot(bizq[i,0],bizq[i,1],'+r') #r
        pl.plot(bdcha[i,0],bdcha[i,1],'+g') #g
#        pl.plot([bizq[i,0],bdcha[i,0]],[bizq[i,1],bdcha[i,1]])
#        normal = bizq[i,1] - bdcha[i,1], - bizq[i,0] + bdcha[i,0]
#        pl.plot([(bizq[i,0]+bdcha[i,0])/2,(bizq[i,0]+bdcha[i,0])/2+normal[0]],\
#        [(bizq[i,1]+bdcha[i,1])/2,(bizq[i,1]+bdcha[i,1])/2 + normal[1]])
#       revisa esto: pinta el barco separado al alargar o acortar la eslora        
#        vertices = np.array([[-1.,-0.25],[-1.,0.25],[-0.25,0.35],[4,0],\
#        [-0.25,-0.35],[-1.,-0.25]])
        
        vertices = np.array([[-bizq[-1,6]/2.,-0.25],[-bizq[-1,6]/2.,0.25],\
        [-0.25,0.35],[bizq[-1,6]/2.,0],[-0.25,-0.35],[-bizq[-1,6]/2.,-0.25]])
        
        rot = np.array([[np.cos(bizq[i,6]),- np.sin(bizq[i,6])],[np.sin(bizq[i,6]),\
        np.cos(bizq[i,6])]])       
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bizq[i,0],bizq[i,1]]       
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
       
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = 'red') #'red'
        pl.gca().add_patch(patchi)

        vertices = np.array([[-bdcha[-1,6]/2.,-0.25],[-bdcha[-1,6]/2.,0.25],\
        [-0.25,0.35],[bdcha[-1,6]/2.,0],[-0.25,-0.35],[-bdcha[-1,6]/2.,-0.25]])        
        rot = np.array([[np.cos(bdcha[i,6]),- np.sin(bdcha[i,6])],[np.sin(bdcha[i,6]),\
        np.cos(bdcha[i,6])]])       
        vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bdcha[i,0],bdcha[i,1]]
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
     
        pathd = Path(vertrot,codes)
        patchd = patches.PathPatch(pathd,facecolor = 'green') #'green'
        pl.gca().add_patch(patchd)
        
        ####################Dibujar cable de arrastre################################
        rot = np.array([[np.cos(bizq[i,6]),- np.sin(bizq[i,6])],[np.sin(bizq[i,6]),\
        np.cos(bizq[i,6])]])  
        popai =  np.dot(rot, np.array([-bizq[-1,6]/2.,0])) + [bizq[i,0],bizq[i,1]]  
        tipi = para[:,0] * cadena[-1,1,0] + cms[:,0]
        
        
        disti = norm(popai - tipi)
        di = disti/cadena[-1,3,0]
        print di
        if di > 1: di = 1
        r = bezier_cvr.bezier4p([[tipi[0]],[tipi[1]]],[[popai[0]],[popai[1]]],1,1,1.5,\
        (1-di) * bizq[i,6]\
         +di * np.arctan2(popai[1]-tipi[1],popai[0] - tipi[0]),\
        (1-di) * np.arctan2(-para[0,0],-para[0,1])\
         +di * np.arctan2(popai[1]-tipi[1],popai[0] - tipi[0]),\
        100)
        bezier_cvr.pintar_bezier(r[0],color = 'b')
        #######################dibujar cable de arrastre derecha#######################
        rot = np.array([[np.cos(bdcha[i,6]),- np.sin(bdcha[i,6])],[np.sin(bdcha[i,6]),\
        np.cos(bdcha[i,6])]])  
        popad =  np.dot(rot, np.array([-bdcha[-1,6]/2.,0])) + [bdcha[i,0],bdcha[i,1]]  
        tipd = - para[:,-1] * cadena[-1,1,0] + cms[:,-1]
        
        
        distd = norm(popad - tipd)
        dd = distd/cadena[-1,3,0]
        print di
        if dd > 1: dd = 1
        r = bezier_cvr.bezier4p([[tipd[0]],[tipd[1]]],[[popad[0]],[popad[1]]],1,1,1.5,\
        (1-dd) * bdcha[i,6]\
         +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
        (1-dd) * np.arctan2(-para[0,0],-para[0,1])\
         +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
        100)
        bezier_cvr.pintar_bezier(r[0],color = 'b')
###############################################################################
        pl.plot(bizq[i,0],bizq[i,1],'+r')
        pl.plot(bdcha[i,0],bdcha[i,1],'+g')
        
        #pl.pause(0.01)        
def graficos(nombref,tipo = 1,pasos=1,ppo=0,fin=np.inf):
    #manejamos el fichero como el culo pero de momento vale.....
    
    if tipo == 1: #two towing ships   
        data = np.load(nombref)
        birec = data['bi']
        bdrec = data['bd']
        cadrec = data['cadena']            
        dibujar(birec,bdrec,cadrec)
        ######para pintar  corrientes (incluido el 13.04.2016)#########################
        #corrientes.vercor(lim,[20,20],corrientes.camp)
        
        ###############################################################################
        pl.figure(2)
        a = []
        for i in zip(bdrec[:-1,2],bdrec[:-1,3]):
            a.append(norm(i))
            
        pl.plot(a,'k')
            
        
        a = []
        for i in zip(birec[:-1,2],birec[:-1,3]):
            a.append(norm(i))
            
        pl.plot(a,'k')
        pl.title('modulo de la velocidad')
        
        a = []
        
        for i in zip(bdrec[:-1,-2],bdrec[:-1,-1]):
            a.append(norm(i)) 
        
        pl.figure(3)
        pl.plot(a,'k')
        
        a = []    
        for i in zip(birec[:-1,-2],birec[:-1,-1]):
            a.append(norm(i)) 
        
        pl.plot(a,'k')
        
        pl.title('Tension en los extremos')
            
        pl.figure(4)
        pl.plot(np.ones(bdrec.shape[0])*0.0)
        pl.plot(bdrec[:-1,6],'k')
        pl.plot(birec[:-1,6],'k') 
        pl.title('rumbo')   
        
        pl.figure(5)
        a = []
        for i in zip(birec[:-1,0],birec[:-1,1],bdrec[:-1,0],bdrec[:-1,1]):
            a.append(norm([i[2]-i[0],i[3]-i[1]]))
            
        pl.plot(a,'k')
        pl.title('distancia entre b')
    else: #a single towing ship, one boom tip tied to the dock and a tanker
        data = np.load(nombref)        
        bdrec = data['bd']
        bdrec[-1,12]= 0 #OJO ANADIDO POR PROBLEMA CON CORRIENTE...en un ejemplo
        cadrec = data['cadena']
        #buque = data['bq']
        buque =np.array([ 86.,16.,-48.,8,-3.14159265])
        if type(buque) is list:
            buque = np.array(buque)
        plan = data['pln']
        try:
            lim =  (data['limites']).tolist()
        except:
            lim = []
         
        
        dibujarpf(bdrec,cadrec,buque,pasos,ppo,fin)
        ######para pintar  corrientes (incluido el 13.04.2016)################
        if lim:
            corrientes.vercor(lim,[20,20],\
            fun = corrientes.campos[int(bdrec[-1,12])])
        ######################################################################
        try:
                        
            a =np.array([norm(i) for i in zip(bdrec[:-1,2],bdrec[:-1,3])])
            offset =\
            np.concatenate(([0],np.cumsum([i[2][-1] for i in plan[:-1]])))
            tmp = np.concatenate([i[2]+k for i,k in zip(plan,offset)])
            vel = np. concatenate([norm(i[1],axis = 0) for i in plan])
            pl.figure(1)
            for i in plan:
                bezier_cvr.pintar_bezier(i[0])
                
            pl.figure(2)                
            pl.plot(tmp,vel)
            pl.plot(tmp,a,'k')
            pl.title('modulo de la velocidad')
        except:
            #para hacerlo compatible con archivos anteriores a may_27_2016
            print 'old file or dubins planning caca'
            
            try:
                for i in plan:
                    bezier_cvr.pintar_bezier(i)
                print('por aqui')    
            except:
                dubing.pintadubing(plan,'k')
                print('por alla')
                pl.figure(5)
                #dubing.pintadubing(plan)
            
            pl.figure(2)
            
            a = []
            
            for i in zip(bdrec[:-1,2],bdrec[:-1,3]):
                a.append(norm(i))
                
            time = np.arange(len(a))
            
            time = 2* time
            pl.plot(time, a,'k')    
            #pl.title('modulo de la velocidad')
#        a =np.array([norm(i) for i in zip(bdrec[:-1,2],bdrec[:-1,3])])
##for i in zip(bdrec[:-1,2],bdrec[:-1,3]):
##    a.append(norm(i))    
#        offset = np.concatenate(([0],np.cumsum([i[2][-1] for i in plan[:-1]])))
#        tmp = np.concatenate([i[2]+k for i,k in zip(plan,offset)])
#        vel = np. concatenate([norm(i[1],axis = 0) for i in plan])
#        pl.plot(tmp,vel)
#        pl.plot(tmp,a,'g')
##for i in plan:
##    pl.plot(i[2] + t0,norm(i[1],axis = 0))
##    t0 = i[2][-1]    
#
#        pl.title('modulo de la velocidad')
        
        a = []
        
        for i in zip(cadrec[:-1,8,0],cadrec[:-1,9,0]):
            a.append(norm(i)) 
        
        pl.figure(3)
        #pl.plot(a,'r')
        
        a = []    
        for i in zip(bdrec[:-1,-1],bdrec[:-1,-2]):
            a.append(norm(i)) 
        print len(a)
        print time.size
        pl.plot(time,a,'k')
        
        #pl.title('Tension en los extremos')
            
        pl.figure(4)
        pl.plot(np.ones(bdrec.shape[0])*0.0)
        pl.plot(bdrec[:-1,6],'k')
        pl.title('rumbo')   
        
#        pl.figure(5)
#        if buque.size:
#            vertices = np.array([[-buque[0]/2.,-buque[1]/2.],\
#            [-buque[0]/2.,buque[1]/2],\
#            [buque[0]/4.,buque[1]/2],\
#            [buque[0]/2,0],\
#            [buque[0]/4.,-buque[1]/2],\
#            [-buque[0]/2.,-buque[1]/2.]])        
#            rot = np.array([[np.cos(buque[4]),- np.sin(buque[4])],[np.sin(buque[4]),\
#            np.cos(buque[4])]])
#            
#            vertrot = np.array([np.dot(rot,j) for j in vertices]) \
#            + [buque[2],buque[3]]
#                               
#            codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
#            Path.CURVE3]
#           
#            pathi = Path(vertrot,codes)
#            patchi = patches.PathPatch(pathi,facecolor = 'black') #'red'
#            pl.gca().add_patch(patchi)
#            
#            pl.plot([5,np.min(vertrot[:,0]) - 5],[0,0],'k')
#        pl.axis('equal')
#        cms = np.array([cadrec[-2,0,:],cadrec[-2,1,:]])
#        para = np.array([cadrec[-2,-2,:],cadrec[-2,-1,:]])
#        pl.plot(cms[0,:],cms[1,:],'bo')
#        pl.hold(True)
#        barrasi = cms + cadrec[-1,1,0] * para
#        barrasd = cms - cadrec[-1,1,0] * para
#        pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],'b')

        #a = []
        #for i in zip(birec[:-1,0],birec[:-1,1],bdrec[:-1,0],bdrec[:-1,1]):
        #    a.append(norm([i[2]-i[0],i[3]-i[1]]))
        #    
        #pl.plot(a,'k')
        #pl.title('distancia entre b')   
        
def dibujarpf(bdcha,cadena,buque = np.array([10.,2.,0.,0.,0.]),pasos=1,ppo=0, fin = np.inf):
    '''dibuja a partir de datos recogidos en arrays de datos tipo barco y cadena
    ver el sistema de preparar matrices para guardar datos dibuja también
    un buque fondeado buque = [eslora, manga, posición, orientación,de     sup'''
    pl.figure(1)
    pl.hold(True)
    if buque.size:
        vertices = np.array([[-buque[0]/2.,-buque[1]/2.],\
        [-buque[0]/2.,buque[1]/2],\
        [buque[0]/4.,buque[1]/2],\
        [buque[0]/2,0],\
        [buque[0]/4.,-buque[1]/2],\
        [-buque[0]/2.,-buque[1]/2.]])        
        rot = np.array([[np.cos(buque[4]),- np.sin(buque[4])],[np.sin(buque[4]),\
        np.cos(buque[4])]])
        
        vertrot = np.array([np.dot(rot,j) for j in vertices]) \
        + [buque[2],buque[3]]
                           
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
       
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = 'green') #'red'
        pl.gca().add_patch(patchi)
        pl.plot([5,np.min(vertrot[:,0]) - 5],[0,0],'m',linewidth = 5)
    pl.axis('equal')
    
    if fin > bdcha.shape[0]-1:
        fin = bdcha.shape[0]-1
        
    for i in range(ppo,fin,pasos):    
        
        cms = np.array([cadena[i,0,:],cadena[i,1,:]])
        para = np.array([cadena[i,-2,:],cadena[i,-1,:]])
        pl.plot(cms[0,:],cms[1,:],'b.')
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
        patchd = patches.PathPatch(pathd,facecolor = 'red') #'green'
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

        

    pl.plot(bdcha[i,0],bdcha[i,1],'+g')
#    pl.figure(5)
#    pl.hold(True)
#    if buque.size:
#        vertices = np.array([[-buque[0]/2.,-buque[1]/2.],\
#        [-buque[0]/2.,buque[1]/2],\
#        [buque[0]/4.,buque[1]/2],\
#        [buque[0]/2,0],\
#        [buque[0]/4.,-buque[1]/2],\
#        [-buque[0]/2.,-buque[1]/2.]])        
#        rot = np.array([[np.cos(buque[4]),- np.sin(buque[4])],[np.sin(buque[4]),\
#        np.cos(buque[4])]])
#        
#        vertrot = np.array([np.dot(rot,j) for j in vertices]) \
#        + [buque[2],buque[3]]
#                           
#        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
#        Path.CURVE3]
#       
#        pathi = Path(vertrot,codes)
#        patchi = patches.PathPatch(pathi,facecolor = 'black') #'red'
#        pl.gca().add_patch(patchi)
#        pl.plot([5,np.min(vertrot[:,0]) - 5],[0,0],'k')
#    pl.axis('equal')
#    cms = np.array([cadena[i,0,:],cadena[i,1,:]])
#    para = np.array([cadena[i,-2,:],cadena[i,-1,:]])
#    pl.plot(cms[0,:],cms[1,:],'bo')
#    pl.hold(True)
#    barrasi = cms + cadena[-1,1,0] * para
#    barrasd = cms - cadena[-1,1,0] * para
#    pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],'k')
#                    

    
    

#    vertices = np.array([[-bdcha[-1,6]/2.,-0.25*bdcha[-1,6]/2],\
#    [-bdcha[-1,6]/2.,0.25*bdcha[-1,6]/2],\
#    [-0.25*bdcha[-1,6]/2,0.35*bdcha[-1,6]/2],[bdcha[-1,6]/2.,0],\
#    [-0.25*bdcha[-1,6]/2,-0.35*bdcha[-1,6]/2],[-bdcha[-1,6]/2.,\
#    -0.25*bdcha[-1,6]/2]])      
#    vertrot = np.array([np.dot(rot,j) for j in vertices]) + [bdcha[i,0],bdcha[i,1]]
#    codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
#    Path.CURVE3]
# 
#    pathd = Path(vertrot,codes)
#    patchd = patches.PathPatch(pathd,facecolor = 'green') #'green'
#    pl.gca().add_patch(patchd)
#     #######################dibujar cable de arrastre derecha#######################
#    rot = np.array([[np.cos(bdcha[i,6]),- np.sin(bdcha[i,6])],[np.sin(bdcha[i,6]),\
#    np.cos(bdcha[i,6])]])  
#    popad =  np.dot(rot, np.array([-bdcha[-1,6]/2.,0])) + [bdcha[i,0],bdcha[i,1]]  
#    tipd = - para[:,-1] * cadena[-1,1,0] + cms[:,-1]
#    
#    
#    distd = norm(popad - tipd)
#    dd = distd/cadena[-1,3,1]
#    #print dd
#    if dd > 1: dd = 1
#    r = bezier_cvr.bezier4p([[tipd[0]],[tipd[1]]],[[popad[0]],[popad[1]]],1,1,1.5,\
#    (1-dd) * bdcha[i,6]\
#     +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
#    (1-dd) * np.arctan2(-para[0,0],-para[0,1])\
#     +dd * np.arctan2(popad[1]-tipd[1],popad[0] - tipd[0]),\
#    100)
#    bezier_cvr.pintar_bezier(r[0],color = 'b')    
#    #pl.pause(0.01)        
#        #hold(False)
