# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 11:12:01 2015
ESte programita simula una bobina que contiene una barrera.
@author: juan
"""
import numpy as np
from numpy.linalg.linalg import norm 
import matplotlib.pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
import bezier_cvr
import boom
class reel:
    """Esta clase define un objeto tipo reel, se trata de una bobina que
    contiene una barrera de contención enrollada: un objeto tipo
    cadena, por defecto la cadena enrollada tiene 10 eslabones de los cuales
    hay 2 que se encuentran desenrollados"""
    def __init__(self, cad = boom.cadena(),out = 2):
        #definimos los parametros de la bobina
        
            #un reel sin boom es como un arbol sin hojas...
            #asi que lo primero que hacemos es crear una cadena asociada al reel
            self.cad = cad #cadena por defecto es de 10 eslabones
            self.empty_weight = 100. #peso de la bobina en vacio
            self.M = self.empty_weight + (self.cad.esl - out) * self.cad.m
            #peso de la bobina cargada 
            self.empty_ra = 1. #radio  de la bobina en vacio
            self.r= self.empty_ra / 2. \
            + np.sqrt(self.empty_ra ** 2 + 4. * self.cad.dh * (self.cad.esl - out)\
            * self.cad.L * 2.) / 2. #radio del reel cargado, el boom se supone
                                  #enrollado entorno a empty_ra
                                  #el parametro L de un boom esta definido como 
                                  #la mitad de su longitud... Suponga además
                                  #que empezamos con un cierto numero de eslabones
                                  #fuera del reel (2 por defecto)
            self.I = self.M * self.r ** 2 / 2  #momento de inercia del reel
            self.alpha = 0. #aceleración angular del reel
            self.w = 0. #velocidad angular del reel
            self.theta = 0. #Cuanto angulo lleva girado. Se acumula mientras
                           #quede bobina
            self.Ar = 100. #resistencia al giro del reel
            self.punto = 3 * np.pi/ 2 # angulo que forma el extremo de la
                               #cadena enrollada con el eje x+
            #definimos como parte del reel, el trozo de cadena que sobresale 
            #del reel pero no constituye todavía un eslabon completo...            
            self.lt = 0. # longitud del tip
            self.tht = 0. #orientación del tip
            self.mt = self.cad.m * self.lt / self.cad.L
            self.It = self.mt * self.lt ** 2 / 24 #momento de 
            #Inercia del tip, se ha elegido tipo varilla... 
            self.wt = 0 #velocidad angular del tip
            self.vt = [0.0,0.0]
            self.At = 0 #resistencia al giro del tip
            
            
            #aqui ya esta todo.... lo que hace falta.. ahora vendría definir
            #una dinamica... Meto aquí una pero con poca fe, ya que no tengo
            #ahora mismo demasiado claro si puedo separarla de la del boom...
            
            
    def movimiento(self,fuerza = np.zeros(2),delta = 5e-4):
        #cosideramos el movimiento del reel y el del trocito de eslabon que
        #se está soltando. Que es siempre parte de l objeto reel
        #fuerza, es la uerza ejercidad en la punta del trocito de eslabón
        #que se está soltando.
        #lo primero que hacemos es contruir 
    ######OJO FUNCION EN CONSTRUCCION########################################
        self       
                                            
    def dibujar(self):
        #pintamos un circulo a lo bruto que es el reel vacío
        fs = np.linspace(0,2*np.pi,50)
        pl.plot(self.empty_ra * np.cos(fs),self.empty_ra * np.sin(fs),'r')
        # y un segundo circulo con el reel lleno
        pl.plot(self.r * np.cos(fs),self.r * np.sin(fs))
        #pintamos el punto...
        pl.plot([0,self.r * np.cos(self.punto)],[0,\
        self.r * np.sin(self.punto)])
        #pintamos el tip
        pl.plot([self.r * np.cos(self.punto), self.r * np.cos(self.punto) \
        -self.lt * np.sin(self.punto)],[self.r * np.sin(self.punto),\
        self.r * np.sin(self.punto) + self.lt * np.cos(self.punto)])
                
                
        
    
