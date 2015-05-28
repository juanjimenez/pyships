# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 22:43:48 2013
Este fichero define un objeto llamado cadena. El objeto contiene una cadena de
logitud variable, y tiene especificadas una serie de propiedades mecanicas de 
las cadenas.
Ademas, incluye funciones para pintar cadenas, moverlas etc.
@author: juan
"""
import numpy as np
from matplotlib import pyplot as pl

class cadena:
    
    
    def __init__(self,n=10):
        
        self.objeto = 'cadena'
        #definicion de variables de estado de la cadena
        
        #numero de elementos en la cadena        
        self.esl = n
        
        #superficie del elemento de cadena
        self.s = 2.0
        
        #grosor del elemento de cadena
        self.q = 0.1 #hinchado???
        
        #grosor del elemento de cadena deshinchado (solo el util cuando la
        #cadena esta en la bobina)
        self.dh = 0.01 
        
        #resitencia al giro del elemento de cadena
        self.A = 10.0
        
        # mitad de longitud del elemento de cadena
        self.L = 0.5
        
        # masa del elemento de cadena
        self.m = 15.0
        
        #Fuerzas aplicadas al los extremos de la cadena. Fext[0], fuerza aplica
        #da en el extremo izquierdo, Fext[1] fuerza aplicada en el extremo dere
        #cho. Si hay barcos tirando de los extremos estas fuerzas no actuan
        self.Fd = np.zeros(2)
        self.Fi = np.zeros(2)
        #factor ml2  = m*L**2 se emplea mucho en los calculos de la dinámica,
        #así que mejor tenerlo definido de antemano
        self.mL2 = self. m * self.L**2
        
        #momento de inercia del elemento de cadena lo distingo de mL2. Me
        #parece que una cadena un poco más complicada de estructura, puede 
        #tener un momento de inercia mas complejo
        #self.I = 0.5 * self.mL2       
        self.I = self.mL2 / 12
        #vectores normales a los elementos de la cadena
        #por defecto todos apuntan en la direccion 'y'
        self.normal = np.array([np.zeros(n),np.ones(n)]) 
        self.para = np.array([-self.normal[1],self.normal[0]]) #vector paralelo 
        #al eslabon
        
        # posicion del centro de masas de los elementos de las cadenas
        self.calcms()
        
        #velocidad de los elementos de la cadena
        #tocada para evitar nans en la division por el
        #modulo de v (esto necesita una revisión más seria)
        self.v = np.zeros(self.normal.shape)
        #modulos de las velocidades. Inicializados a 1 para evitar divisiones 
        #entre cero. Valdría la pena repasar las ecuaciones y buscar una forma
        #elegante de manejar esto...
        self.vmod = np.ones(self.esl)
        
        # Aceleracion de los elementos de la cadena
        self.a = self.v.copy()        
        
        # perturbaciones debidas a factores externos
        self.per = self.v.copy()
        
        
        #aceleracion angular de los elementos de la cadena
        self.alfa = np.zeros(self.normal.shape[1])
        
        #velocidad angular de los elementos de la cadenabi.movimiento
        self.w = self.alfa.copy()

        #termino independiente para calculo de tensiones
        #self.B = np.zeros(2 * self.esl +2)
        #matriz de coeficientes para el calculo de la tensiones
        #Si tenemos n eslabones tendremos n+1 tension, incluyendo las
        #de los extremos. como cada tensión tiene componentes x e y,
        #en total tenemos 2 * + 2 tensiones que calcular.... de ahí el 
        #tamaño de B y M
        #self.M = np.zeros((2 * self.esl +2, 2 * self.esl +2))
        #vector de tensiones aplicadas a los eslabones
        #hacemos hueco para dos tensiones por eslabon, aun a sabiendas de que
        #la continuidad de los eslabones hace que la tension entre dos eslbones
        #consecutivos sea la misma. LO hacemos así, porque permite intercalar
        #barcos en mitad de la cadena con gran facilidad... 
        self.T = np.zeros([self.esl,4])

    def resetear(self,primero = np.zeros(2)):
        #se resetea una cadena con un numero de eslabones ya definido en el
        #constructor
        self.s = 2.0 #superficie del elemento de cadena
        self.q = 0.1 #grosor del elemento de cadena
        self.A = 1.0 #10.0 #resitencia al giro del elemento de cadena
        self.L = 0.5 # mitad de longitud del elemento de cadena
        self.m = 10.0 # masa del elemento de cadena
        #Fuerzas aplicadas al los extremos de la cadena. Fext[0], fuerza aplica
        #da en el extremo izquierdo, Fext[1] fuerza aplicada en el extremo dere
        #cho. Si hay barcos tirando de los extremos estas fuerzas no actuan
        self.Fd = np.zeros(2)
        self.Fi = np.zeros(2)
        #factor ml2  = m*L**2 se emplea mucho en los calculos de la dinamica,
        #asi que mejor tenerlo definido de antemano
        self.mL2 = self. m * self.L**2
        #momento de inercia del elemento de cadena lo distingo de mL2. Me
        #parece que una cadena un poco más complicada de estructura, puede 
        #tener un momento de inercia mas complejo
        self.I = 0.5 * self.mL2       
        #vectores normales a los elementos de la cadena
        #por defecto todos apuntan en la direccion 'y'
        self.normal = np.array([np.zeros(self.esl),np.ones(self.esl)]) 
        self.para = np.array([-self.normal[1],self.normal[0]]) #vector paralelo
        #al eslabon
        
        # posicion del centro de masas de los elementos de las 
        #cadenas
        self.calcms(primero)
        #velocidad de los elementos de la cadena
        self.v = np.zeros(self.normal.shape)
        #modulos de las velocidades. Inicializados a 1 para evitar divisiones 
        #entre cero. Valdría la pena repasar las ecuaciones y buscar una forma
        #elegante de manejar esto...
        self.vmod = np.ones(self.esl)
        # Aceleracion de los elementos de la cadena
        self.a = self.v.copy()
        # perturbaciones debidas a factores externos. (no incluidas de momento
        #en los cálculos)
        self.per = self.v.copy()
        #aceleracion angular de los elementos de la cadena
        self.alfa = np.zeros(self.normal.shape[1])
        #velocidad angular de los elementos de la cadena
        self.w = self.alfa.copy()
        #termino independiente para calculo de tensiones
        #self.B = np.zeros(2 * self.esl +2)
        #matriz de coeficientes para el calculo de la tensiones
        #Si tenemos n eslabones tendremos n+1 tension, incluyendo las
        #de los extremos. como cada tensión tiene componentes x e y,
        #en total tenemos 2 * + 2 tensiones que calcular.... de ahí el 
        #tamaño de B y M
        #self.M = np.zeros((2 * self.esl +2, 2 * self.esl +2))
        #hacemos hueco para dos tensiones por eslabon, aun a sabiendas de que
        #la continuidad de los eslabones hace que la tension entre dos eslbones
        #consecutivos sea la misma. LO hacemos así, porque permite intercalar
        #barcos en mitad de la cadena con gran facilidad...        
        self.T = np.zeros([self.esl,4])
        
    def listar(self):        
        print "caracteristicas cadena \n"   
        print 'esl', self.esl 
        print 's', self.s
        print 'q', self.q
        print 'A', self.A
        print 'L', self.L
        print 'm', self.m
        print 'Fd', self.Fd
        print 'Fi', self.Fi
        print 'mL2', self.mL2 
        print 'I', self.I       
        print 'normal', self.normal 
        print 'para', self.para        
        print 'cms', self.cms
        print 'v', self.v
        print 'vmod', self.vmod
        print 'a', self.a       
        print 'per', self.per
        print 'alfa', self.alfa
        print 'w', self.w
        print 'T', self.T
        
    
    def calcms(self,primero = np.zeros(2),ultimo = np.zeros(2), orden = 1):
        '''Esta funcion calcula la posicion de los centros de masa de todos
        los elementos de la cadena a partir de la posicion del primero, si la
        variable orden = 1, o a partir del último si la variable
        orden = -1, y de la orientacion del resto, dada por la matriz que 
        contiene los vectores normales de todos los elementos'''
        self.cms = np.zeros(self.normal.shape)
        
        if orden == 1:
            self.cms[0,0] = primero[0]
            self.cms[1,0] = primero[1]
            for i in range(1,self.cms.shape[1]):
                #de primero (izquierda) a ultimo (derecha)
                #calculamos las posiciones de los centros de los elementos,
                # a partir de los vectores normales a ellos, las componentes son 
                #las mismas intercambiando x e y y cambiando x de signo para 
                #obterner un vector unitario en la direccion del elemento.
                #midiendo las direcciones de izquierda a derecha
                #sumamos la posición del elemento anterior. 
                self.cms[0,i] = self.cms[0,i-1] + self.L * self.normal[1,i-1] +\
                self.L *self. normal[1,i]
                self.cms[1,i] = self.cms[1,i-1] - self.L * self.normal[0,i-1] -\
                self.L * self.normal[0,i]
            self.ord = 1

        else:
            self.cms[0,-1] = ultimo[0]
            self.cms[1,-1] = ultimo[1]
            for i in range(self.cms.shape[1]-2,-1,-1):
                
                #de ultimo (derecha) a primero (izquierda)
                #calculamos las posiciones de los centros de los elementos,
                # a partir de los vectores normales a ellos, las componentes son 
                #las mismas intercambiando x e y y cambiando x de signo para 
                #obterner un vector unitario en la direccion del elemento.
                #midiendo las direcciones de izquierda a derecha
                #sumamos la posición del elemento anterior. 
                self.cms[0,i] = self.cms[0,i+1] - self.L * self.normal[1,i+1] -\
                self.L * self.normal[1,i]
                self.cms[1,i] = self.cms[1,i+1] + self.L * self.normal[0,i+1] +\
                self.L * self.normal[0,i]
                
            self.ord = -1
            
    def movimiento(self, delta = 5e-4):            
            ''' a partir de la matriz de coeficientes y del vector de términos
            independientes para el calculo de tensiones en la cadena,
            calculamos el movimiento de las cadenas para un paso de integración 
            Delta. (esto está hecho muy cutre usando Euler, Más adelante vere-
            mos la manera de emplea un par encajado como por ejemplo
            Dorman-price)'''
            
            #obtenemos las tensiones en los extremos de los eslabones:
            #self.T = np.linalg.linalg.solve(M,B)
            
            #diferencia entre tensiones aplicadas alos extremos de cada eslabón
            dT = np.transpose(self.T[:,2:] - self.T[:,0:2])
            sT = np.transpose(self.T[:,2:] + self.T[:,0:2])
             #self.T[1:-2:2]])
            #y a partir de ellas las aceleraciones de cada eslabón            
            self.a = (dT - (np.abs(np.sum(self.v*self.normal,0)) * self.s\
            + np.abs(np.sum(self.v*self.para,0)) * self.q) * self.v/self.vmod)\
            / self.m\
            #-2*np.array([-self.v[1],self.v[0]])*self.w#término coriolis
            
            #y la aceleración angular de los eslabones
            self.alfa =(np.sum(sT * self.normal,0) * self.L - self.A * self.w)\
            /self.I             
            #a partir de aqui y del paso de tiempo delta, avanzamos las
            #velocidades
            
            #primero la velocidad lineal del cm
            self.v += self.a * delta
            self.vmod = np.sqrt(sum(self.v**2,0))
            #odio buscar los que tienen modulo cero y sustituirlos...
            #es lento y sucio y no me gusta
            #hay que buscar una codificación mejor de las ecuaciones..
            self.vmod[np.nonzero(self.vmod == 0)] = 1
            
            #y despues la angular en torno a el...
            self.w += self.alfa * delta 
            
            #modificamos la orientación de los eslabones
            
            cthsl = self.normal * np.cos(self.w * delta)
            sthsl = self.normal * np.sin(self.w * delta)
            
            self.normal[:] = [cthsl[0,:] - sthsl[1,:],
            cthsl[1,:] + sthsl[0,:]]
            self.para[:] = [-self.normal[1,:], self.normal[0,:]]
            
            #y la posicion de los centros de masa
            self.cms += self.v * delta 
            
    def dibujar(self):
           
            pl.plot(self.cms[0,:],self.cms[1,:],'o')
            pl.hold(True)
            barrasi = self.cms + self.L * self.para 
            barrasd = self.cms - self.L * self.para
            pl.plot([barrasi[0,:],barrasd[0,:]],[barrasi[1,:],barrasd[1,:]],'k')
            
            