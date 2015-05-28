# -*- coding: utf-8 -*-
"""
Created on Thu Sep 05 20:57:04 2013
Barco solo es un ficherito que define las variables y las ecuaciones necesarias 
para simular el movimiento de un barco suelto.

Empezamos muy guarro, ya iremos mejorando.
Básicamente la clase puede implementar dos tipos de barco, segun el sistema de 
gobierno que empleen. En el primer tipo: 'waterjet', se supone que el sistema de 
propulsión gira con el sistema de gobierno: tipico de propolsion con waterjet o
motores fuera borda. En el segundo caso se supone que el sistema de propulsion
es fijo, apunta siempre hacia la popa y el barco es gobernado haciendo resisten
cia con un timon o similar. El los dos casos, el gobierno afecta a la velocidad
de avante del barquito, aunque de manera distinta. 
etc.
@author: juan
"""

import numpy as np
from numpy.linalg.linalg import norm 
import matplotlib.pyplot as pl
from matplotlib.path import Path
import matplotlib.patches as patches
import bezier_cvr
#import barchain_matrix_ini
class barco:
            
    def __init__(self, tipo = 'waterjet'):
        
    #Definiciones de variables de estado del barco
       self.pb = np.array([0.0,0.0]) #posicion x,y
       self.vb = np.array([0.0,0.0]) #velocidad vx, vy
       self.ab = np.array([0.0,0.0]) # aceleracion ax, ay
       self.alfab = 0.0 #aceleracion angular del barco
       self.wb = 0.0 #velocidad de giro del varvo en torno a su centro de masas
       self.theta = np.pi/2.0 #orientacion del barco respecto al eje x
       self.Fm = 0.0; #fuerza aplicada al barco
       self.setl = 0.0 #valor de consigna motor
       self.setw = 0.0 #valor de consigna timon
       self.interrv = 0.0 #integral del error de velocidad
       self.interrumbo = 0.0 #integral del error de rumbo
       self.thewj = 0.0 #orientacion del waterject o timon
       self.thewjmax = np.pi/6.0 # Angulo maximo del waterject o timon
       self.Ac = 10.0 #pdistancia del waterjet al centro de giro del barco???        
       self.M = 0.0 #momento de giro aplicado al barco por la fuerza del waterjet
       self.mb = 1000.0 #masa del barco en kg
       self.mul = 1000.0 #resistencia al avance (lineal con la velocidad)
       self.Ib = 1000.0 #momento de inercia del barco respecto a centro de masas
       self.mua = 1000.0 #resitencia al giro proprcional a v angular (gui~nada, yaw???
       self.ls = 2.0 #eslora del barco ??
       self.mut = 10000.0 #resitencia al avance en la direccion de deriva (sway)
       self.pmax = 50.0*735.0 #Fuerza maxima motor
       self.Tsf = 0.1 #constante de tiempo del motor (fuerza)
       self.Tsw = 0.1 #constante de tiempo waterject
       self.Ab = 0.1 #
       self.Aw = 0.1 #
       self.cfr = 10 #factor de fv = repeat([0.,5.],[1,puntos.shape[0]])orma del timon. (solo es aplicable a barcos
                     # con timon)
       self.link = 0 #este sería el eslabon de la cadena al que esta enganchado
       self.tipo = tipo
       #segun el tipo de barco se le asigna un sistema de propulsion y unos
       #valores a las constantes del controlador...
       if self.tipo == 'waterjet':
           self.propulsion = wj_prop
           self.kvpid = [100000.0, 10000.5, 10000.5] #Constantes del controlador
                                              #pid de velocidad. ojo el orden
                                              #de la lista es Kp, Kd, kI
           self.krpid = [8, 0.5, 0.05] #Constantes del controlador pid de rumbo
       else:
           self.propulsion = tim_prop
           self.kvpid = [10000.0, 1000.5, 1000.5] #Constantes del controlador
                                              #pid de velocidad. ojo el orden
                                              #de la lista es Kp, Kd, kI
           self.krpid = [5, 0.1, 0.01] #Constantes del controlador pid de rumbo
           
               
    def cstb(self,rt):
        '''calcula cambio de sistema de referencia de tierra a barco rt debe ser
        un vector de coordenadas en el sistema de referencia tierra la funcion
        devuelve las coordenadas del mismo punto medidas desde el sistema de 
        referencia barco'''
        tc = np.cos(self.theta)
        ts = np.sin(self.theta)
        rb = np.dot(np.array([[tc,ts],[-ts,tc]]),rt)
        return(rb)
        
    def csbt(self,rb): 
        '''calcula cambio de sistema de referencia de barco a tierra rb debe ser
        un vector de coordenadas en el sistema de referencia barco la funcion
        devuelve las coordenadas del mismo punto medidas desde el sistema de 
        referencia tierra'''
        tc = np.cos(self.theta)
        ts = np.sin(self.theta)
        rt = np.dot(np.array([[tc,-ts],[ts,tc]]),rb)
        return(rt)
        
    def dibujar(self):
        '''dibuja el barquito en su posición y orientación actual'''
        pl.plot(self.pb[0],self.pb[1],'+')
        vertices = np.array([[-1.,-0.25],[-1.,0.25],[-0.25,0.35],[1,0],\
        [-0.25,-0.35],[-1.,-0.25]])
        rot = np.array([[np.cos(self.theta),- np.sin(self.theta)],[np.sin(self.theta),\
        np.cos(self.theta)]])  
        vertrot = np.array([np.dot(rot,j) for j in vertices]) +\
        [self.pb[0],self.pb[1]]       
        codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
        Path.CURVE3]
        pathi = Path(vertrot,codes)
        patchi = patches.PathPatch(pathi,facecolor = 'blue')
        pl.gca().add_patch(patchi)
        
    def movimientotst(self,delta = 5e-4,extf = [0,0],dfcm = [0,0] ):
        '''Calcula el movimiento de un barco. Es igual que movimiento salvo 
        que se ha forzado un valor constante para la Fm. Solo
        sirve para probar si el modelo funciona, comparando con el de
        matlab esta funcion deberia ser borrada con el tiempo. La version
        actual es de 10 de junio de 2014
        y la presencia de una fuerza exterior, por defecto la fuerza exterior
        es cero y el punto de aplicacion el centro de masas tambien. En caso 
        de que se a~nada una fuerza exterior hay que a~nadir tambien el punto
        de aplicacion (IMPORTANTE: el punto de aplicacion debe darse en el
        sistema de coordenadas Barco). 
        Delta representa el paso de integracion(tiempo) empleado para calcular
        el movimiento del barco. '''
        
#        self.ab[0] = (self.Fm * np.cos(self.theta)\
#        - self.mul * np.cos(self.theta) \
#        * (self.vb[1] * np.sin(self.theta)\
#        + self.vb[0] * np.cos(self.theta))\
#        - self.mut * self.ls * np.sin(self.theta)\
#        * (self.vb[0] * np.sin(self.theta)\
#        - self.vb[1] * np.cos(self.theta)) + extf[0])/ self.mb
#        #eje y
#        self.ab[1] = (self.Fm * np.sin(self.theta)\
#        - self.mul * np.sin(self.theta) \
#        *(self.vb[1] * np.sin(self.theta)\
#        + self.vb[0] * np.cos(self.theta))\
#        - self.mut * self.ls * np.cos(self.theta)\
#        * (-self.vb[0] * np.sin(self.theta)\
#        +self.vb[1] * np.cos(self.theta)) + extf[1])  / self.mb
#        #velocidad del barco
#        self.vb[0] = self.ab[0] * delta + self.vb[0]
#        self.vb[1] = self.ab[1] * delta + self.vb[1] 
#        #posicion del barco.
#        self.pb[0] = self.vb[0] * delta + self.pb[0]
#        self.pb[1] = self.vb[1] * delta + self.pb[1]
#        #acelaracion angular
#        self.alfab = (self.M - self.mua * self.ls * self.wb - dfcm[1]\
#        * (extf[0] * np.cos(self.theta) +extf[1] * np.sin(self.theta)) +\
#        dfcm[0]\
#        * (extf[1] * np.cos(self.theta) - extf[0] * np.sin(self.theta)))/self.Ib
        self.propulsion(self,extf,dfcm)
        #velocidad del barco
        self.vb[0] = self.ab[0] * delta + self.vb[0]
        self.vb[1] = self.ab[1] * delta + self.vb[1] 
        #posicion del barco.
        self.pb[0] = self.vb[0] * delta + self.pb[0]
        self.pb[1] = self.vb[1] * delta + self.pb[1]
        #velocidad angular
        self.wb += self.alfab * delta 
        #rumbo
        self.theta += self.wb * delta
        
    def movimiento(self,delta = 5e-4,extf = [0,0],dfcm = [0,0] ):
        '''Calcula el movimiento de un barco. Bajo su propia fuerza impulsora y 
        la presencia de una fuerza exterior, por defecto la fuerza exterior
        es cero y el punto de aplicacion el centro de masas tambien. En caso 
        de que se a~nada una fuerza exterior hay que a~nadir tambien el punto
        de aplicacion (IMPORTANTE: el punto de aplicacion debe darse en el
        sistema de coordenadas Barco). 
        Delta representa el paso de integracion(tiempo) empleado para calcular
        el movimiento del barco'''
        
        #Vamos a suponer que la actuacción, tanto para incrementar la fuerza
        #coo para mover el waterjet o timon o lo que sea, responde como un
        #sistema de primer orden. A la vez, saturamos la actuacción a un valor
        #máximo para la fuerza y a un valor +thetawj y -thetawj para la orienta
        #cion del timón (referida al eje del barco)
        
        #saturamos la actuacion, (odio el if else y la cpu dando saltos ;)
        self.setl = ((self.setl >= 0) & (self.setl <= self.pmax)) * self.setl\
        + self.pmax * (self.setl >self.pmax )        
        
        #supongo que el timon es perfecto y simetrico....
        self.setw = \
        ((self.setw <= self.thewjmax) & (self.setw >= -self.thewjmax))\
        * self.setw\
        + self.thewjmax * (self.setw >self.thewjmax ) \
        - self.thewjmax * (self.setw < -self.thewjmax )
        
        #ajustamos los valores con un sistema de primer orden
        self.Fm = self.setl - self.Tsf*self.Fm *delta       
        self.thewj = self.setw - self.Tsw * self.thewj * delta
        
        self.propulsion(self,extf,dfcm)
                      
        #velocidad del barco
        self.vb[0] += self.ab[0] * delta
        self.vb[1] += self.ab[1] * delta 
        #posicion del barco.
        self.pb[0] += self.vb[0] * delta
        self.pb[1] += self.vb[1] * delta
        
        #velocidad angular
        self.wb += self.alfab * delta 
        #rumbo
        self.theta += self.wb * delta
        #pasamos el rumbo a angulos entre pi y -pi
        self.theta = self.theta -np.fix(self.theta/2./np.pi)\
        -2 * np.pi *((self.theta -np.fix(self.theta/2./np.pi)) > np.pi)\
        +2. * np.pi * ((self.theta -np.fix(self.theta/2./np.pi)) <= -np.pi)
        
    def controlador(self,c_rumbo=0, c_velocidad=0,delta=0):
        '''Esta funcion implementa un controlador básico que recibe consignas
        de rumbo y velocidad de avance para el barco, delta es el paso de inte
        gración que se está empleando en la simulación.'''
        
        #control de velocidad        
        errorv = c_velocidad - (np.cos(self.theta) * self.vb[0] +\
        np.sin(self.theta) * self.vb[1])
        
        self.interrv = self.interrv + errorv * delta
         
        
        self.setl = self.kvpid[0] * errorv  - \
        self.kvpid[1] * (np.cos(self.theta) * self.ab[0] \
        + np.sin(self.theta)) * self.ab[1] + \
        self.kvpid[2] + self.interrv
        
        #Control de rumbo
        
        errumbo = c_rumbo - self.theta \
        -2. * np.pi * ((c_rumbo - self.theta) > np.pi)\
        +2. * np.pi * ((c_rumbo - self.theta) <= -np.pi)
        
        self.interrumbo = self.interrumbo + errumbo * delta \
        * (np.abs(self.setw) < np.abs(self.thewjmax)) #antiwindup del integral        
        self.setw = -self.krpid[0] * errumbo + \
        self.krpid[1] * self.wb - \
        self.krpid[2] * self.interrumbo
        
        print 'c_rumbo', c_rumbo, 'rumbo', self.theta
        print 'errumbo = ', errumbo, ',  ', 'errorv = ', errorv, '\n'
        
    def planificador(self, puntos, rumbos, veloc = None, tiempos = None,
    ndat =10.  ):
        ''' Este planificador une con curvas de bezier de grado tres los puntos 
        de paso suministrados en el array puntos,de modo que el barco quede 
        orientado de acurdo con los angulos suminstrados en rumbos. si se le
        asignan velocidades a los puntos de paso, la curva de becier se ajustara
        para que el barco tome dicha velocidades, por ultimo si se asigna tiempos
        se calculara la trayectoria para que el barco pase por los puntos de 
        paso en los tiempos establecidos. De momento el planificador no tiene 
        en cuenta los límites impuestos por la propia dinamica del barco,  lo 
        cual hace que sea posible planificar trayectorias irrealizables... '''
        
        if veloc is None:
             #si no se especifican las velocidades por los puntos de paso,
             #supondremos que el barco parte del reposo y que las velocidades
             #se fijan a 5 m/s. Por supuesto estos valores por defecto
             #podrían cambiarse. Quizá lo más cómodo sería definirlos como
             #atributos del barco...
             veloc = np.repeat([0.,5.],[1,puntos.shape[0]-1])
        if tiempos is None:
            #si no se especifican los tiempos de paso, estos se obtienen a 
            #partir de las velocidades medias y la distancia entre los puntos.
            #Esto puede suponer que el barco tenga que acelerar o desacelerar 
            #mas de la cuenta en el tramo entre puntos de paso... pero pa 
            #probar nos vale
            
            tiempos = 2 * norm(puntos[1:] - puntos[:-1],axis =1).squeeze() / \
            (veloc[1:] + veloc[:-1])
            
            print tiempos
            
        lstplan = []
            
        for i,j,k,l,m,n,o in zip(puntos[:-1],puntos[1:],veloc[:-1],
        veloc[1:],tiempos, rumbos[:-1], rumbos[1:]):
                    #suponemos que el barco parte del repose en el primer punto
                    #elegimos una velocidad constate, 5m/s fija, en realidad
            print i, '\n'
            print j, '\n'
            print k, '\n'
            print l, '\n'
            print m, '\n'
            
            lstplan.append(bezier_cvr.bezier4p(i,j,k,l,m,n,o,ndat))
        
        return lstplan
            

    def resetear(self):
        self.pb = np.array([0.0,0.0]) #posicion x,y
        self.vb = np.array([0.0,0.0]) #velocidad vx, vy
        self.ab = np.array([0.0,0.0]) # aceleracion ax, ay
        self.alfab = 0.0 #adeleracion angular del barco
        self.wb = 0.0 #velocidad de giro del varvo en torno a su centro de masas
        self.theta = np.pi/2.0 #orientacion del barco respecto al eje x
        self.Fm = 0.0; #fuerza aplicada al barco
        self.setl = 0.0 
        self.setw = 0.0
        self.thewj = 0.0 #orientacion del waterject
        self.thewjmax = np.pi/6.0 # Angulo maximo del waterject
        self.Ac = 10.0 #pdistancia del waterjet al centro de giro del barco???        
        self.M = 0.0 #momento de giro aplicado al barco por la fuerza del waterjet
        self.mb = 1000.0 #masa del barco en kg
        self.mul = 1000.0 #resistencia al avance (lineal con la velocidad)
        self.Ib = 1000.0 #momento de inercia del barco respecto a centro de masas
        self.mua = 1000.0 #resitencia al giro proprcional a v angular (gui~nada, yaw???
        self.ls = 2.0 #eslora del barco ??
        self.mut = 10000.0 #resitencia al avance en la direccion de deriva (sway)
        self.pmax = 50.0*735.0 #potencia aplicada a waterjet)
        self.pmaxw = 50.0*735.0 #potencia maxima waterjet
        self.Ab = 0.1 #
        self.Aw = 0.1 #
        
    def listar(self):
        print "Caracteristicas:\n"
        print 'pb = posicion', self.pb 
        print 'vbc= velocidad', self.vb
        print 'ab = aceleracion', self.ab
        print 'alfab = aceleracion angular', self.alfab
        print 'wb = velocidad angular', self.wb 
        print 'theta = rumbo', self.theta
        print 'Fm = fuerza aplicada', self.waterjet
        print 'setl = consigna motor', self.setl 
        print 'setw = consigna rumbo', self.setw
        print 'thewj = orientacion waterjet', self.thewj
        print 'thewjmax = tope de giro waterjet', self.thewjmax
        print 'Ac = ', self.Ac        
        print 'M = momento del barco', self.M
        print 'mb = masa del barco', self.mb
        print 'mul coef. resistencia avance', self.mul
        print 'Ib = momento inercia', self.Ib 
        print 'mua = coef. resistencia al giro', self.mua 
        print 'ls = eslora', self.ls
        print 'mut = coef resistencia a desplazaniento lateral', self.mut
        print 'pmax = potencia maxima aplicable al motor', self.pmax
        print 'pmaxw = potencia maxima aplicable al giro',self.pmaxw
        print 'Ab', self.Ab
        print 'Aw', self.Aw
        
        

def trazar(barco,delta,rumbo,vel,tf,paso):
    ''' funcion para obtener la trasyectoria de un barco 
        delta: paso de integracion
        tf: tiempo final
        paso: intervalo para dibujar paso debe ser >= que delta 
        Tenemos que explicar mas cosas en esta función porque admite un
        controlador que podría no estar...'''
    
    t = np.arange(0.0,tf,delta)
    pl.subplot(2,2,1)
    pl.hold(True)
    pl.axis('equal')
    pl.xlabel('x')
    pl.ylabel('y')
    pl.subplot(2,2,2)
    pl.hold(True)
    pl.xlabel('t')
    pl.ylabel('rumbo')
    pl.subplot(2,2,3)
    pl.hold(True)
    pl.xlabel('t')
    pl.ylabel('velocidad')
    pl.subplot(2,2,4)
    pl.hold(True)
    pl.xlabel('t')
    pl.ylabel('setl/1000(s)y setw/0.5(t)')
    tp = t[0]
    thp = barco.theta
    vp = barco.vb
    bslp = barco.setl
    bsvp = barco.setw   
    for i in range(np.size(t)):
        barco.controlador(rumbo,vel,delta)
        barco.movimiento(delta = delta)
        
        if i%paso == 0:
            pl.subplot(2,2,1)
            barco.dibujar()
#            
#            pl.plot(barco.pb[0],barco.pb[1],'+')
#            vertices = np.array([[-1.,-0.25],[-1.,0.25],[-0.25,0.35],[1,0],\
#            [-0.25,-0.35],[-1.,-0.25]])
#            rot = np.array([[np.cos(barco.theta),- np.sin(barco.theta)],[np.sin(barco.theta),\
#            np.cos(barco.theta)]])  
#            vertrot = np.array([np.dot(rot,j) for j in vertices]) +\
#            [barco.pb[0],barco.pb[1]]       
#            codes = [Path.MOVETO,Path.LINETO,Path.CURVE3,Path.CURVE3,Path.CURVE3,\
#            Path.CURVE3]
#            pathi = Path(vertrot,codes)
#            patchi = patches.PathPatch(pathi,facecolor = 'blue')
#            pl.gca().add_patch(patchi)
           
            pl.subplot(2,2,2)
            pl.plot([tp,t[i]],[thp,barco.theta],'b')
            
            pl.subplot(2,2,3)
            pl.plot([tp,t[i]],[np.cos(thp) * vp[0] +
            np.sin(thp) * vp[1],np.cos(barco.theta) * barco.vb[0] +
            np.sin(barco.theta) * barco.vb[1]],'r')
            
            pl.subplot(2,2,4)
            pl.plot([tp,t[i]],[bslp/10000,barco.setl/10000],'k')
            pl.plot([tp,t[i]],[bsvp/0.5,barco.setw/0.5],'g')
            #pl.plot(t[i],barco.Fm,'o')
            
            #pl.pause(0.01)
            
            tp = t[i]    
            thp = barco.theta
            vp = barco.vb.copy()
            bslp = barco.setl
            bsvp = barco.setw 
        

def wj_prop(barco, extf, dfcm):
    #funcion empleada si se emplea un waterjet, motor fuera borda o similar
    #para obtener los valores de las aceleraciones lineal y angular,
    barco.ab[0] = (barco.Fm * np.cos(barco.theta + barco.thewj)\
    - barco.mul * np.cos(barco.theta) \
    * (barco.vb[1] * np.sin(barco.theta)\
    + barco.vb[0] * np.cos(barco.theta))\
    - barco.mut * barco.ls * np.sin(barco.theta)\
    * (barco.vb[0] * np.sin(barco.theta)\
    - barco.vb[1] * np.cos(barco.theta)) + extf[0])/ barco.mb
    
    #eje y 
    barco.ab[1] = (barco.Fm * np.sin(barco.theta+barco.thewj)\
    - barco.mul * np.sin(barco.theta) \
    *(barco.vb[1] * np.sin(barco.theta)\
    + barco.vb[0] * np.cos(barco.theta))\
    - barco.mut * barco.ls * np.cos(barco.theta)\
    * (-barco.vb[0] * np.sin(barco.theta)\
    +barco.vb[1] * np.cos(barco.theta)) + extf[1])  / barco.mb
    
    #Momento generado por el waterjet sobrel eje del barco
    barco.M = - barco.Ac * barco.Fm * np.sin(barco.thewj)
    
    #aceleración angular        
    barco.alfab = (barco.M - barco.mua * barco.ls * barco.wb - dfcm[1]\
    * (extf[0] * np.cos(barco.theta) +extf[1] * np.sin(barco.theta)) +\
    dfcm[0]\
    * (extf[1] * np.cos(barco.theta)\
    - extf[0] * np.sin(barco.theta)))/barco.Ib
    
def tim_prop(barco,extf,dfcm):
    #funcion empleada si se emplea un timon
    #para obtener los valores de las aceleraciones lineal y angular,
    barco.ab[0] = (barco.Fm * np.cos(barco.theta)\
    - barco.mul * np.cos(barco.theta) \
    * (barco.vb[1] * np.sin(barco.theta)\
    + barco.vb[0] * np.cos(barco.theta))\
    - barco.mut * barco.ls * np.sin(barco.theta)\
    * (barco.vb[0] * np.sin(barco.theta)\
    - barco.vb[1] * np.cos(barco.theta)) + extf[0]
    + barco.cfr*(-barco.vb[0]*np.sin(barco.theta +barco.thewj)
    + barco.vb[1]*np.cos(barco.theta +barco.thewj))\
    * np.sin(barco.theta +barco.thewj)) / barco.mb
    
    #eje y 
    barco.ab[1] = (barco.Fm * np.sin(barco.theta)\
    - barco.mul * np.sin(barco.theta) \
    * (barco.vb[1] * np.sin(barco.theta)\
    + barco.vb[0] * np.cos(barco.theta))\
    - barco.mut * barco.ls * np.cos(barco.theta)\
    * (-barco.vb[0] * np.sin(barco.theta)\
    + barco.vb[1] * np.cos(barco.theta)) + extf[1]
    - barco.cfr*(-barco.vb[0]*np.sin(barco.theta +barco.thewj)
    + barco.vb[1]*np.cos(barco.theta +barco.thewj))\
    * np.cos(barco.theta +barco.thewj)) / barco.mb
    
    #Momento generado por el waterjet sobrel eje del barco
    barco.M =  barco.Ac * barco.cfr * \
    (-barco.vb[0]*np.sin(barco.theta +barco.thewj)
    + barco.vb[1]*np.cos(barco.theta +barco.thewj))
    
    #aceleración angular        
    barco.alfab = (barco.M - barco.mua * barco.ls * barco.wb - dfcm[1]\
    * (extf[0] * np.cos(barco.theta) +extf[1] * np.sin(barco.theta)) +\
    dfcm[0]\
    * (extf[1] * np.cos(barco.theta)\
    - extf[0] * np.sin(barco.theta)))/barco.Ib
    