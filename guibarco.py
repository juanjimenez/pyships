# -*- coding: utf-8 -*-
"""
Created on Sat May 24 10:57:38 2014
Eso  es un intento de una ventana para crear y ajustar los parametros de un 
barco... de los que creabarcosolo.py
Es un primer intento...
@author: juan
"""

import guidata
_app = guidata.qapplication()

import guidata.dataset.datatypes as dt
import guidata.dataset.dataitems as di

class Processing(dt.DataSet):
    """
    Titulillo
    y lo que escribo aquí lo comenta en algún sitio
    """
    a = di.FloatItem('posicion x', default = 0.0)
    b = di.FloatItem('position y', default = 0.0)
    floatarray = di.FloatArrayItem("Float array", default=ones( (50, 5), float),
                                format=" %.2e ").set_pos(col=1)
param = Processing()

param.edit()