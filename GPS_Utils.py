# -*- coding: utf-8 -*-
"""
Created on Wed May 10 14:27:48 2017
Ete fichero de momento copia los ficheros en c que compila el barquito para
hacer calculos de posición a partir de coordenas GPS. Segun nos aclaremos lo
iremos haciendo mehon...
@author: juan

_______________________________________________________________________________
DATOS MISTERIOSOS CONSIGNADOS AL PRINCIPIO DE utils.h

EARTH_RADIUS 6371000 #middle Earth radius 6370 km
MADRID_ALTITUDE 655  #Madrid ALTITUDE 655m

m1 = 111132.92       #latitude calculation term 1
m2 = -559.82         #latitude calculation term 2
m3 = 1.175           #latitude calculation term 3
m4 = -0.0023         #latitude calculation term 4

p1 = 111412.84       #longitude calculation term 1
p2 = -93.5           #longitude calculation term 2
p3 = 0.118           #longitude calculation term 3
_______________________________________________________________________________
"""

#locations are implemented as a dictionary 
#each list contains [centre_x,centre_y] [longitude,latitude]
locations = dict()
#Pista atletismo complu
locations['PISTA'] =  [-3.705378,40.375407]
#estanque tierno galvan
locations['TIERNO'] = [-3.682995,40.386326]
#estanque de pradolongo
locations['PRADOL'] = [-3.705378,40.375407]
#embalse del atazar
locations['ATAZAR'] = [-3.535306,40.917227]
#More can be added... when needed

working_area_rad = 1500 #to be discover what the hell is this

#define ahora unos fatores para la latitud y la longitud que se obtendrían
#a partir de la localizacion. Le doy forma funciional aunque creo que es la 
#tipica cosa que se usa en la inicializaci-on y nunca más...

def latlon_factor(lat):
    lat_factor = 111132.92 - 559.82 * cos(2*lat)
    lon_factor = 