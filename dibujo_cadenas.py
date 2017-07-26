# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 11:25:02 2014
Esto solo se uso una vez para hacer el dibujo que aparece en el report de 
modelo de arrastre. Supone que antes se ha corrido el programa iterador.py
con las siguientes adaptaciones (los numeros al principio representan lineas
de codigo en el programa iterador.py: 

14 delta = 0.01
.
20 cd5 = boom.cadena(6)
.
52 tam = 165 #aqui definimos de momento el tamaño del experimento...
53 step =165
.
74 bi.controlador(pi/3,5,delta)
.   
76 bd.controlador(pi/3,5,delta) 
.
108 if (i+1)%step == 0:

Se trata de usar una cadena con poquitos eslabones para que la cosa quede
curiosa... el step de dibujo se hace igual al tamaño y el la condicion de
dibujo se toma en i+1 para que solo dibuje una posicion del experimento (la 
ultima).
EL codigo de este programa pinta los vectores posiciones de un par de eslabones
asi sus vectores 'paralelo' y 'perpendicular' usasdos en el modelo y los rotula
NI que decir tiene que si se cambian las condiciones del experimento los 
rotulos saldran descolocados
@author: juan
"""


arrow(cd5.cms[0,2],cd5.cms[1,2],0.5*cd5.normal[0,2],0.5*cd5.normal[1,2],\
length_includes_head=True,color='k')

arrow(cd5.cms[0,2],cd5.cms[1,2],0.5*cd5.para[0,2],0.5*cd5.para[1,2],\
length_includes_head=True,color='k')
arrow(0,0,cd5.cms[0,1],cd5.cms[1,1],length_includes_head=True,color='k')
arrow(0,0,cd5.cms[0,2],cd5.cms[1,2],length_includes_head=True,color='k')
arrow(cd5.cms[0,1],cd5.cms[1,1],0.5*cd5.normal[0,1],0.5*cd5.normal[1,1],\
length_includes_head=True,color='k')
arrow(cd5.cms[0,1],cd5.cms[1,1],0.5*cd5.para[0,1],0.5*cd5.para[1,1],\
length_includes_head=True,color='k')
#axis([0,6,0,7])
text(1,0.32,r'$\vec{r}_{i+1}$')
text(0.8,0.72,r'$\vec{r}_{i}$')


text(1.6,1.16,r'$\vec{n}_{i}$',rotation=180/pi*arctan(cd5.normal[1,1]/cd5.normal[0,1]))
text(1.23,1.25,r'$\vec{p}_{i}$',rotation=180/pi*arctan(cd5.para[1,1]/cd5.para[0,1])+90)
text(1.85,0.75,r'$\vec{n}_{i+1}$',rotation=180/pi*arctan(cd5.normal[1,2]/cd5.normal[0,2]))
text(1.6,0.58,r'$\vec{p}_{i+1}$',rotation=180/pi*arctan(cd5.para[1,2]/cd5.para[0,2]))
xlabel('x')
ylabel('y')
axis([0,3.5,0,3.5])
gca().axes.get_xaxis().set_ticklabels([])
gca().axes.get_yaxis().set_ticklabels([])