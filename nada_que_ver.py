# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 17:41:34 2015
TONTERIAS DE DIBUJITOS NO tiene que ver con los barcos
@author: juan
"""
import matplotlib.animation as animation
i = [[0,0],[0,0]]


#primer paso
index = int(round(rand()))

i[index] = [i[index][1],2 * round(rand()) - 1]
fig1 = figure()
plot(i[0],i[1])

hold(True)

for j in range(15000):
    if index == 0:
        i[index][0] = i[index][1] 
        index = 1
        i[index] = [i[index][1],i[index][1] + 2 * round(rand()) - 1]
    else:
        i[index][0] = i[index][1]
        index = 0
        i[index] = [i[index][1],i[index][1] + 2 * round(rand()) - 1]
        
    #print i

    plot(i[0],i[1])
    pause(10e-6)
     
   