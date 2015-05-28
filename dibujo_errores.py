# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 11:19:57 2015
Usada para dibujar errores... Solo se uso ese día para ver que errores salian 

@author: juan
"""
figure()
a = gca()
a.axes.set_yscale('log')
hold(True)
plot(arange(6), abs(maxerr01[0]), 'o')
plot(arange(6), abs(maxerr001[0]), '^r')
plot(arange(6), abs(maxerr0001[0]), 'vg')
plot(arange(6), abs(maxerr00001[0]),'sk')
grid(which = 'both')
xlim(-0.1,7)
legend(('$\Delta t = 0.01s,\\ n = 10^4$','$\Delta t = 0.001s,\\ n = 10^5$'\
,'$\Delta t = 10^{-4}s,\\ n = 10^6$','$\Delta t = 10^{-5}s,\\ n = 10^7$'),\
numpoints = 1,fontsize = 11)

figure()
a = gca()
a.axes.set_yscale('log')
hold(True)
plot(arange(6), abs(maxerr01[1]), 'o')
plot(arange(6), abs(maxerr001[1]), '^r')
plot(arange(6), abs(maxerr0001[1]), 'vg')
plot(arange(6), abs(maxerr00001[1]),'sk')
grid(which = 'both')
xlim(-0.1,7)
legend(('$\Delta t = 0.01s,\\ n = 10^4$','$\Delta t = 0.001s,\\ n = 10^5$'\
,'$\Delta t = 10^{-4}s,\\ n = 10^6$','$\Delta t = 10^{-5}s,\\ n = 10^7$'),\
numpoints = 1, fontsize = 11)

figure()
a = gca()
a.axes.set_yscale('log')
hold(True)
plot(arange(6), sqrt(maxerr01[1]**2 + maxerr01[0]**2), 'o')
plot(arange(6), sqrt(maxerr001[1]**2 + maxerr001[0]**2), '^r')
plot(arange(6), sqrt(maxerr0001[1]**2 + maxerr0001[0]**2), 'vg')
plot(arange(6), sqrt(maxerr00001[1]**2 + maxerr00001[0]**2),'sk')
grid(which = 'both')
xlim(-0.1,7)
legend(('$\Delta t = 0.01s,\\ n = 10^4$','$\Delta t = 0.001s,\\ n = 10^5$'\
,'$\Delta t = 10^{-4}s,\\ n = 10^6$','$\Delta t = 10^{-5}s,\\ n = 10^7$'),\
numpoints = 1, fontsize = 11)