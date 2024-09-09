import math as m
import sys
#import matplotlib.pyplot as plt
import numpy as np
# Save constants
a = (3.0*3.7)**0.5
b = 6
l = 2
f = (1.11*0.64)**0.5
Rmin = 0.01
rcut = 2.25
emix = 1.2
# Find X,Y and Z from rcut constraint

N = 500
r = np.linspace(Rmin,rcut,num=N)
V = []
F = []
for i in r:
     V.append(emix*m.exp(a*(1-i/f) - b*(1-i/f)**l*m.log(i/f)))
     F.append(-emix*m.exp(a*(1-i/f) - b*(1-i/f)**l*m.log(i/f))*(-a/f - b*(1-i/f)**l/i + b/f*l*(1-i/f)**(l-1)*m.log(i/f)))


#plt.plot(r,V)
#plt.show()
#plt.plot(r,F)
#plt.show()

NAME = "YK" + "Mix" + ".pot"

G = open(NAME,'w')
G.write('# DATE: 2022-03-08 UNTIS: lj CONTRIBUTOR: Luis Nieves\n')
G.write('# YK potential due to YK 1976, a = '+str(a)+', f = ' + str(f) +' and eps_mix = ' +str(emix) +' for mixed interaction\n')
G.write('\n')
G.write('YKMix\n')
G.write('N '+ str(N) +'\n')
G.write('\n')
for i in range(0,N):
	G.write(str(i+1) + ' '+ str(r[i]) + ' ' + str(V[i]) + ' ' +str(F[i]) +'\n')
G.close()
