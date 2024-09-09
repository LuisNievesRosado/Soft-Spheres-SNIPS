import math as m
#import matplotlib.pyplot as plt
import numpy as np
# Save constants
a = 3.7
b = 6
l = 2
f = 0.64
Rmin = 0.01
rcut = 2.25

# Find X,Y and Z from rcut constraint

N = 500
r = np.linspace(Rmin,rcut,num=N)
V = []
F = []
for i in r:
     V.append(m.exp(a*(1-i/f) - b*(1-i/f)**l*m.log(i/f)))
     F.append(-m.exp(a*(1-i/f) - b*(1-i/f)**l*m.log(i/f))*(-a/f - b*(1-i/f)**l/i + b/f*l*(1-i/f)**(l-1)*m.log(i/f)))
#i = 2.25
#TEST = A*i**(-n) + l1*(1-m.tanh(k1*(i-d1))) + l2*(1-m.tanh(k2*(i-d2))) + X*i**2 + Y*i + Z	


#plt.plot(r,V)
#plt.show()
#plt.plot(r,F)
#plt.show()

NAME = "YK" + "BCC" + ".pot"

G = open(NAME,'w')
G.write('# DATE: 2022-03-08 UNTIS: lj CONTRIBUTOR: Luis Nieves\n')
G.write('# YK potential due to YK 1976, a = ' + str(a) + 'and f = ' +str(f) +  ' converges to BCC\n')
G.write('\n')
G.write('YKBCC\n')
G.write('N '+ str(N) +'\n')
G.write('\n')
for i in range(0,N):
	G.write(str(i+1) + ' '+ str(r[i]) + ' ' + str(V[i]) + ' ' +str(F[i]) +'\n')
G.close()
