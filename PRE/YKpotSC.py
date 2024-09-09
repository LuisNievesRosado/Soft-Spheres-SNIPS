import math as m
import matplotlib.pyplot as plt
import numpy as np
# Save constants
a = 3.0
b = 6
l = 2
Rmin = 0.01
rcut = 2.25
f = 1.11

# Find X,Y and Z from rcut constraint

N = 500
r = np.linspace(Rmin,rcut,num=N)
V = []
F = []
V2 = []
for i in r:
     V.append(m.exp(a*(1-i/f) - b*(1-i/f)**l*m.log(i/f)))
     F.append(-m.exp(a*(1-i/f) - b*(1-i/f)**l*m.log(i/f))*(-a/f - b*(1-i/f)**l/i + b/f*l*(1-i/f)**(l-1)*m.log(i/f)))
#i = 2.25
#TEST = A*i**(-n) + l1*(1-m.tanh(k1*(i-d1))) + l2*(1-m.tanh(k2*(i-d2))) + X*i**2 + Y*i + Z	


#a = 5

#for i in r:
#     V2.append(m.exp(a*(1-i) - b*(1-i)**l*m.log(i)))
#ep = 7.5

#Ve = ep*(r-1)**2

#plt.plot(r,Ve,label='harmonic')
#plt.plot(r,V,label='SC')
#plt.plot(r,V2,label='hex')
#plt.legend()
#axes = plt.gca()

#axes.set_ylim([0,10])
#plt.show()
#plt.plot(r,F)
#plt.show()

NAME = "YK" + "SC" + ".pot"

G = open(NAME,'w')
G.write('# DATE: 2022-03-08 UNTIS: lj CONTRIBUTOR: Luis Nieves\n')
G.write('# YK potential due to YK 1976, a = ' + str(a) + 'and f = ' +str(f) +  ' converges to SC.\n')
G.write('\n')
G.write('YKSC\n')
G.write('N '+ str(N) +'\n')
G.write('\n')
for i in range(0,N):
	G.write(str(i+1) + ' '+ str(r[i]) + ' ' + str(V[i]) + ' ' +str(F[i]) +'\n')
G.close()
