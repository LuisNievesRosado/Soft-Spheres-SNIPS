import pyscal.core as pc
import pyscal.crystal_structures as pcs
import os
from pyscal.trajectory import Trajectory
import matplotlib.pyplot as plt
import numpy as np
import sys as s
import statistics



BL1 = 1.0
BL2 = 2.0
TL1 = 39.15
TL2 = 40.15



traj = Trajectory('SS2.lammpstrj')
print('Trajectory read complete')

if not os.path.exists('LAYER1'):
	os.makedirs('LAYER1')
	
if not os.path.exists('LAYER2'):
	os.makedirs('LAYER2')

if not os.path.exists('VBOT1'):
	os.makedirs('VBOT1')
	
if not os.path.exists('VTOP1'):
	os.makedirs('VTOP1')

if not os.path.exists('VBOT2'):
	os.makedirs('VBOT2')

if not os.path.exists('VTOP2'):
	os.makedirs('VTOP2')	
		
if not os.path.exists('HISTO'):
	os.makedirs('HISTO')		

if not os.path.exists('TIME'):
	os.makedirs('TIME')		


TT = traj.nblocks
N = traj.natoms

if ANTYPE == 'ini':
	ANRANGE = range(TT-2,TT-1)
elif ANTYPE == 'all':
	ANRANGE = range(0,TT-1,5)


count = 0

L1_ISVO = []
L1_SVPS = []

L2_ISVO = []
L2_SVPS = []

L1P_ISVO = []
L1P_SVPS = []

L2P_ISVO = []
L2P_SVPS = []

for j in ANRANGE:
	print('Processing timestep ' + str(j+1) + ' of ' + str(str(TT)))
	count = count + 1
	SL = traj[j].to_dict()
	
	
	BOXX1 = SL[0]["box"][0][0]
	BOXX2 = SL[0]["box"][0][1]
	
	BOXY1 = SL[0]["box"][1][0]
	BOXY2 = SL[0]["box"][1][1]
	
	BOXZ1 = SL[0]["box"][2][0]
	BOXZ2 = SL[0]["box"][2][1]
	
	
	ID = SL[0]["atoms"]["id"]
	XS = SL[0]["atoms"]["xs"]
	YS = SL[0]["atoms"]["ys"]
	ZS = SL[0]["atoms"]["zs"]
	TY = SL[0]["atoms"]["type"]
	
	X = []
	Y = []
	Z = []
	XO = []
	YO = []
	ZO = []
	
	for i in range(0,N):
#		X.append([ID[i],XS[i]*(BOXX2-BOXX1)+BOXX1])
#		Y.append([ID[i],YS[i]*(BOXY2-BOXY1)+BOXY1])
#		Z.append([ID[i],ZS[i]*(BOXZ2-BOXZ1)+BOXZ1])
		XO.append([ID[i],TY[i],XS[i]*(BOXX2-BOXX1)])
		YO.append([ID[i],TY[i],YS[i]*(BOXY2-BOXY1)])
		ZO.append([ID[i],TY[i],ZS[i]*(BOXZ2-BOXZ1)])    
	
#	X.sort()
#	Y.sort()
#	Z.sort()
	
	XO.sort()
	YO.sort()
	ZO.sort()
	
	
	
	
	ZH1 = []
	ZH2 = []

	
	for i in range(0,N):
		if int(ZO[i][1]) == 1:
			ZH1.append(ZO[i][2])
		elif int(ZO[i][1]) == 2:
			ZH2.append(ZO[i][2])


	
	BOT1 = []
	BOT2 = []
	TOP1 = []
	TOP2 = []
	
	for i in range(0,N):
		if ZO[i][2] < BL1:
			BOT1.append([int(XO[i][0]),int(XO[i][1]),XO[i][2],YO[i][2],ZO[i][2]])
		elif ZO[i][2] > BL1 and ZO[i][2] < BL2:
			BOT2.append([int(XO[i][0]),int(XO[i][1]),XO[i][2],YO[i][2],ZO[i][2]])
		elif ZO[i][2] > TL1 and ZO[i][2] < TL2:
			TOP2.append([int(XO[i][0]),int(XO[i][1]),XO[i][2],YO[i][2],ZO[i][2]])
		elif ZO[i][2] > TL2:
			TOP1.append([int(XO[i][0]),int(XO[i][1]),XO[i][2],YO[i][2],ZO[i][2]])
		
	plt.figure()
	plt.hist(ZH1,bins = 100,color='m',label='ISV')
	plt.hist(ZH2,bins = 100,color='r',alpha = 0.5,label='SVPS')
	plt.axvspan(0,2.0,alpha=0.5,color='k')
	plt.axvspan(BOXZ2-BOXZ1-2.0,BOXZ2-BOXZ1,alpha=0.5,color='k')
	plt.axvspan(0,0.5,alpha=0.5,color='k')
	plt.axvspan(BOXZ2-BOXZ1-0.5,BOXZ2-BOXZ1,alpha=0.5,color='k')
	plt.axvline(x=BL1,color='k',linestyle='--')
	plt.axvline(x=BL2,color='k',linestyle='--')
	plt.axvline(x=TL1,color='k',linestyle='--')
	plt.axvline(x=TL2,color='k',linestyle='--')
	plt.legend()
	plt.ylabel('N')
	plt.xlabel('Z Position')
	plt.savefig('HISTO\Histo'+str(count)+'.png',bbox_inches='tight')
	if ANTYPE == 'ini':
		plt.show()
	plt.close()
	
		
	#Compute distribution
	
	countO = 0
	countS = 0
	
	for i in BOT1:
		if i[1] == 1:
			countO = countO + 1
		elif i[1] == 2:
			countS = countS + 1	
	for i in TOP1:
		if i[1] == 1:
			countO = countO + 1
		elif i[1] == 2:
			countS = countS + 1	

	L1_ISVO.append(countO)	
	L1_SVPS.append(countS)

	L1P_ISVO.append(100*countO/(countS+countO))
	L1P_SVPS.append(100*countS/(countS+countO))	


	countI = 0
	countO = 0
	
	for i in BOT2:
		if i[1] == 1:
			countO = countO + 1
		elif i[1] == 2:
			countS = countS + 1	
	for i in TOP2:
		if i[1] == 1:
			countO = countO + 1
		elif i[1] == 2:
			countS = countS + 1	

	L2_ISVO.append(countO)	
	L2_SVPS.append(countS)

	L2P_ISVO.append(100*countO/(countS+countO))
	L2P_SVPS.append(100*countS/(countS+countO))	

		
			
	
	
	
	
	### Output for voronoi code
	G1 = open('VBOT1/BotL1_'+str(count)+'.inp','w')
	
	G1.write(str(BOXX2-BOXX1)+' 0.0 0.0\n')
	G1.write('0.0 '+ str(BOXY2-BOXY1)+ ' 0.0\n')
	G1.write('0.0 0.0 0.0\n')
	for i in range(0,len(BOT1)):
		G1.write(str(BOT1[i][2]) + ' ' + str(BOT1[i][3]) + ' 0.0\n')
	G1.close()	
	
	G2 = open('VBOT2/BotL2_'+str(count)+'.inp','w')
	
	G2.write(str(BOXX2-BOXX1)+' 0.0 0.0\n')
	G2.write('0.0 '+ str(BOXY2-BOXY1)+ ' 0.0\n')
	G2.write('0.0 0.0 0.0\n')
	for i in range(0,len(BOT2)):
		G2.write(str(BOT2[i][2]) + ' ' + str(BOT2[i][3]) + ' 0.0\n')
	G2.close()
	
	G3 = open('VTOP1/TopL1_'+str(count)+'.inp','w')
	
	G3.write(str(BOXX2-BOXX1)+' 0.0 0.0\n')
	G3.write('0.0 '+ str(BOXY2-BOXY1)+ ' 0.0\n')
	G3.write('0.0 0.0 0.0\n')
	for i in range(0,len(TOP1)):
		G3.write(str(TOP1[i][2]) + ' ' + str(TOP1[i][3]) + ' 0.0\n')
	G3.close()
	
	G4 = open('VTOP2/TopL2_'+str(count)+'.inp','w')
	
	G4.write(str(BOXX2-BOXX1)+' 0.0 0.0\n')
	G4.write('0.0 '+ str(BOXY2-BOXY1)+ ' 0.0\n')
	G4.write('0.0 0.0 0.0\n')
	for i in range(0,len(TOP2)):
		G4.write(str(TOP2[i][2]) + ' ' + str(TOP2[i][3]) + ' 0.0\n')
	G4.close()
	
	
###  Output as data file ###

	
	N1 = len(BOT1)+len(TOP1)
	G1 = open('LAYER1\Layer1_'+str(count)+'.data','w')
	
	G1.write('LAMMPS Description \n')
	G1.write('\n')
	G1.write(str(N1) + ' atoms \n')
	G1.write('0 bonds \n')
	G1.write('0 angles \n')
	G1.write('\n')
	G1.write('3 atom types \n')
	G1.write('0 bond types \n')
	G1.write('0 angle types \n')
	G1.write('\n')
	G1.write(str(BOXX1) + ' ' + str(BOXX2) + ' xlo xhi \n')
	G1.write(str(BOXY1) + ' ' + str(BOXY2) + ' ylo yhi \n')
	G1.write(str(BOXZ1) + ' ' + str(BOXZ2) + ' zlo zhi \n')
	G1.write('\n')
	G1.write('Masses \n')
	G1.write('\n')
	G1.write('1 1.0 \n')
	G1.write('2 1.0 \n')
	G1.write('3 1.0 \n')
	G1.write('\n')
	G1.write('Atoms \n')
	G1.write('\n')
	for i in range(0,len(BOT1)):
		G1.write(str(i + 1) + ' 1 ' + str(BOT1[i][1]) + ' 0 ' + str(BOT1[i][2])+ ' ' + str(BOT1[i][3])+ ' ' + str(BOT1[i][4]) + '\n')
	for i in range(0,len(TOP1)):
		G1.write(str(i + 1) + ' 1 ' + str(TOP1[i][1]) + ' 0 ' + str(TOP1[i][2])+ ' ' + str(TOP1[i][3])+ ' ' + str(TOP1[i][4]) + '\n')		
	G1.close()
	
	N2 = len(BOT2)+len(TOP2)
	G2 = open('LAYER2/Layer2_'+str(count)+'.data','w')
	G2.write('LAMMPS Description \n')
	G2.write('\n')
	G2.write(str(N2) + ' atoms \n')
	G2.write('0 bonds \n')
	G2.write('0 angles \n')
	G2.write('\n')
	G2.write('3 atom types \n')
	G2.write('0 bond types \n')
	G2.write('0 angle types \n')
	G2.write('\n')
	G2.write(str(BOXX1) + ' ' + str(BOXX2) + ' xlo xhi \n')
	G2.write(str(BOXY1) + ' ' + str(BOXY2) + ' ylo yhi \n')
	G2.write(str(BOXZ1) + ' ' + str(BOXZ2) + ' zlo zhi \n')
	G2.write('\n')
	G2.write('Masses \n')
	G2.write('\n')
	G2.write('1 1.0 \n')
	G2.write('2 1.0 \n')
	G2.write('3 1.0 \n')
	G2.write('\n')
	G2.write('Atoms \n')
	G2.write('\n')
	for i in range(0,len(BOT2)):
		G2.write(str(i + 1) + ' 1 ' + str(BOT2[i][1]) + ' 0 ' + str(BOT2[i][2])+ ' ' + str(BOT2[i][3])+ ' ' + str(BOT2[i][4]) + '\n')
	for i in range(0,len(TOP2)):
		G2.write(str(i + 1) + ' 1 ' + str(TOP2[i][1]) + ' 0 ' + str(TOP2[i][2])+ ' ' + str(TOP2[i][3])+ ' ' + str(TOP2[i][4]) + '\n')		
	G2.close()


plt.figure()
plt.plot(L1_ISVO,color='k',linewidth=2.0,label='ISVO')
plt.plot(L1_SVPS,color='r',linewidth=2.0,label='SVPS')
plt.legend()
plt.title('Distribution of Particles in Layer 1')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L1.png',bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(L2_ISVO,color='k',linewidth=2.0,label='ISVO')
plt.plot(L2_SVPS,color='r',linewidth=2.0,label='SVPS')
plt.legend()
plt.title('Distribution of Particles in Layer 2')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L2.png',bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(L1P_ISVO,color='k',linewidth=2.0,label='ISVO')
plt.plot(L1P_SVPS,color='r',linewidth=2.0,label='SVPS')
plt.legend()
plt.title('Percentage Distribution of Particles in Layer 1')
plt.ylabel('Percentage')
plt.xlabel('Time')
plt.savefig('TIME/L1P.png',bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(L2P_ISVO,color='k',linewidth=2.0,label='ISVO')
plt.plot(L2P_SVPS,color='r',linewidth=2.0,label='SVPS')
plt.legend()
plt.title('Percentage Distribution of Particles in Layer 2')
plt.ylabel('Percentage')
plt.xlabel('Time')
plt.savefig('TIME/L2P.png',bbox_inches='tight')
plt.close()

G1 = open('TIME/L1_ISVO.dat','w')
for i in L1_ISVO:
	G1.write(str(i) + '\n')
G1.close()

G1 = open('TIME/L2_ISVO.dat','w')
for i in L2_ISVO:
	G1.write(str(i) + '\n')
G1.close()

G1 = open('TIME/L1P_ISVO.dat','w')
for i in L1P_ISVO:
	G1.write(str(i) + '\n')
G1.close()

G1 = open('TIME/L2P_ISVO.dat','w')
for i in L2P_ISVO:
	G1.write(str(i) + '\n')
G1.close()

G1 = open('TIME/L1_SVPS.dat','w')
for i in L1_SVPS:
	G1.write(str(i) + '\n')
G1.close()

G1 = open('TIME/L2_SVPS.dat','w')
for i in L2_SVPS:
	G1.write(str(i) + '\n')
G1.close()

G1 = open('TIME/L1P_SVPS.dat','w')
for i in L1P_SVPS:
	G1.write(str(i) + '\n')
G1.close()

G1 = open('TIME/L2P_SVPS.dat','w')
for i in L2P_SVPS:
	G1.write(str(i) + '\n')
G1.close()