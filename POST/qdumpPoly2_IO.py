import pyscal.core as pc
import pyscal.crystal_structures as pcs
import os
from pyscal.trajectory import Trajectory
import matplotlib.pyplot as plt
import numpy as np
import sys as s
import statistics
#N = int(s.argv[1])
#traj = Trajectory(INFILE)

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

count = 0

L1_ISVO = []
L1_ISVO1 = []
L1_ISVO2 = []
L1_ISVO3 = []
L1_ISVO4 = []
L1_ISVO5 = []
L1_ISVO6 = []
L1_ISVO7 = []
L1_ISVO8 = []
L1_ISVO9 = []
L1_ISV = []
L1_ISV1 = []
L1_ISV2 = []
L1_ISV3 = []
L1_ISV4 = []
L1_ISV5 = []
L1_ISV6 = []
L1_ISV7 = []
L1_ISV8 = []
L1_ISV9 = []


L2_ISVO = []
L2_ISVO1 = []
L2_ISVO2 = []
L2_ISVO3 = []
L2_ISVO4 = []
L2_ISVO5 = []
L2_ISVO6 = []
L2_ISVO7 = []
L2_ISVO8 = []
L2_ISVO9 = []
L2_ISV = []
L2_ISV1 = []
L2_ISV2 = []
L2_ISV3 = []
L2_ISV4 = []
L2_ISV5 = []
L2_ISV6 = []
L2_ISV7 = []
L2_ISV8 = []
L2_ISV9 = []


L1P_ISVO = []
L1P_ISVO1 = []
L1P_ISVO2 = []
L1P_ISVO3 = []
L1P_ISVO4 = []
L1P_ISVO5 = []
L1P_ISVO6 = []
L1P_ISVO7 = []
L1P_ISVO8 = []
L1P_ISVO9 = []
L1P_ISV = []
L1P_ISV1 = []
L1P_ISV2 = []
L1P_ISV3 = []
L1P_ISV4 = []
L1P_ISV5 = []
L1P_ISV6 = []
L1P_ISV7 = []
L1P_ISV8 = []
L1P_ISV9 = []


L2P_ISVO = []
L2P_ISVO1 = []
L2P_ISVO2 = []
L2P_ISVO3 = []
L2P_ISVO4 = []
L2P_ISVO5 = []
L2P_ISVO6 = []
L2P_ISVO7 = []
L2P_ISVO8 = []
L2P_ISVO9 = []
L2P_ISV = []
L2P_ISV1 = []
L2P_ISV2 = []
L2P_ISV3 = []
L2P_ISV4 = []
L2P_ISV5 = []
L2P_ISV6 = []
L2P_ISV7 = []
L2P_ISV8 = []
L2P_ISV9 = []


for j in range(0,TT-1,5):
#for j in range(TT-2,TT-1):
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
	ZH3 = []
	
	for i in range(0,N):
		if int(ZO[i][1]) == 1:
			ZH1.append(ZO[i][2])
		elif int(ZO[i][1]) == 2:
			ZH1.append(ZO[i][2])
		elif int(ZO[i][1]) == 3:
			ZH1.append(ZO[i][2])
		elif int(ZO[i][1]) == 4:
			ZH1.append(ZO[i][2])
		elif int(ZO[i][1]) == 5:
			ZH1.append(ZO[i][2])			
		elif int(ZO[i][1]) == 6:
			ZH1.append(ZO[i][2])	
		elif int(ZO[i][1]) == 7:
			ZH1.append(ZO[i][2])
		elif int(ZO[i][1]) == 8:
			ZH1.append(ZO[i][2])			
		elif int(ZO[i][1]) == 9:
			ZH1.append(ZO[i][2])	
		elif int(ZO[i][1]) == 10:
			ZH2.append(ZO[i][2])
		elif int(ZO[i][1]) == 11:
			ZH2.append(ZO[i][2])			
		elif int(ZO[i][1]) == 12:
			ZH2.append(ZO[i][2])	
		elif int(ZO[i][1]) == 13:
			ZH2.append(ZO[i][2])	
		elif int(ZO[i][1]) == 14:
			ZH2.append(ZO[i][2])
		elif int(ZO[i][1]) == 15:
			ZH2.append(ZO[i][2])			
		elif int(ZO[i][1]) == 16:
			ZH2.append(ZO[i][2])	
		elif int(ZO[i][1]) == 17:
			ZH2.append(ZO[i][2])
		elif int(ZO[i][1]) == 18:
			ZH2.append(ZO[i][2])			
			
			
	BL1 = 1.0
	BL2 = 2.0
	TL1 = 44.2
	TL2 = 45.2
	
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
	plt.hist(ZH1,bins = 100,color='m',label='ISVO')
	plt.hist(ZH2,bins = 100,color='r',alpha = 0.5,label='ISV')
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
	#plt.show()
	plt.close()
	
		
	#Compute distribution
	
	countO = 0
	countO1 = 0
	countO2 = 0
	countO3 = 0
	countO4 = 0
	countO5 = 0
	countO6 = 0
	countO7 = 0
	countO8 = 0
	countO9 = 0
	countI = 0
	countI1 = 0
	countI2 = 0
	countI3 = 0
	countI4 = 0
	countI5 = 0
	countI6 = 0
	countI7 = 0
	countI8 = 0
	countI9 = 0
	
	for i in BOT1:
		if i[1] == 1:
			countO = countO + 1	
			countO1 = countO1 + 1	
		elif i[1] == 2:
			countO = countO + 1	
			countO2 = countO2 + 1		
		elif i[1] == 3:
			countO = countO + 1	
			countO3 = countO3 + 1		
		elif i[1] == 4:
			countO = countO + 1	
			countO4 = countO4 + 1		
		elif i[1] == 5:
			countO = countO + 1	
			countO5 = countO5 + 1			
		elif i[1] == 6:
			countO = countO + 1	
			countO6 = countO6 + 1	
		elif i[1] == 7:
			countO = countO + 1	
			countO7 = countO7 + 1	
		elif i[1] == 8:
			countO = countO + 1	
			countO8 = countO8 + 1		
		elif i[1] == 9:
			countO = countO + 1	
			countO9 = countO9 + 1		
		elif i[1] == 10:
			countI = countI + 1	
			countI1 = countI1 + 1		
		elif i[1] == 11:
			countI = countI + 1	
			countI2 = countI2 + 1			
		elif i[1] == 12:
			countI = countI + 1	
			countI3 = countI3 + 1	
		elif i[1] == 13:
			countI = countI + 1	
			countI4 = countI4 + 1	
		elif i[1] == 14:
			countI = countI + 1	
			countI5 = countI5 + 1		
		elif i[1] == 15:
			countI = countI + 1	
			countI6 = countI6 + 1		
		elif i[1] == 16:
			countI = countI + 1	
			countI7 = countI7 + 1		
		elif i[1] == 17:
			countI = countI + 1	
			countI8 = countI8 + 1			
		elif i[1] == 18:
			countI = countI + 1	
			countI9 = countI9 + 1				
	for i in TOP1:
		if i[1] == 1:
			countO = countO + 1	
			countO1 = countO1 + 1	
		elif i[1] == 2:
			countO = countO + 1	
			countO2 = countO2 + 1		
		elif i[1] == 3:
			countO = countO + 1	
			countO3 = countO3 + 1		
		elif i[1] == 4:
			countO = countO + 1	
			countO4 = countO4 + 1		
		elif i[1] == 5:
			countO = countO + 1	
			countO5 = countO5 + 1			
		elif i[1] == 6:
			countO = countO + 1	
			countO6 = countO6 + 1	
		elif i[1] == 7:
			countO = countO + 1	
			countO7 = countO7 + 1	
		elif i[1] == 8:
			countO = countO + 1	
			countO8 = countO8 + 1		
		elif i[1] == 9:
			countO = countO + 1	
			countO9 = countO9 + 1		
		elif i[1] == 10:
			countI = countI + 1	
			countI1 = countI1 + 1		
		elif i[1] == 11:
			countI = countI + 1	
			countI2 = countI2 + 1			
		elif i[1] == 12:
			countI = countI + 1	
			countI3 = countI3 + 1	
		elif i[1] == 13:
			countI = countI + 1	
			countI4 = countI4 + 1	
		elif i[1] == 14:
			countI = countI + 1	
			countI5 = countI5 + 1		
		elif i[1] == 15:
			countI = countI + 1	
			countI6 = countI6 + 1		
		elif i[1] == 16:
			countI = countI + 1	
			countI7 = countI7 + 1		
		elif i[1] == 17:
			countI = countI + 1	
			countI8 = countI8 + 1			
		elif i[1] == 18:
			countI = countI + 1	
			countI9 = countI9 + 1		

	L1_ISVO.append(countO)	
	L1_ISV.append(countI)
	
	L1_ISVO1.append(countO1)	
	L1_ISV1.append(countI1)

	L1_ISVO2.append(countO2)	
	L1_ISV2.append(countI2)	
	
	L1_ISVO3.append(countO3)	
	L1_ISV3.append(countI3)	

	L1_ISVO4.append(countO4)	
	L1_ISV4.append(countI4)

	L1_ISVO5.append(countO5)	
	L1_ISV5.append(countI5)	
	
	L1_ISVO6.append(countO6)	
	L1_ISV6.append(countI6)

	L1_ISVO7.append(countO7)	
	L1_ISV7.append(countI7)

	L1_ISVO8.append(countO8)	
	L1_ISV8.append(countI8)	
	
	L1_ISVO9.append(countO9)	
	L1_ISV9.append(countI9)	
	
	
	L1P_ISVO1.append(100*countO1/(countO+countI))
	L1P_ISV1.append(100*countI1/(countO+countI))

	L1P_ISVO2.append(100*countO2/(countO+countI))
	L1P_ISV2.append(100*countI2/(countO+countI))

	L1P_ISVO3.append(100*countO3/(countO+countI))
	L1P_ISV3.append(100*countI3/(countO+countI))	
	
	L1P_ISVO4.append(100*countO4/(countO+countI))
	L1P_ISV4.append(100*countI4/(countO+countI))

	L1P_ISVO5.append(100*countO5/(countO+countI))
	L1P_ISV5.append(100*countI5/(countO+countI))

	L1P_ISVO6.append(100*countO6/(countO+countI))
	L1P_ISV6.append(100*countI6/(countO+countI))	

	L1P_ISVO7.append(100*countO7/(countO+countI))
	L1P_ISV7.append(100*countI7/(countO+countI))

	L1P_ISVO8.append(100*countO8/(countO+countI))
	L1P_ISV8.append(100*countI8/(countO+countI))

	L1P_ISVO9.append(100*countO9/(countO+countI))
	L1P_ISV9.append(100*countI9/(countO+countI))	
	
	L1P_ISVO.append(100*countO/(countO+countI))
	L1P_ISV.append(100*countI/(countO+countI))	


	countO = 0
	countO1 = 0
	countO2 = 0
	countO3 = 0
	countO4 = 0
	countO5 = 0
	countO6 = 0
	countO7 = 0
	countO8 = 0
	countO9 = 0
	countI = 0
	countI1 = 0
	countI2 = 0
	countI3 = 0
	countI4 = 0
	countI5 = 0
	countI6 = 0
	countI7 = 0
	countI8 = 0
	countI9 = 0
	
	for i in BOT2:
		if i[1] == 1:
			countO = countO + 1	
			countO1 = countO1 + 1	
		elif i[1] == 2:
			countO = countO + 1	
			countO2 = countO2 + 1		
		elif i[1] == 3:
			countO = countO + 1	
			countO3 = countO3 + 1		
		elif i[1] == 4:
			countO = countO + 1	
			countO4 = countO4 + 1		
		elif i[1] == 5:
			countO = countO + 1	
			countO5 = countO5 + 1			
		elif i[1] == 6:
			countO = countO + 1	
			countO6 = countO6 + 1	
		elif i[1] == 7:
			countO = countO + 1	
			countO7 = countO7 + 1	
		elif i[1] == 8:
			countO = countO + 1	
			countO8 = countO8 + 1		
		elif i[1] == 9:
			countO = countO + 1	
			countO9 = countO9 + 1		
		elif i[1] == 10:
			countI = countI + 1	
			countI1 = countI1 + 1		
		elif i[1] == 11:
			countI = countI + 1	
			countI2 = countI2 + 1			
		elif i[1] == 12:
			countI = countI + 1	
			countI3 = countI3 + 1	
		elif i[1] == 13:
			countI = countI + 1	
			countI4 = countI4 + 1	
		elif i[1] == 14:
			countI = countI + 1	
			countI5 = countI5 + 1		
		elif i[1] == 15:
			countI = countI + 1	
			countI6 = countI6 + 1		
		elif i[1] == 16:
			countI = countI + 1	
			countI7 = countI7 + 1		
		elif i[1] == 17:
			countI = countI + 1	
			countI8 = countI8 + 1			
		elif i[1] == 18:
			countI = countI + 1	
			countI9 = countI9 + 1				
	for i in TOP2:
		if i[1] == 1:
			countO = countO + 1	
			countO1 = countO1 + 1	
		elif i[1] == 2:
			countO = countO + 1	
			countO2 = countO2 + 1		
		elif i[1] == 3:
			countO = countO + 1	
			countO3 = countO3 + 1		
		elif i[1] == 4:
			countO = countO + 1	
			countO4 = countO4 + 1		
		elif i[1] == 5:
			countO = countO + 1	
			countO5 = countO5 + 1			
		elif i[1] == 6:
			countO = countO + 1	
			countO6 = countO6 + 1	
		elif i[1] == 7:
			countO = countO + 1	
			countO7 = countO7 + 1	
		elif i[1] == 8:
			countO = countO + 1	
			countO8 = countO8 + 1		
		elif i[1] == 9:
			countO = countO + 1	
			countO9 = countO9 + 1		
		elif i[1] == 10:
			countI = countI + 1	
			countI1 = countI1 + 1		
		elif i[1] == 11:
			countI = countI + 1	
			countI2 = countI2 + 1			
		elif i[1] == 12:
			countI = countI + 1	
			countI3 = countI3 + 1	
		elif i[1] == 13:
			countI = countI + 1	
			countI4 = countI4 + 1	
		elif i[1] == 14:
			countI = countI + 1	
			countI5 = countI5 + 1		
		elif i[1] == 15:
			countI = countI + 1	
			countI6 = countI6 + 1		
		elif i[1] == 16:
			countI = countI + 1	
			countI7 = countI7 + 1		
		elif i[1] == 17:
			countI = countI + 1	
			countI8 = countI8 + 1			
		elif i[1] == 18:
			countI = countI + 1	
			countI9 = countI9 + 1		

	L2_ISVO.append(countO)	
	L2_ISV.append(countI)

	L2_ISVO1.append(countO1)	
	L2_ISV1.append(countI1)

	L2_ISVO2.append(countO2)	
	L2_ISV2.append(countI2)	

	L2_ISVO3.append(countO3)	
	L2_ISV3.append(countI3)	

	L2_ISVO4.append(countO4)	
	L2_ISV4.append(countI4)

	L2_ISVO5.append(countO5)	
	L2_ISV5.append(countI5)	

	L2_ISVO6.append(countO6)	
	L2_ISV6.append(countI6)

	L2_ISVO7.append(countO7)	
	L2_ISV7.append(countI7)

	L2_ISVO8.append(countO8)	
	L2_ISV8.append(countI8)	

	L2_ISVO9.append(countO9)	
	L2_ISV9.append(countI9)	
	
	
	L2P_ISVO1.append(100*countO1/(countO+countI))
	L2P_ISV1.append(100*countI1/(countO+countI))

	L2P_ISVO2.append(100*countO2/(countO+countI))
	L2P_ISV2.append(100*countI2/(countO+countI))

	L2P_ISVO3.append(100*countO3/(countO+countI))
	L2P_ISV3.append(100*countI3/(countO+countI))	

	L2P_ISVO4.append(100*countO4/(countO+countI))
	L2P_ISV4.append(100*countI4/(countO+countI))

	L2P_ISVO5.append(100*countO5/(countO+countI))
	L2P_ISV5.append(100*countI5/(countO+countI))

	L2P_ISVO6.append(100*countO6/(countO+countI))
	L2P_ISV6.append(100*countI6/(countO+countI))	

	L2P_ISVO7.append(100*countO7/(countO+countI))
	L2P_ISV7.append(100*countI7/(countO+countI))

	L2P_ISVO8.append(100*countO8/(countO+countI))
	L2P_ISV8.append(100*countI8/(countO+countI))

	L2P_ISVO9.append(100*countO9/(countO+countI))
	L2P_ISV9.append(100*countI9/(countO+countI))	

	L2P_ISVO.append(100*countO/(countO+countI))
	L2P_ISV.append(100*countI/(countO+countI))	

	
	
	
	
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
	G1.write('18 atom types \n')
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
	G1.write('4 1.0 \n')
	G1.write('5 1.0 \n')
	G1.write('6 1.0 \n')
	G1.write('7 1.0 \n')
	G1.write('8 1.0 \n')
	G1.write('9 1.0 \n')
	G1.write('10 1.0 \n')
	G1.write('11 1.0 \n')
	G1.write('12 1.0 \n')
	G1.write('13 1.0 \n')
	G1.write('14 1.0 \n')
	G1.write('15 1.0 \n')
	G1.write('16 1.0 \n')
	G1.write('17 1.0 \n')
	G1.write('18 1.0 \n')

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
	G2.write('18 atom types \n')
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
	G2.write('4 1.0 \n')
	G2.write('5 1.0 \n')
	G2.write('6 1.0 \n')
	G2.write('7 1.0 \n')
	G2.write('8 1.0 \n')
	G2.write('9 1.0 \n')
	G2.write('10 1.0 \n')
	G2.write('11 1.0 \n')
	G2.write('12 1.0 \n')
	G2.write('13 1.0 \n')
	G2.write('14 1.0 \n')
	G2.write('15 1.0 \n')
	G2.write('16 1.0 \n')
	G2.write('17 1.0 \n')
	G2.write('18 1.0 \n')
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
plt.plot(L1_ISV,color='r',linewidth=2.0,label='ISV')
plt.legend()
plt.title('Distribution of Particles in Layer 1')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L1.png',bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(L1_ISVO1,linewidth=2.0,label='ISVO1')
plt.plot(L1_ISVO2,linewidth=2.0,label='ISVO2')
plt.plot(L1_ISVO3,linewidth=2.0,label='ISVO3')
plt.plot(L1_ISVO4,linewidth=2.0,label='ISVO4')
plt.plot(L1_ISVO5,linewidth=2.0,label='ISVO5')
plt.plot(L1_ISVO6,linewidth=2.0,label='ISVO6')
plt.plot(L1_ISVO7,linewidth=2.0,label='ISVO7')
plt.plot(L1_ISVO7,linewidth=2.0,label='ISVO8')
plt.plot(L1_ISVO8,linewidth=2.0,label='ISVO9')
plt.legend()
plt.title('Distribution of ISVO Particles in Layer 1')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L1ISVO.png',bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(L1_ISV1,linewidth=2.0,label='ISV1')
plt.plot(L1_ISV2,linewidth=2.0,label='ISV2')
plt.plot(L1_ISV3,linewidth=2.0,label='ISV3')
plt.plot(L1_ISV4,linewidth=2.0,label='ISV4')
plt.plot(L1_ISV5,linewidth=2.0,label='ISV5')
plt.plot(L1_ISV6,linewidth=2.0,label='ISV6')
plt.plot(L1_ISV7,linewidth=2.0,label='ISV7')
plt.plot(L1_ISV8,linewidth=2.0,label='ISV8')
plt.plot(L1_ISV9,linewidth=2.0,label='ISV9')
plt.legend()
plt.title('Distribution of ISV Particles in Layer 1')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L1ISV.png',bbox_inches='tight')
plt.close()


plt.figure()
plt.plot(L2_ISVO1,linewidth=2.0,label='ISVO1')
plt.plot(L2_ISVO2,linewidth=2.0,label='ISVO2')
plt.plot(L2_ISVO3,linewidth=2.0,label='ISVO3')
plt.plot(L2_ISVO4,linewidth=2.0,label='ISVO4')
plt.plot(L2_ISVO5,linewidth=2.0,label='ISVO5')
plt.plot(L2_ISVO6,linewidth=2.0,label='ISVO6')
plt.plot(L2_ISVO7,linewidth=2.0,label='ISVO7')
plt.plot(L2_ISVO8,linewidth=2.0,label='ISVO8')
plt.plot(L2_ISVO9,linewidth=2.0,label='ISVO9')
plt.legend()
plt.title('Distribution of ISVO Particles in Layer 2')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L2ISVO.png',bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(L2_ISV1,linewidth=2.0,label='ISV1')
plt.plot(L2_ISV2,linewidth=2.0,label='ISV2')
plt.plot(L2_ISV3,linewidth=2.0,label='ISV3')
plt.plot(L2_ISV4,linewidth=2.0,label='ISV4')
plt.plot(L2_ISV5,linewidth=2.0,label='ISV5')
plt.plot(L2_ISV6,linewidth=2.0,label='ISV6')
plt.plot(L2_ISV7,linewidth=2.0,label='ISV7')
plt.plot(L2_ISV8,linewidth=2.0,label='ISV8')
plt.plot(L2_ISV9,linewidth=2.0,label='ISV9')
plt.legend()
plt.title('Distribution of ISV Particles in Layer 2')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L2ISV.png',bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(L1P_ISVO1,linewidth=2.0,label='ISVO1')
plt.plot(L1P_ISVO2,linewidth=2.0,label='ISVO2')
plt.plot(L1P_ISVO3,linewidth=2.0,label='ISVO3')
plt.plot(L1P_ISVO4,linewidth=2.0,label='ISVO4')
plt.plot(L1P_ISVO5,linewidth=2.0,label='ISVO5')
plt.plot(L1P_ISVO6,linewidth=2.0,label='ISVO6')
plt.plot(L1P_ISVO7,linewidth=2.0,label='ISVO7')
plt.plot(L1P_ISVO8,linewidth=2.0,label='ISVO8')
plt.plot(L1P_ISVO9,linewidth=2.0,label='ISVO9')
plt.legend()
plt.title('Percentage Distribution of ISVO Particles in Layer 1')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L1ISVO.png',bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(L1P_ISV1,linewidth=2.0,label='ISV1')
plt.plot(L1P_ISV2,linewidth=2.0,label='ISV2')
plt.plot(L1P_ISV3,linewidth=2.0,label='ISV3')
plt.plot(L1P_ISV4,linewidth=2.0,label='ISV4')
plt.plot(L1P_ISV5,linewidth=2.0,label='ISV5')
plt.plot(L1P_ISV6,linewidth=2.0,label='ISV6')
plt.plot(L1P_ISV7,linewidth=2.0,label='ISV7')
plt.plot(L1P_ISV8,linewidth=2.0,label='ISV8')
plt.plot(L1P_ISV9,linewidth=2.0,label='ISV9')
plt.legend()
plt.title('Percentage Distribution of ISV Particles in Layer 1')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L1ISV.png',bbox_inches='tight')
plt.close()


plt.figure()
plt.plot(L2P_ISVO1,linewidth=2.0,label='ISVO1')
plt.plot(L2P_ISVO2,linewidth=2.0,label='ISVO2')
plt.plot(L2P_ISVO3,linewidth=2.0,label='ISVO3')
plt.plot(L2P_ISVO4,linewidth=2.0,label='ISVO4')
plt.plot(L2P_ISVO5,linewidth=2.0,label='ISVO5')
plt.plot(L2P_ISVO6,linewidth=2.0,label='ISVO6')
plt.plot(L2P_ISVO7,linewidth=2.0,label='ISVO7')
plt.plot(L2P_ISVO8,linewidth=2.0,label='ISVO8')
plt.plot(L2P_ISVO9,linewidth=2.0,label='ISVO9')
plt.legend()
plt.title('Percentage Distribution of ISVO Particles in Layer 2')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L2ISVO.png',bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(L2P_ISV1,linewidth=2.0,label='ISV1')
plt.plot(L2P_ISV2,linewidth=2.0,label='ISV2')
plt.plot(L2P_ISV3,linewidth=2.0,label='ISV3')
plt.plot(L2P_ISV4,linewidth=2.0,label='ISV4')
plt.plot(L2P_ISV5,linewidth=2.0,label='ISV5')
plt.plot(L2P_ISV6,linewidth=2.0,label='ISV6')
plt.plot(L2P_ISV7,linewidth=2.0,label='ISV7')
plt.plot(L2P_ISV8,linewidth=2.0,label='ISV8')
plt.plot(L2P_ISV9,linewidth=2.0,label='ISV9')
plt.legend()
plt.title('Percentage Distribution of ISV Particles in Layer 2')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L2ISV.png',bbox_inches='tight')
plt.close()


plt.figure()
plt.plot(L2_ISVO,color='k',linewidth=2.0,label='ISVO')
plt.plot(L2_ISV,color='r',linewidth=2.0,label='ISV')
plt.legend()
plt.title('Distribution of Particles in Layer 2')
plt.ylabel('Number of Particles')
plt.xlabel('Time')
plt.savefig('TIME/L2.png',bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(L1P_ISVO,color='k',linewidth=2.0,label='ISVO')
plt.plot(L1P_ISV,color='r',linewidth=2.0,label='ISV')
plt.legend()
plt.title('Percentage Distribution of Particles in Layer 1')
plt.ylabel('Percentage')
plt.xlabel('Time')
plt.savefig('TIME/L1P.png',bbox_inches='tight')
plt.close()

plt.figure()
plt.plot(L2P_ISVO,color='k',linewidth=2.0,label='ISVO')
plt.plot(L2P_ISV,color='r',linewidth=2.0,label='ISV')
plt.legend()
plt.title('Percentage Distribution of Particles in Layer 2')
plt.ylabel('Percentage')
plt.xlabel('Time')
plt.savefig('TIME/L2P.png',bbox_inches='tight')
plt.close()



G1 = open('TIME/L1_ISVOdiv.dat','w')
for i in range(0,len(L1_ISVO1)):
	G1.write(str(L1_ISVO1[i]) + ' ' + str(L1_ISVO2[i]) + ' ' + str(L1_ISVO3[i]) + ' '  +
			str(L1_ISVO4[i]) + ' ' + str(L1_ISVO5[i]) + ' ' + str(L1_ISVO6[i]) + ' ' +
			str(L1_ISVO7[i]) + ' ' + str(L1_ISVO8[i]) + ' ' + str(L1_ISVO9[i]) + ' ' + '/n')
G1.close()

G1 = open('TIME/L1_ISVdiv.dat','w')
for i in range(0,len(L1_ISV1)):
	G1.write(str(L1_ISV1[i]) + ' ' + str(L1_ISV2[i]) + ' ' + str(L1_ISV3[i]) + ' ' +  
			str(L1_ISV4[i]) + ' ' + str(L1_ISV5[i]) + ' ' + str(L1_ISV6[i]) + ' '  +
			str(L1_ISV7[i]) + ' ' + str(L1_ISV8[i]) + ' ' + str(L1_ISV9[i]) + ' ' + '/n')
G1.close()


G1 = open('TIME/L2_ISVOdiv.dat','w')
for i in range(0,len(L2_ISVO1)):
	G1.write(str(L2_ISVO1[i]) + ' ' + str(L2_ISVO2[i]) + ' ' + str(L2_ISVO3[i]) + ' ' + 
			str(L2_ISVO4[i]) + ' ' + str(L2_ISVO5[i]) + ' ' + str(L2_ISVO6[i]) + ' ' +
			str(L2_ISVO7[i]) + ' ' + str(L2_ISVO8[i]) + ' ' + str(L2_ISVO9[i]) + ' ' + '/n')
G1.close()

G1 = open('TIME/L2_ISVdiv.dat','w')
for i in range(0,len(L2_ISV1)):
	G1.write(str(L2_ISV1[i]) + ' ' + str(L2_ISV2[i]) + ' ' + str(L2_ISV3[i]) + ' '  +
			str(L2_ISV4[i]) + ' ' + str(L2_ISV5[i]) + ' ' + str(L2_ISV6[i]) + ' ' +
			str(L2_ISV7[i]) + ' ' + str(L2_ISV8[i]) + ' ' + str(L2_ISV9[i]) + ' ' + '/n')
G1.close()


G1 = open('TIME/L1_ISVO.dat','w')
for i in L1_ISVO:
	G1.write(str(i) + '/n')
G1.close()

G1 = open('TIME/L2_ISVO.dat','w')
for i in L2_ISVO:
	G1.write(str(i) + '/n')
G1.close()

G1 = open('TIME/L1P_ISVO.dat','w')
for i in L1P_ISVO:
	G1.write(str(i) + '/n')
G1.close()

G1 = open('TIME/L2P_ISVO.dat','w')
for i in L2P_ISVO:
	G1.write(str(i) + '/n')
G1.close()

G1 = open('TIME/L1_ISV.dat','w')
for i in L1_ISV:
	G1.write(str(i) + '/n')
G1.close()

G1 = open('TIME/L2_ISV.dat','w')
for i in L2_ISV:
	G1.write(str(i) + '/n')
G1.close()

G1 = open('TIME/L1P_ISV.dat','w')
for i in L1P_ISV:
	G1.write(str(i) + '/n')
G1.close()

G1 = open('TIME/L2P_ISV.dat','w')
for i in L2P_ISV:
	G1.write(str(i) + '/n')
G1.close()
