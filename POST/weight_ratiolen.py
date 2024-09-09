###########################UPDATED 3/3/2022 ######################################
############# Prajwal Bangalore Prakash ##################
# Minor changes Luis Nieves March 4th, 2022
# Minor changes Luis Nieves March 10th, 2022
# Minor changes Luis Nieves March 21st, 2022
import numpy as np
import matplotlib.pyplot as plt
import sys as s
import os
import glob
import math
import random 
import statistics
from PIL import Image
from heapq import nsmallest
#import cv2

FILEIN = s.argv[1] # LN34
OUTN = s.argv[2] # LN34
CUTOFF = 30 # LN34

OUTDATA = [] # LN34
files=[FILEIN] # LN34
flip=[]
x_1=[]
y_1=[]
z_1=[]
c= 0
for g in range(len(files)):
	f1=open(files[g])
  
	for line in f1:
		c= c+1
		x_1.append(float(line.split()[0]))
		y_1.append(float(line.split()[1]))
		z_1.append(float(line.split()[2]))

L= float(x_1[0])
HN= float(x_1[0])/2.0
#print(HN, 'HN') # LN34
HNx=[]
HNy=[]

HNx.append(HN)
HNy.append(HN)

HNx.append(-HN)
HNy.append(HN)

HNx.append(-HN)
HNy.append(-HN)

HNx.append(HN)
HNy.append(-HN)

HNx.append(HN)
HNy.append(HN)

#print(HNx, HNy) # LN34

files=['coordvoron.inp']
flip=[]
x=[]
y=[]
z=[]
c= 0
for g in range(len(files)):
	f1=open(files[g])
  
	for line in f1:
		c= c+1
		x.append(float(line.split()[0]))
		y.append(float(line.split()[1]))
		z.append(float(line.split()[2]))
#print(len(x), 'N value') # LN34
files=['neighborlist.dat']
flip=[]
s1=[]
s2=[]
s3=[]
c= 0
for g in range(len(files)):
	f1=open(files[g])
  
	for line in f1:
		c= c+1
		s1.append((line.split()[0]))
		s2.append((line.split()[1]))
		s3.append(float(line.split()[2]))

count, bins, ignored = plt.hist(s3, 25, edgecolor='black', linewidth=1.2)
plt.tick_params(labelsize=20)
#plt.savefig('edgelength.png',bboxsize='tight') # LN34

files=['neighborlistpoints.dat']
flip=[]
x1=[]
x2=[]
x3=[]
c= 0
for g in range(len(files)):
	f1=open(files[g])
  
	for line in f1:
		c= c+1
		x1.append((line.split()[0]))
		x2.append(float(line.split()[1]))
		x3.append(float(line.split()[2]))

#print(len((s3)), 'len of s') # LN34
#print(len(x1), 'lenx1') # LN34
#count1=0
#s1new=[]
#s2new=[]
#s3new=[]
#for i in range(len(s1)):
#	if(s3[i]>10.0):
#		count1=count1+1
#		s1new.append(s1[i])
#		s2new.append(s2[i])
#		s3new.append(s3[i])

#print(s3new[0], 's1new')
#count, bins, ignored = plt.hist(s3new, 25, edgecolor='black', linewidth=1.2)
#plt.show() # LN


files=['neighnum.dat']
neighnum=[]
for g in range(len(files)):
	f1=open(files[g])
  
	for line in f1:
		c= c+1
		neighnum.append((line.split()[0]))

bins_list = [2, 3, 4, 5, 6, 7, 8, 9, 10]
plt.figure()
plt.bar(*np.unique(neighnum, return_counts=True))
plt.title("Number of sides")
plt.tick_params(labelsize=20)
#plt.savefig('neighnum.png',bboxsize='tight')# LN34
#plt.show() # LN34



#print(x3, 'x3')

f=open("ratiolen.dat",'w')
f1=open("angle.dat",'w')
N_list = list(range(1, len(x)+1))
#print(N_list, 'N_list')
lenrat=[]
angle=[]
neighnew_proj=[]
neighnew_len=[]

figure = plt.figure()
figure.set_size_inches(20, 20)
ax = plt.gca()
ax.scatter(x, y, s=10, c="k")
c4=0
c6=0
cdef=0
dist=[]
dist_def=[]
dist_ISV=[]
dist_SVPS=[]
num_def=[]
part1_len_4=[]
partf1=[]

partf2=[]
f15=open("ISVneigh_"+OUTN+".dat",'w')     #LN310
f19=open("ISV_"+OUTN+".dat",'w')          #LN310
f16=open("SVPSneigh_"+OUTN+".dat",'w')    #LN310
f20=open("SVPS_"+OUTN+".dat",'w')         #LN310
f17=open("defect_"+OUTN+".dat",'w')       #LN310
f18=open("defectneigh_"+OUTN+".dat",'w')  #LN310
part1n4=[]
part2n4=[]

part1n6=[]
part2n6=[]
dist_4=[]
part1ndef=[]
part2ndef=[]

for k in range(len(x)):
    xvert=[]
    yvert=[]
    partneigh1=[]
    partneigh2=[]
    c= 0
    cnew_proj=0
    cnew_len=0
    part1=[]
    part2=[]
    dist_check=[]
    for i in range(len(s1)):
        
        #part1=[]
        #part2=[]
#        if(int(s1[i])==int(N_list[k])):
		#part1=[]
		#part2=[]
        if(int(s1[i])==int(N_list[k])):
#            print(s2[i])
            c=c+1
            #print(c, neighnum[1], i)
#            if(int(c)==int(neighnum[k])):
            if(int(c)==int(neighnum[k])):
                #print(s3[i], s3[i+1-int(neighnum[k])], 'first and last')
                ratio1= float(s3[i])/float(s3[i+1-int(neighnum[k])])
                #print(ratio1, 'ratio1')
                ratio2= float(1)/float(ratio1)
                
                vect1x= x2[i]-x[k]
                vect1y= x3[i]-y[k]
                denom1=np.sqrt(vect1x**2.0+vect1y**2.0)
                vect1x= float(x2[i]-x[k])/denom1
                vect1y= float(x3[i]-y[k])/denom1
                vect2x= x2[i+1-int(neighnum[k])]-x[k]
                vect2y= x3[i+1-int(neighnum[k])]-y[k]
                denom2=np.sqrt(vect2x**2.0+vect2y**2.0)
                vect2x= float(x2[i+1-int(neighnum[k])]-x[k])/denom2
                vect2y= float(x3[i+1-int(neighnum[k])]-y[k])/denom2
                
                proj= (vect1x*vect2x)+(vect1y*vect2y)
                #proj=0.0
                angle.append(math.degrees(math.acos((proj))))
                ne=int(s2[i])
                DX= (x[k]-x[ne-1])
                DY= (y[k]-y[ne-1])
                #DX= DX- L*round(float(DX/L))
                #DY= DY- L*round(float(DY/L))
                var= (DX)**2.0+(DY)**2.0
                dist_check.append(np.sqrt(var))
                part1.append(k+1)
                part2.append(ne)
                if(math.degrees(math.acos((proj)))>CUTOFF): # LN34
                    #if(math.degrees(math.acos((proj)))<100):
                    #ne= int(s2[i])
                    #print(ne)
                    #DX= (x[k]-x[ne-1])
                    #DY= (y[k]-y[ne-1])
                    #DX= DX- L*round(float(DX/L))
                    #DY= DY- L*round(float(DY/L))
                    #var= (DX)**2.0+(DY)**2.0
                    #dist.append(np.sqrt(var))
                    #part1.append(k+1)
                    #part2.append(ne)
                    partneigh1.append(int(k)+1)
                    partneigh2.append(ne)
                    cnew_proj=cnew_proj+1
                    xvert.append(x2[i])
                    yvert.append(x3[i])
                
                #print(ratio2, 'ratio2')
                ratio= min(ratio1,ratio2)
                
                if(ratio> 0.10):
                    cnew_len=cnew_len+1
                    
                f.write("%04.18e\n"%(ratio))
                f1.write("%04.18e\n"%(math.degrees(math.acos((proj)))))
                lenrat.append(ratio)				
            else:
                ratio1= float(s3[i])/float(s3[i+1])
                ratio2= float(1.0)/float(ratio1)
                
                vect1x= x2[i]-x[k]
                vect1y= x3[i]-y[k]
                denom1=np.sqrt(vect1x**2.0+vect1y**2.0)
                vect1x= float(x2[i]-x[k])/denom1
                vect1y= float(x3[i]-y[k])/denom1
                vect2x= x2[i+1]-x[k]
                vect2y= x3[i+1]-y[k]
                denom2=np.sqrt(vect2x**2.0+vect2y**2.0)
                vect2x= float(x2[i+1]-x[k])/denom2
                vect2y= float(x3[i+1]-y[k])/denom2
                
                proj= (vect1x*vect2x)+(vect1y*vect2y)
                #proj=0.0
                angle.append(math.degrees(math.acos((proj))))
                ne=int(s2[i])
                DX= (x[k]-x[ne-1])
                DY= (y[k]-y[ne-1])
                #DX= DX- L*round(float(DX/L))
                #DY= DY- L*round(float(DY/L))
                var= (DX)**2.0+(DY)**2.0
                dist_check.append(np.sqrt(var))
                part1.append(k+1)
                part2.append(ne)
                if(math.degrees(math.acos((proj)))>CUTOFF): # LN34
                    #if(math.degrees(math.acos((proj)))<100):
                    #ne= int(s2[i])
                    #print(ne)
                    #DX= (x[k]-x[ne-1])
                    #DY= (y[k]-y[ne-1])
                    #DX= DX- L*round(float(DX/L))
                    #DY= DY- L*round(float(DY/L))
                    #var= (DX)**2.0+(DY)**2.0
                    #dist.append(np.sqrt(var))
                    #part1.append(k+1)
                    #part2.append(ne)
                    #print(part1)
                    #print(part2)
                    partneigh1.append(int(k)+1)
                    partneigh2.append(ne)
                    cnew_proj=cnew_proj+1
                    xvert.append(x2[i])
                    yvert.append(x3[i])
                
                ratio= min(ratio1,ratio2)
                
                if(ratio> 0.10):
                    cnew_len=cnew_len+1
                f.write("%04.18e\n"%(ratio))
                f1.write("%04.18e\n"%(math.degrees(math.acos((proj)))))
                lenrat.append(ratio)
    neighnew_proj.append(cnew_proj)
    neighnew_len.append(cnew_len)
    if(int(cnew_proj)==int(4)):
		#part1_len_4.append(len(part1))
#		for m in range(len(part1)):
#			DX= (x[part1[m]]-x[part2[m]])
#			DY= (y[part1[m]]-y[part2[m]])
#			DX= DX- L*round(float(DX/L))
#			DY= DY- L*round(float(DY/L))
#			dist_ISV.append(np.sqrt((DX)**2.0+(DY)**2.0))
        ax.fill(xvert, yvert, facecolor='red', edgecolor='black', linewidth=3.0, alpha= 0.4)
        OUTDATA.append([x[k]+L/2,y[k]+L/2,1]) # LN34
        c4=c4+1
        #dist_4.append(nsmallest(3, dist_check)[0])
        #dist_4.append(nsmallest(3, dist_check)[1])
        #dist_4.append(nsmallest(3, dist_check)[2])
        #dist_4.append(nsmallest(3, dist_check)[3])
        f19.write("%s\t%s\n"%(int(k)+1,int(N_list[k])))
        #part1n4.append()
        #print(len(dist_check), 'dist')
        #print('hello')
        #idx=np.argpartition(dist_check, int(cnew_proj))
        idx= np.argsort(dist_check)[:int(cnew_proj)]
        #print('hello')
        #part1n4.append(part1[idx[0]])
        #part1n4.append(part1[idx[0]])
        for m in range(int(cnew_proj)):
            #print(m)
            #print(part1[m])
            #print(part2[m])
            DX= (x[part1[m]-1]-x[part2[m]-1])
            DY= (y[part1[m]-1]-y[part2[m]-1])
            #DX= DX - L*round(float(DX/L))
            #DY= DY - L*round(float(DY/L))
            part1n4.append(part1[idx[m]])
            part2n4.append(part2[idx[m]])
            dist_ISV.append(np.sqrt((DX)**2.0+(DY)**2.0))
            #if(np.sqrt((DX)**2.0+(DY)**2.0)>60):
            #if(np.sqrt((DX)**2.0+(DY)**2.0)>50 and np.sqrt((DX)**2.0+(DY)**2.0)<60):
            #    ax.fill(xvert, yvert, facecolor='red', edgecolor='black', linewidth=3.0)
                
        #print(part2, 'part2')
        #print(part1, 'part1')
        #for m in range(len(part1)):
		#	DX= (x[part1[m]]-x[part2[m]])
			
        #part1_len_4.append(len(part2))
    elif(int(cnew_proj)==int(6)):
        #ax.scatter(x[k], y[k], s=10, c="k")
        ax.fill(xvert, yvert, facecolor='blue', edgecolor='black', linewidth=3.0, alpha= 0.4)
        OUTDATA.append([x[k]+L/2,y[k]+L/2,2]) # LN34
        c6=c6+1
        f20.write("%s\t%s\n"%(int(k)+1,int(N_list[k])))
        idx= np.argsort(dist_check)[:int(cnew_proj)]
        for m in range(int(cnew_proj)):
            DX= (x[part1[m]-1]-x[part2[m]-1])
            DY= (y[part1[m]-1]-y[part2[m]-1])
            #DX= DX - L*round(float(DX/L))
            #DY= DY - L*round(float(DY/L))
            part1n6.append(part1[idx[m]])
            part2n6.append(part2[idx[m]])
            dist_SVPS.append(np.sqrt((DX)**2.0+(DY)**2.0))
        #ax.fill(xvert, yvert, facecolor='none', edgecolor='black', linewidth=1)
    else:
	
		#for m in range(len(part1)):
        #dist_def.append(dist)
        num_def.append(int(cnew_proj))
        #partf1.append(part1)
        #partf2.append(part2)
        f17.write("%s\t%s\n"%(int(k)+1,int(N_list[k])))
        idx= np.argsort(dist_check)[:int(cnew_proj)]
        for m in range(int(cnew_proj)):
            DX= (x[part1[m]-1]-x[part2[m]-1])
            DX= (y[part1[m]-1]-y[part2[m]-1])
            #DX= DX - L*round(float(DX/L))
            #DY= DY - L*round(float(DY/L))
            part1ndef.append(part1[idx[m]])
            part2ndef.append(part2[idx[m]])
            dist_def.append(np.sqrt((DX)**2.0+(DY)**2.0))
        #ax.scatter(x[k], y[k], s=10, c="k")
        ax.fill(xvert, yvert, facecolor='none', edgecolor='black', linewidth=3.0, alpha= 0.4)
        OUTDATA.append([x[k]+L/2,y[k]+L/2,3]) # LN
        cdef=cdef+1
        #ax.fill(xvert, yvert, facecolor='none', edgecolor='black', linewidth=1)


#ax.scatter(x, y, s=10, c="k")
plt.plot(HNx,HNy,'-k', linewidth= 5.0)
#plt.savefig('Modified.png',bboxsize='tight')
plt.axis('off')
plt.yticks(rotation=180)
plt.xlim(-HN,HN)
plt.ylim(-HN,HN)
plt.tick_params(labelsize=30)
#cover = resizeimage.resize_cover(image, [200, 200])
#cover.save('voronoifinal.png', dpi=1000)
plt.savefig('voronoifinal'+OUTN+ '.png', dpi=600, bboxsize='tight') # LN
#plt.show() # LN34

# LN34
FOUT = open('OUTDATA'+OUTN+'.dat','w')
for i in OUTDATA:
	FOUT.write(str(i[0]) + ' ' + str(i[1]) + ' ' + str(i[2]) + '\n')
FOUT.close()	

# LN34
FOUT2 = open('GENDATA'+OUTN+'.dat','w')
FOUT2.write(str(float(c4)/float(len(x)))+'\n')
FOUT2.write(str(float(c6)/float(len(x)))+'\n')
FOUT2.write(str(float(cdef)/float(len(x))))
FOUT2.close()

#print(part1_len_4,'part1_len4') # LN

f15.close()
f19.close()
f16.close()
f20.close()
f17.close()
f18.close()
#image=Image.open('voronoifinal.png')

#imageBox = image.getbbox()
#cropped=image.crop(imageBox)
#cropped.save('L_2d_cropped.png')
#print(float(c4)/float(len(x)), 'square fraction') # LN
#print(float(c6)/float(len(x)), 'hexagonal fraction') # LN
#print(float(cdef)/float(len(x)), 'defect fraction') # LN
#vect=1.0
#vect1x= x2[i]-x[k]
#vect1y= x3[i]-y[k]
#denom1= np.sqrt(1.0)
#vect2x= x2[i+1]-x[k]
#vect2y= x3[i+1]-y[k]
#print(dist_def, 'dist_def')


f27=open("defect_dist_"+OUTN+".dat",'w')

for i in range(len(dist_def)):
		f27.write("%s\n"%((dist_def[i])))

f27.close()
f28=open("ISV_dist_"+OUTN+".dat",'w')

for i in range(len(dist_ISV)):
		f28.write("%s\n"%((dist_ISV[i])))

f28.close()
f29=open("SVPS_dist_"+OUTN+".dat",'w')

for i in range(len(dist_SVPS)):
		f29.write("%s\n"%((dist_SVPS[i])))

f29.close()
####################################### use the output from here for the final plot (like-like neighbors) #######################
f30=open("pore4neigh_"+OUTN+".dat",'w')

for i in range(len(part1n4)):
		f30.write("%s\t%s\n"%((part1n4[i],part2n4[i])))

f30.close()
f31=open("pore4dist_"+OUTN+".dat",'w')
f50=open("pore4neigh_final_"+OUTN+".dat",'w')


dist_4=[]
for m in range(len(part1n4)):
    c=0
    for k in range(len(part1n4)):
        
        if(part2n4[m]==part1n4[k]):
            c=1#I found a similar neighbor pore
    if(c==1):
        DX= (x[part1n4[m]-1]-x[part2n4[m]-1])
        DY= (y[part1n4[m]-1]-y[part2n4[m]-1])
        #DX= DX - L*round(float(DX/L))
        #DY= DY - L*round(float(DY/L))
        if(np.sqrt((DX)**2.0+(DY)**2.0)<70):
            dist_4.append(np.sqrt((DX)**2.0+(DY)**2.0))
            f31.write("%s\n"%((np.sqrt((DX)**2.0+(DY)**2.0))))
            f50.write("%s\t%s\n"%(part1n4[m], part2n4[m]))

f31.close()
f50.close()



plt.figure()
count, bins, ignored = plt.hist(dist_4, 15, edgecolor='blue', density=True, linewidth=1.2)
plt.tick_params(labelsize=20)
#plt.xlim(0,70) #LN 321
plt.savefig('ISV_dist4_'+OUTN+'.png',bboxsize='tight') #LN310

#plt.show()  #LN310

xplot=[]
yplot=[]
for m in range(len(part1n4)):
    c=0
    for k in range(len(part1n4)):
        
        if(part2n4[m]==part1n4[k]):
            c=1#I found a similar neighbor pore
    if(c==1):
        DX= (x[part1n4[m]-1]-x[part2n4[m]-1])
        DY= (y[part1n4[m]-1]-y[part2n4[m]-1])
        #DX= DX - L*round(float(DX/L))
        #DY= DY - L*round(float(DY/L))
        if(np.sqrt((DX)**2.0+(DY)**2.0)<35):
            xplot.append(x[part1n4[m]-1])
            yplot.append(y[part1n4[m]-1])
            xplot.append(x[part2n4[m]-1])
            yplot.append(y[part2n4[m]-1])
            #dist_4.append(np.sqrt((DX)**2.0+(DY)**2.0))
            #f31.write("%s\n"%((np.sqrt((DX)**2.0+(DY)**2.0))))
            #f50.write("%s\t%s\n"%(part1n4[m], part2n4[m]))

figure = plt.figure()
figure.set_size_inches(20, 20)
ax = plt.gca()
ax.scatter(xplot, yplot, s=10, c="g")
plt.plot(HNx,HNy,'-k', linewidth= 5.0)
#plt.savefig('Modified.png',bboxsize='tight')
plt.axis('off')
plt.savefig('pointcheck_'+OUTN+'.png', dpi=600, bboxsize='tight') #LN310
#plt.show() #LN310


f32=open("pore6neigh_"+OUTN+".dat",'w') #LN310

for i in range(len(part1n6)):
		f32.write("%s\t%s\n"%((part1n6[i],part2n6[i])))

f33=open("pore6dist_"+OUTN+".dat",'w') #LN310
f43=open("pore6neigh_final_"+OUTN+".dat",'w') #LN310

#for i in range(len(part1n4)):
#		f31.write("%s\t%s\n"%((part1n4[i],part2n4[i])))
dist_6=[]
for m in range(len(part1n6)):
    c=0
    for k in range(len(part1n6)):
        
        if(part2n6[m]==part1n6[k]):
            c=1#I found a similar neighbor pore
    if(c==1):
        DX= (x[part1n6[m]-1]-x[part2n6[m]-1])
        DY= (y[part1n6[m]-1]-y[part2n6[m]-1])
        #DX= DX - L*round(float(DX/L))
        #DY= DY - L*round(float(DY/L))
        if(np.sqrt((DX)**2.0+(DY)**2.0)<70):
            dist_6.append(np.sqrt((DX)**2.0+(DY)**2.0))
            f33.write("%s\n"%((np.sqrt((DX)**2.0+(DY)**2.0))))
            f43.write("%s\t%s\n"%((part1n6[m], part2n6[m])))

plt.figure()
count, bins, ignored = plt.hist(dist_6, 15, edgecolor='blue', density=True, linewidth=1.2)
plt.tick_params(labelsize=20)
plt.savefig('SVPS_dist6_'+OUTN+'.png',bboxsize='tight') #LN310
#plt.show()  #LN310

f34=open("poredefneigh_"+OUTN+".dat",'w') #LN310

for i in range(len(part1ndef)):
		f34.write("%s\t%s\n"%((part1ndef[i],part2ndef[i])))

f35=open("poredefdist_"+OUTN+".dat",'w') #LN310

#for i in range(len(part1n4)):
#		f31.write("%s\t%s\n"%((part1n4[i],part2n4[i])))
dist_def_1=[]
for m in range(len(part1ndef)):
    c=0
    for k in range(len(part1ndef)):
        
        if(part2ndef[m]==part1ndef[k]):
            c=1#I found a similar neighbor pore
    if(c==1):
        DX= (x[part1ndef[m]-1]-x[part2ndef[m]-1])
        DY= (y[part1ndef[m]-1]-y[part2ndef[m]-1])
        #DX= DX - L*round(float(DX/L))
        #DY= DY - L*round(float(DY/L))
        if(np.sqrt((DX)**2.0+(DY)**2.0)<70):
            dist_def_1.append(np.sqrt((DX)**2.0+(DY)**2.0))
            f35.write("%s\n"%((np.sqrt((DX)**2.0+(DY)**2.0))))

plt.figure()
count, bins, ignored = plt.hist(dist_def_1, 15, edgecolor='blue', density=True, linewidth=1.2)
plt.tick_params(labelsize=20)
plt.savefig('def_'+OUTN+'.png',bboxsize='tight') #LN310
#plt.show() #LN310


##########################################################################################################################




#for m in range(len(part1)):
#            DX= (x[part1[m]-1]-x[part2[m]-1])
#            DY= (y[part1[m]-1]-y[part2[m]-1])
#            DX= DX - L*round(float(DX/L))
#            DY= DY - L*round(float(DY/L))
#            dist_SVPS.append(np.sqrt((DX)**2.0+(DY)**2.0))
plt.figure()
count, bins, ignored = plt.hist(dist_def, 15, edgecolor='blue', density=True, linewidth=1.2)
plt.tick_params(labelsize=20)
plt.savefig('defect_dist_'+OUTN+'.png',bboxsize='tight') #LN310
#plt.show()
#print(len(count))
#print(len(bins))
#print(ignored)


plt.figure()
count, bins, ignored = plt.hist(dist_ISV,15, edgecolor='blue', density=True, linewidth=1.2)
plt.tick_params(labelsize=20)
plt.savefig('ISV_dist_'+OUTN+'.png',bboxsize='tight') #LN310
#plt.show()   #LN310


plt.figure()
count, bins, ignored = plt.hist(dist_SVPS, 15, edgecolor='blue', density=True, linewidth=1.2)
plt.tick_params(labelsize=20)
plt.savefig('SVPS_dist_'+OUTN+'.png',bboxsize='tight') #LN310
#plt.show()   #LN310

plt.figure()
count, bins, ignored = plt.hist(lenrat, 15, edgecolor='black', linewidth=1.2)
#plt.show() #LN310
#print(lenrat, len(lenrat), 'lenrat')

plt.figure()
count, bins, ignored = plt.hist(angle, 15, edgecolor='black', linewidth=1.2)
plt.tick_params(labelsize=20)
plt.savefig('Angle_'+OUTN+'.png',bboxsize='tight') #LN310
#plt.show() #LN310
#print(angle, len(angle), 'angle')
#for i in len(neighnum):

f3=open("flag_"+OUTN+".dat",'w') #LN310
c= 0
for i in range(len(neighnew_proj)):
    c= c+1
    if(int(neighnew_proj[i])==int(4)):
        f3.write("%s\n"%(1))
    elif(int(neighnew_proj[i])==int(6)):
        f3.write("%s\n"%(3))
    else:
        f3.write("%s\n"%(2))

#print(len(neighnew_proj), 'length of angle neigh') #LN310
#print(c, 'length of angle neigh') #LN310

f4= open("neighnew_proj_"+OUTN+".dat",'w') #LN310

for i in range(len(neighnew_proj)):
	f4.write("%s\n"%(neighnew_proj[i]))

bins_list = [2, 3, 4, 5, 6, 7, 8, 9, 10]
plt.figure()
plt.bar(*np.unique(neighnew_proj, return_counts=True))
plt.title("Angle cutoff")
plt.tick_params(labelsize=20)
plt.savefig('neighnew_proj_'+OUTN+'.png',bboxsize='tight') #LN310
#plt.show()  #LN310



bins_list = [2, 3, 4, 5, 6, 7, 8, 9, 10]
plt.figure()
plt.bar(*np.unique(neighnew_len, return_counts=True))
plt.title("ratio cutoff")
plt.tick_params(labelsize=20)
plt.savefig('neighnew_len_'+OUTN+'.png',bboxsize='tight') #LN310
#plt.show()     #LN310
