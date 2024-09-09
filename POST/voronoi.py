###########################UPDATED 3/3/2022 ######################################
############# Prajwal Bangalore Prakash ##################
# Minor chages Luis Nieves March 4th, 2022
import freud
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys as s  # LN

########### just the points ############################

fname = s.argv[1] # LN

fpath1 = fname
data1 = []
with open(fpath1) as infile:
	for line in infile:
		data1.append([ float(x) for x in line.split() ])

data1 = np.array(data1)
L= float(data1[0][0])
print(L)
#print(data1[3:], 'data1')
points = (data1[3:]-(L/2.0))
#points = data1[3:]
print(max(points[:,0]), 'max x coordinate')
print(max(points[:,1]), 'max y coordinate')
points[:,2]= 0.0
#print(points, 'points')
f=open("coordvoron.inp",'w')
for i in range(len(points)):
	f.write("%04.18e\t%04.18e\t%04.18e\n"%(points[i,0], points[i,1], points[i,2]))
#points= np.array([[-0.5, -0.5, 0], [0.5, -0.5, 0], [-0.5, 0.5, 0], [0.5, 0.5, 0]])

plt.scatter(points[:, 0], points[:, 1])
plt.title("Points")
plt.gca().set_aspect("equal")
#plt.savefig('points.png',bboxsize='tight') 
#plt.show()
####################### Box Length #######################
#L = 2
box = freud.box.Box.square(L)
voro = freud.locality.Voronoi()
cells = voro.compute((box, points)).polytopes
print(len(cells), 'length of cells')
print(cells[0][0], 'cells')
print(cells[0][1], 'cells')
print(cells[0][2], 'cells')
print(len(cells[0]), 'len cells')

print(cells[3716][0][0], 'cells 3716')

vertices=[]

cells_np = np.vstack((cells))
print(np.shape(cells_np))
length = []
#for i in range(0,1):
for i in range(len(cells)):
    length.append(np.repeat(i+1,np.shape(cells[i])[0]))
length = np.hstack((length))[:,np.newaxis]
all_arr = np.hstack((length,cells_np[:,:2]))
np.savetxt("neighborlistpoints.dat",all_arr,fmt = '%i %f %f')


for i in range(len(cells)):

    for j in range(len(cells[i])):
        vertices.append(cells[i][j])

print(vertices[0][0], 'vertices')
print(vertices[0][1], 'vertices')
print(vertices[0][2], 'vertices')
print(vertices[0], 'vertices')
print(vertices[1], 'vertices')
print(vertices[2], 'vertices')
print(len(vertices), 'len of vertices')
plt.figure()
ax = plt.gca()
ax.scatter(points[:, 0], points[:, 1], s=12, c="k")
#voro.plot(ax=ax, cmap= "RdBu")
voro.plot(ax=ax, cmap= "jet")
ax.scatter(points[:, 0], points[:, 1], s=12, c="k")
#plt.xlim(-2500, 2500)
#plt.show() # LN

#figure = plt.figure()
#figure.set_size_inches(20, 20)
#ax = plt.gca()
#ax.scatter(points[3659,0],points[3659,1], s=10, c="k")
#ax.scatter(points[1363,0],points[1363,1],s=10,c="b")
#ax.scatter(points[1571,0],points[1571,1],s=10,c="b")
#ax.scatter(points[2752,0],points[2752,1],s=10,c="b")
#ax.scatter(points[2773,0],points[2773,1],s=10,c="b")
#ax.scatter(points[2845,0],points[2845,1],s=10,c="b")
#ax.scatter(points[3907,0],points[3907,1],s=10,c="b")
#ax.fill(cells[3659][:,0], cells[3659][:,1], facecolor='red', edgecolor='black', linewidth=3.0, alpha= 0.4)
#plt.show()

cell_len=[]

for i in range(len(cells)):
	cell_len.append(float(len(cells[i])))
#print(cell_len)
bins_list = [2, 3, 4, 5, 6, 7, 8, 9, 10]
plt.figure()
plt.bar(*np.unique(cell_len, return_counts=True))
plt.title("Number of sides")
plt.tick_params(labelsize=20)
plt.savefig('number of sides.png',bboxsize='tight')
#plt.show()
#ax.scatter(points[:, 0], points[:, 1], s=12, c="k")
#plt.savefig('voronoi.png',bboxsize='tight') 
#plt.show()

#plt.figure()
#plt.hist(voro.volumes)
#plt.title("Voronoi cell volumes")
#plt.show()
#print(voro.volumes, 'voro volume')
nlist = voro.nlist
#m= property nlist
#nlist1= freud.locality.NeighborQueryResult.toNeighborList()
######## neighbor lists
print(nlist[1,:], 'nlist')
print(len(nlist[:]), 'len nilist')
nlist = np.asarray(nlist) 
nlist += 1
unique, counts = np.unique(nlist[:,0], return_counts=True)
count_array = np.asarray((unique, counts)).T
for i in range(len(cells)):
    if(len(cells[i])!=count_array[i,1]):
        temp = len(cells[i])
        cells[i] = np.unique(cells[i].round(decimals=4),axis=0)
        if(len(cells[i])<temp):
            print("Point removed because of overlap")
        else:
            print("Something is wrong. It will not work")
cells_np = np.vstack((cells))
print(np.shape(cells_np))
length = []
#for i in range(0,1):
for i in range(len(cells)):
    length.append(np.repeat(i+1,np.shape(cells[i])[0]))
length = np.hstack((length))[:,np.newaxis]
all_arr = np.hstack((length,cells_np[:,:2]))
np.savetxt("neighborlistpoints.dat",all_arr,fmt = '%i %f %f')

weights = np.asarray(voro.nlist.weights)[:,np.newaxis]
nlist = np.hstack((nlist,weights))        
np.savetxt("neighborlist.dat",nlist,fmt="%i %i %f")
print(np.shape(nlist),np.shape(cells_np),np.shape(voro.nlist.weights), 'nlist shape')

#print((nlist[:]),'lenght of nlist' )

# f=open("neighborlist.dat",'w')
# for i in range(0,len(nlist[:])):
	# part1=nlist[i,0]
	# part2=nlist[i,1]
	# dist= np.sqrt((points[part1,0]-points[part2,0])**2.0+(points[part1,1]-points[part2,1])**2.0)
	# dist=dist-(L*np.round(dist/L))
	# f.write("%s\t%s\t%04.18e\t%04.18e\n"%(int(nlist[i,0]+1.0),int(nlist[i,1]+1.0),voro.nlist.weights[i], abs(dist)/voro.nlist.weights[i]))
	# #f.write(str(nlist[i,0]+1.0),str(nlist[i,1]+1.0),float(voro.nlist.weights[i]))

# for i in range(0,len(nlist[:])):
    # if(int((nlist[i,0]+1.0))==3717):
        # print(vertices[i][0], vertices[i][1], 'vertex check')

   
f=open("neighnum.dat",'w')
for i in range(len(cells)):
#	for k in range(len(cells[i])):
	f.write("%s\n"%(int(len(cells[i]))))

	
# ###### corresponding neighbor weights ########
# print(voro.nlist.weights[0], 'weights')
# print(len(voro.nlist.weights), 'len weights')



# line_data = np.asarray(
    # [[points[i], points[i] + box.wrap(points[j] - points[i])] for i, j in nlist]
# )[:, :, :2]
# line_collection = matplotlib.collections.LineCollection(line_data, alpha=0.2)
#print(line_data, 'line collection')
#plt.figure()
#ax = plt.gca()
#voro.plot(ax=ax, cmap="jet")
#ax.add_collection(line_collection)
#ax.scatter(points[:, 0], points[:, 1], s=12, c="k")
#plt.show()










##############hexagonal lattice ###############
############## generating coordinates for hexagonal lattice points ############
#def hexagonal_lattice(rows=3, cols=3, noise=0, seed=None):
#    if seed is not None:
#        np.random.seed(seed)
#    # Assemble a hexagonal lattice
#    points = []
#    for row in range(rows * 2):
#        for col in range(cols):
#            x = (col + (0.5 * (row % 2))) * np.sqrt(3)
#            y = row * 0.5
#            points.append((x, y, 0))
#    points = np.asarray(points)
#    points += np.random.multivariate_normal(
#        mean=np.zeros(3), cov=np.eye(3) * noise, size=points.shape[0]
#    )
#    # Set z=0 again for all points after adding Gaussian noise
#    points[:, 2] = 0

    # Wrap the points into the box
#    box = freud.box.Box(Lx=cols * np.sqrt(3), Ly=rows, is2D=True)
#    points = box.wrap(points)
#    return box, points

# Compute the Voronoi diagram and plot
#box, points = hexagonal_lattice()
#voro = freud.locality.Voronoi()
#voro.compute((box, points))
#voro

# Compute the Voronoi diagram
#box, points = hexagonal_lattice(rows=12, cols=8, noise=0.03, seed=2)
#voro = freud.locality.Voronoi()
#voro.compute((box, points))

# Plot Voronoi with points and a custom cmap
#plt.figure()
#ax = plt.gca()
#voro.plot(ax=ax, cmap="RdBu")
#ax.scatter(points[:, 0], points[:, 1], s=2, c="k")
#plt.show()


############# volume of the Voronoi cells ##############

#plt.hist(voro.volumes)
#plt.title("Voronoi cell volumes")
#plt.show()

###########The Voronoi class also computes a freud.locality.NeighborList, where particles are neighbors if they share an edge in the Voronoi diagram. 
#The NeighborList effectively represents the bonds in the Delaunay triangulation. 
#The neighbors are weighted by the length (in 2D) or area (in 3D) between them. 
#The neighbor weights are stored in voro.nlist.weights


#nlist = voro.nlist
#line_data = np.asarray(
#    [[points[i], points[i] + box.wrap(points[j] - points[i])] for i, j in nlist]
#)[:, :, :2]
#line_collection = matplotlib.collections.LineCollection(line_data, alpha=0.2)
#plt.figure()
#ax = plt.gca()
#voro.plot(ax=ax, cmap="RdBu")
#ax.add_collection(line_collection)
#plt.show()

#print(voro.nlist.weights)

#files=['Grout80.dat']
#flip=[]
#x=[]
#y=[]
#c= 0
#for g in range(len(files)):
#	f1=open(files[g])
# 
#	for line in f1:
#		c= c+1
#		x.append(float(line.split()[0]))
#		y.append(float(line.split()[1]))

#dia=[]

#for i in range(len(s)):
#	dia.append(2.0*np.sqrt(s[i]/np.pi))

#f1= open("dia.dat", 'w')

#for i in range(len(dia)):
#print(dia)
#	f1.write("%04.18e\n"%(dia[i]))
#msq=0
#m=0
#for i in range(len(dia)):
#	msq=msq+(dia[i])**2.0
#	m=m+dia[i]
#msq=msq/len(dia)
#m=m/len(dia)

#poly=(msq/(m)**2.0)-1
#print(np.sqrt(poly), 'polydispersity')


#a1=np.arange(1.3,1.85,0.005)
#mu = 1.556
#sigma= 0.089
#print(sigma/mu)
#count, bins, ignored = plt.hist(dia, 15, edgecolor='black', linewidth=1.2)
#plt.plot(a1, 1/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (a1 - mu)**2 / (2 * sigma**2) ),linewidth=2, color='r')