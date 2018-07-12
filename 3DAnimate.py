import matplotlib.pylab as plt
import numpy as np
import sys
from matplotlib.pyplot import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from itertools import product, combinations

rep=False	# Animation repeats?
# Files to open
File  = '27Jun_fullneg'
Dat   = File+'_Data.dat'				# The Plot Data
Signs = File+'_Signs.dat'			# Signs of the Masses

# Open the output file 
f = open(Dat,'r')
# read all lines
line1 = f.readline()		# Parameters
lines = f.readlines()	# Data
#Close the file
f.close()

paras = line1.split()
# Number of bodies to plot
N=int(paras[0])
# Distance along x-axis to edge of box from origin
L=float(paras[1])/2.
# Plot range
lim=1.1*L
#Number of data points
dp = len(lines)

# Define the colors to use in the plot, using the signs of the masses
f = open(Signs,'r')
read = f.readline()
f.close
signs = [float(read.split()[i]) for i in xrange(N)]
colors=[0 for i in xrange(N)]
pos=cm.Blues(np.linspace(0,1,N))
neg=cm.Reds(np.linspace(0,1,N))
for i in xrange(N):
	if (signs[i]>0):
		colors[i] = 'b'
	else:
		colors[i] = 'r'

# Read data and make plots
# Read coordinates into 2D list, plot with different color for each body
# This should work automatically, with user only changing value of N
x=[[0 for j in xrange(dp)] for i in xrange(N)]
y=[[0 for j in xrange(dp)] for i in xrange(N)]
z=[[0 for j in xrange(dp)] for i in xrange(N)]
t=[0 for i in xrange(dp)]

j=0
for line in lines:
	data = line.split()
	for i in xrange(N):
		x[i][j]=float(data[3*i])
		y[i][j]=float(data[3*i+1])
		z[i][j]=float(data[3*i+2])
	j=j+1

# Setup the plot figure
fig=plt.figure(figsize=(12,12))
ax=plt.axes(xlim=(-lim,lim),ylim=(-lim,lim),zlim=(-lim,lim),projection='3d')
ax.view_init(elev=10,azim=60)

# Draw the Box
r = [-L,L]
for s,e in combinations(np.array(list(product(r, r, r))), 2):
	if np.sum(np.abs(s-e)) == r[1]-r[0]:
		ax.plot3D(*zip(s,e),color='b',linestyle='dotted')

plots = [ax.plot([],[],[],marker='.',color=colors[i])[0] for i in xrange(N)]
def init():
	for plot in plots:
		plot.set_data([],[])
		plot.set_3d_properties([])
	return plots
def animate(j):
	i=0
	k=1*j
	for plot in plots:
		plot.set_data([x[i][k]],[y[i][k]])
		plot.set_3d_properties([z[i][k]])
		i=i+1
	return plots	
anim = animation.FuncAnimation(fig,animate,init_func=init,frames=dp,interval=0,blit=True)

# Label the axes
ax.set_xlabel('x [Mpc]')
ax.set_ylabel('y [Mpc]')
ax.set_zlabel('z [Mpc]')

# Add a title to the plot
plt.title(str(N)+'-Body Simulation',size=18)

# Set Plot Legend
pos_mass = plt.Line2D([0],[0],linestyle="none",c='b',marker='.',label='+ve Mass')
neg_mass = plt.Line2D([0],[0],linestyle="none",c='r',marker='.',label='-ve Mass')
ax.legend(handles=[pos_mass,neg_mass],bbox_to_anchor=(0.,1.075),loc=2,borderaxespad=0.,numpoints=1)

# Show the plot
plt.show()
