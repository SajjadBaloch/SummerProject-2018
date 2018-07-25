import matplotlib.pylab as plt
import sys
from matplotlib import animation

interval=10	# Interval between frames in ms
sel=1			# To skip data points (to speed up animation)
rep=False		# Animation repeats?
Multi=False	# Multiple copies of the central box?
# Files to open
File  = '24-07-18_N50Test'
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
if (Multi):
	lim=3.1*L
else:
	lim=1.1*L
#Number of data points
dp = len(lines)

# Define the colors to use in the plot, using the signs of the masses
f = open(Signs,'r')
read = f.readline()
f.close
signs = [float(read.split()[i]) for i in xrange(N)]
colors=[0 for i in xrange(N)]
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

for j,line in zip(xrange(dp),lines):
	if (sel>1) and (j%sel==0):
		continue
	data = line.split()
	for i in xrange(N):
		x[i][j]=float(data[3*i])
		y[i][j]=float(data[3*i+1])

if (sel>1):
	for j in xrange(dp/sel):
		for i in xrange(N):
			x[i].remove(0.)
			y[i].remove(0.)

# Setup the plot figure
fig=plt.figure(figsize=(12,12))
ax=plt.axes(xlim=(-lim,lim),ylim=(-lim,lim))

plots = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
plots2 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
plots3 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
plots4 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
plots5 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
plots6 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
plots7 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
plots8 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
plots9 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]

def init():
	for plot in plots:
		plot.set_data([],[])
	return plots
# Center Box
def animate(j):
	for i,plot in zip(xrange(N),plots):
			plot.set_data([x[i][j]],[y[i][j]])
	return plots
# Lower Center Box
def animate2(j):
	for i,plot in zip(xrange(N),plots2):
			plot.set_data([x[i][j]],[y[i][j]-2*L])
	return plots2
# Lower Right Box
def animate3(j):
	for i,plot in zip(xrange(N),plots3):
			plot.set_data([x[i][j]+2*L],[y[i][j]-2*L])
	return plots3
# Lower Left Box
def animate4(j):
	for i,plot in zip(xrange(N),plots4):
			plot.set_data([x[i][j]-2*L],[y[i][j]-2*L])
	return plots4
# Center Right Box
def animate5(j):
	for i,plot in zip(xrange(N),plots5):
			plot.set_data([x[i][j]+2*L],[y[i][j]])
	return plots5
# Center Left Box
def animate6(j):
	for i,plot in zip(xrange(N),plots6):
			plot.set_data([x[i][j]-2*L],[y[i][j]])
	return plots6
# Upper Center Box
def animate7(j):
	for i,plot in zip(xrange(N),plots7):
			plot.set_data([x[i][j]],[y[i][j]+2*L])
	return plots7
# Upper Left Box
def animate8(j):
	for i,plot in zip(xrange(N),plots8):
			plot.set_data([x[i][j]-2*L],[y[i][j]+2*L])
	return plots8
# Upper Right Box
def animate9(j):
	for i,plot in zip(xrange(N),plots9):
			plot.set_data([x[i][j]+2*L],[y[i][j]+2*L])
	return plots9
	
# Center Box
anim  = animation.FuncAnimation(fig,animate,init_func=init,frames=dp,interval=interval,repeat=rep)
if (Multi):
	# Lower Center Box
	anim2 = animation.FuncAnimation(fig,animate2,init_func=init,frames=dp,interval=interval,repeat=rep)
	# Lower Right Box
	anim3 = animation.FuncAnimation(fig,animate3,init_func=init,frames=dp,interval=interval,repeat=rep)
	# Lower Left Box
	anim4 = animation.FuncAnimation(fig,animate4,init_func=init,frames=dp,interval=interval,repeat=rep)
	# Center Right Box
	anim5 = animation.FuncAnimation(fig,animate5,init_func=init,frames=dp,interval=interval,repeat=rep)
	# Center Left Box
	anim6 = animation.FuncAnimation(fig,animate6,init_func=init,frames=dp,interval=interval,repeat=rep)
	# Upper Center Box
	anim7 = animation.FuncAnimation(fig,animate7,init_func=init,frames=dp,interval=interval,repeat=rep)
	# Upper Left Box
	anim8 = animation.FuncAnimation(fig,animate8,init_func=init,frames=dp,interval=interval,repeat=rep)
	# Upper Right Box
	anim9 = animation.FuncAnimation(fig,animate9,init_func=init,frames=dp,interval=interval,repeat=rep)

# Label the axes
plt.xlabel('x [Mpc]')
plt.ylabel('y [Mpc]')

# Add a title to the plot
plt.title(str(N)+'-Body Simulation',size=18)

# Set Plot Legend
pos_mass = plt.scatter([],[],color='b',marker='.',label='+ve Mass')
neg_mass = plt.scatter([],[],color='r',marker='.',label='-ve Mass')
ax.legend(handles=[pos_mass,neg_mass],bbox_to_anchor=(0.,1.075),loc=2,borderaxespad=0.)

# Show the plot
plt.show()
