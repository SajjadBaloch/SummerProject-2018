import matplotlib.pylab as plt
import sys
from matplotlib import animation

interval=2	# Interval between frames in ms
rep=False	# Animation repeats?
Multi=True	# Multiple copies of the central box?
Save=False	# Save the animation?

# Files to open
File  = 'Plot'
Dat   = File+'_Data.dat'			# The Plot Data
Signs = File+'_Signs.dat'			# Signs of the Masses
Video = File+'_Animation.mp4'		# The video to save (only save at home)

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
# Dynamical timescale
tdyn=float(paras[2])
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
t=[0 for i in xrange(dp)]

for j,line in zip(xrange(dp),lines):
	data = line.split()
	for i in xrange(N):
		x[i][j]=float(data[3*i])
		y[i][j]=float(data[3*i+1])
	t[j]=float(data[-4])/tdyn

# Setup the plot figure
fig=plt.figure(figsize=(12,12))
ax=plt.axes(xlim=(-lim,lim),ylim=(-lim,lim))

plots  = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
if (Multi):
	plots2 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
	plots3 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
	plots4 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
	plots5 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
	plots6 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
	plots7 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
	plots8 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]
	plots9 = [ax.plot([],[],marker='.',color=colors[i])[0] for i in xrange(N)]

time = ax.text(0.8,1.02,'',transform=ax.transAxes,fontsize=14)

if (Multi):
	def init():
		for plot,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9 in	zip(plots,plots2,plots3,plots4,plots5,plots6,plots7,plots8,plots9):
			plot.set_data([],[])
			plot2.set_data([],[])
			plot3.set_data([],[])
			plot4.set_data([],[])
			plot5.set_data([],[])
			plot6.set_data([],[])
			plot7.set_data([],[])
			plot8.set_data([],[])
			plot9.set_data([],[])
		time.set_text('')
		return plots,plots2,plots3,plots4,plots5,plots6,plots7,plots8,plots9,time
	def animate(j):
		for i,plot,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9 in zip(xrange(N),plots,plots2,plots3,plots4,plots5,plots6,plots7,plots8,plots9):
			plot.set_data([x[i][j]],[y[i][j]])
			plot2.set_data([x[i][j]],[y[i][j]-2*L])
			plot3.set_data([x[i][j]+2*L],[y[i][j]-2*L])
			plot4.set_data([x[i][j]-2*L],[y[i][j]-2*L])
			plot5.set_data([x[i][j]+2*L],[y[i][j]])
			plot6.set_data([x[i][j]-2*L],[y[i][j]])
			plot7.set_data([x[i][j]],[y[i][j]+2*L])
			plot8.set_data([x[i][j]-2*L],[y[i][j]+2*L])
			plot9.set_data([x[i][j]+2*L],[y[i][j]+2*L])
		time.set_text("Time = "+str(t[j])[:5]+"tdyn")
		return plots,plots2,plots3,plots4,plots5,plots6,plots7,plots8,plots9,time
else:
	def init():
		for plot in plots:
			plot.set_data([],[])
		time.set_text('')
		return plots,time
	def animate(j):
		for i,plot in zip(xrange(N),plots):
				plot.set_data([x[i][j]],[y[i][j]])
		time.set_text("Time = "+str(t[j])[:5]+"tdyn")
		return plots,time
	
anim  = animation.FuncAnimation(fig,animate,init_func=init,frames=dp,interval=interval,repeat=rep)

# Label the axes
plt.xlabel('x [Mpc]')
plt.ylabel('y [Mpc]')

# Add a title to the plot
plt.title(str(N)+'-Body Simulation',size=18)

# Set Plot Legend
pos_mass = plt.scatter([],[],color='b',marker='.',label='+ve Mass')
neg_mass = plt.scatter([],[],color='r',marker='.',label='-ve Mass')
ax.legend(handles=[pos_mass,neg_mass],bbox_to_anchor=(0.,1.075),loc=2,borderaxespad=0.)

if (Save):
	# Save the animation
	anim.save('Animations/'+Video,fps=60)
	print '=====File '+Video+' Saved====='
else:
	# Show the plot
	plt.show()
