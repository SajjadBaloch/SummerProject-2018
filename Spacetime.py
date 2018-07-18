import matplotlib.pylab as plt
import sys
from matplotlib import rc

# Save the plot?
save=False
# Plot the mirror copies of the data?
Multi=True
# Filenames
File  = 'Plot'
Dat   = File+'_Data'				# The Plot Data
Signs = File+'_Signs'			# Signs of the Masses
Image = File+'_SpacetimePlot'	# Name of the Image to Save

# Customisation variables
s=0.5			# Marker size for spacial plots
ms=0.06		# Marker size for spacetime plots
ml=0.06		# Marker edgeline size for spacetime plots
lw=0.2		# Line widths
txt=6 		# Plot label text size
lblpad=1		# Plot label padding
tick=4		# Tick text size
tickh=4		# Tick size
pad=0.7		# Tick padding

# Open the output file 
f = open(Dat+'.dat','r')
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
	lim=3.2*L
else:
	lim=1.1*L
# Dynamical timescale
tdyn=float(paras[2])
#Number of data points
dp = len(lines)

# Define the colors to use in the plot, using the signs of the masses
f = open(Signs+'.dat','r')
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
x =[[0 for i in xrange(dp)] for j in xrange(N)]
xp=[[0 for i in xrange(dp)] for j in xrange(N)]
xm=[[0 for i in xrange(dp)] for j in xrange(N)]
y =[[0 for i in xrange(dp)] for j in xrange(N)]
yp=[[0 for i in xrange(dp)] for j in xrange(N)]
ym=[[0 for i in xrange(dp)] for j in xrange(N)]
z =[[0 for i in xrange(dp)] for j in xrange(N)]
zp=[[0 for i in xrange(dp)] for j in xrange(N)]
zm=[[0 for i in xrange(dp)] for j in xrange(N)]
t =[0 for i in xrange(dp)]

print '=====Reading Data====='
for j,line in zip(xrange(dp),lines):
	data = line.split()
	for i in xrange(N):
		x[i][j]=float(data[3*i])
		y[i][j]=float(data[3*i+1])
		z[i][j]=float(data[3*i+2])
		if (Multi):
			xp[i][j]=x[i][j]+2*L
			xm[i][j]=x[i][j]-2*L
			yp[i][j]=y[i][j]+2*L
			ym[i][j]=y[i][j]-2*L
			zp[i][j]=z[i][j]+2*L
			zm[i][j]=z[i][j]-2*L
		t[j]=float(data[-(N+2)])/tdyn
print '======Data Read======='
print 'ttot =',t[-1],'tdyn'
print '===Generating Plots==='

tlim = max(t)+0.1	
# Customisation
plt.rc('path',simplify=True)
plt.rc('savefig',dpi=300,format='png')
plt.rc('axes',linewidth=lw,labelpad=lblpad)
plt.rc('xtick',labelsize=tick,direction='inout')
plt.rc('xtick.major',size=tickh,width=lw,pad=pad)
plt.rc('ytick',labelsize=tick,direction='inout')
plt.rc('ytick.major',size=tickh,width=lw,pad=pad)

# Set up the plot
plt.figure(dpi=200)

nrows=3;ncols=4
# Set up plot figures
ax1=plt.subplot2grid((nrows,ncols),(0,0),rowspan=1,colspan=1)
ax2=plt.subplot2grid((nrows,ncols),(1,0),rowspan=1,colspan=1,sharex=ax1,sharey=ax1)
ax3=plt.subplot2grid((nrows,ncols),(2,0),rowspan=1,colspan=1,sharex=ax1,sharey=ax1)
ax4=plt.subplot2grid((nrows,ncols),(0,1),rowspan=1,colspan=2)
ax5=plt.subplot2grid((nrows,ncols),(1,1),rowspan=1,colspan=2,sharex=ax4,sharey=ax4)
ax6=plt.subplot2grid((nrows,ncols),(2,1),rowspan=1,colspan=2,sharex=ax4,sharey=ax4)
ax7=plt.subplot2grid((nrows,ncols),(0,3),rowspan=1,colspan=1,sharex=ax1,sharey=ax1)
ax8=plt.subplot2grid((nrows,ncols),(1,3),rowspan=1,colspan=1,sharex=ax1,sharey=ax1)
ax9=plt.subplot2grid((nrows,ncols),(2,3),rowspan=1,colspan=1,sharex=ax1,sharey=ax1)

# Plot the data
for i in xrange(N):
	if (i % 5 == 0):
		percent=float(i)/float(N)*100.
		print str(percent)+'%'
	
	# Initial Positions
	ax1.scatter(x[i][0],y[i][0],s=s,color=colors[i],marker='.',linewidths=s)
	ax2.scatter(x[i][0],z[i][0],s=s,color=colors[i],marker='.',linewidths=s)
	ax3.scatter(y[i][0],z[i][0],s=s,color=colors[i],marker='.',linewidths=s)

	# Spacetime Plots
	ax4.scatter(t,x[i],s=ms,color=colors[i],marker='.',linewidths=ml)
	ax5.scatter(t,y[i],s=ms,color=colors[i],marker='.',linewidths=ml)
	ax6.scatter(t,z[i],s=ms,color=colors[i],marker='.',linewidths=ml)

	# Final Positions
	ax7.scatter(x[i][-1],y[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
	ax8.scatter(x[i][-1],z[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
	ax9.scatter(y[i][-1],z[i][-1],s=s,color=colors[i],marker='.',linewidths=s)

	if (Multi):	
		ax1.scatter(x[i][0],yp[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax1.scatter(x[i][0],ym[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax1.scatter(xp[i][0],y[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax1.scatter(xm[i][0],y[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax1.scatter(xp[i][0],yp[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax1.scatter(xp[i][0],ym[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax1.scatter(xm[i][0],yp[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax1.scatter(xm[i][0],ym[i][0],s=s,color=colors[i],marker='.',linewidths=s)

		ax2.scatter(x[i][0],zp[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax2.scatter(x[i][0],zm[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax2.scatter(xp[i][0],z[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax2.scatter(xm[i][0],z[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax2.scatter(xp[i][0],zp[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax2.scatter(xp[i][0],zm[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax2.scatter(xm[i][0],zp[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax2.scatter(xm[i][0],zm[i][0],s=s,color=colors[i],marker='.',linewidths=s)

		ax3.scatter(y[i][0],zp[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax3.scatter(y[i][0],zm[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax3.scatter(yp[i][0],z[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax3.scatter(ym[i][0],z[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax3.scatter(yp[i][0],zp[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax3.scatter(yp[i][0],zm[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax3.scatter(ym[i][0],zp[i][0],s=s,color=colors[i],marker='.',linewidths=s)
		ax3.scatter(ym[i][0],zm[i][0],s=s,color=colors[i],marker='.',linewidths=s)

		ax4.scatter(t,xp[i],s=ms,color=colors[i],marker='.',linewidths=ml)
		ax4.scatter(t,xm[i],s=ms,color=colors[i],marker='.',linewidths=ml)

		ax5.scatter(t,yp[i],s=ms,color=colors[i],marker='.',linewidths=ml)
		ax5.scatter(t,ym[i],s=ms,color=colors[i],marker='.',linewidths=ml)

		ax6.scatter(t,zp[i],s=ms,color=colors[i],marker='.',linewidths=ml)
		ax6.scatter(t,zm[i],s=ms,color=colors[i],marker='.',linewidths=ml)

		ax7.scatter(x[i][-1],yp[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax7.scatter(x[i][-1],ym[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax7.scatter(xp[i][-1],y[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax7.scatter(xm[i][-1],y[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax7.scatter(xp[i][-1],yp[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax7.scatter(xp[i][-1],ym[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax7.scatter(xm[i][-1],yp[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax7.scatter(xm[i][-1],ym[i][-1],s=s,color=colors[i],marker='.',linewidths=s)

		ax8.scatter(x[i][-1],zp[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax8.scatter(x[i][-1],zm[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax8.scatter(xp[i][-1],z[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax8.scatter(xm[i][-1],z[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax8.scatter(xp[i][-1],zp[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax8.scatter(xp[i][-1],zm[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax8.scatter(xm[i][-1],zp[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax8.scatter(xm[i][-1],zm[i][-1],s=s,color=colors[i],marker='.',linewidths=s)

		ax9.scatter(y[i][-1],zp[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax9.scatter(y[i][-1],zm[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax9.scatter(yp[i][-1],z[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax9.scatter(ym[i][-1],z[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax9.scatter(yp[i][-1],zp[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax9.scatter(yp[i][-1],zm[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax9.scatter(ym[i][-1],zp[i][-1],s=s,color=colors[i],marker='.',linewidths=s)
		ax9.scatter(ym[i][-1],zm[i][-1],s=s,color=colors[i],marker='.',linewidths=s)

# Set axes labels and plot ranges
ax1.set_xlabel(r'x(t=0) [Mpc]',size=txt)
ax1.set_ylabel(r'y(t=0) [Mpc]',size=txt)
ax2.set_xlabel(r'x(t=0) [Mpc]',size=txt)
ax2.set_ylabel(r'z(t=0) [Mpc]',size=txt)
ax3.set_xlabel(r'y(t=0) [Mpc]',size=txt)
ax3.set_ylabel(r'z(t=0) [Mpc]',size=txt)

ax4.set_xlabel(r't/tdyn',size=txt)
ax4.set_ylabel(r'x(t) [Mpc]',size=txt)
ax5.set_xlabel(r't/tdyn',size=txt)
ax5.set_ylabel(r'y(t) [Mpc]',size=txt)
ax6.set_xlabel(r't/tdyn',size=txt)
ax6.set_ylabel(r'z(t) [Mpc]',size=txt)

ax7.set_xlabel(r'x(t=tfin) [Mpc]',size=txt)
ax7.set_ylabel(r'y(t=tfin) [Mpc]',size=txt)
ax8.set_xlabel(r'x(t=tfin) [Mpc]',size=txt)
ax8.set_ylabel(r'z(t=tfin) [Mpc]',size=txt)
ax9.set_xlabel(r'y(t=tfin) [Mpc]',size=txt)
ax9.set_ylabel(r'z(t=tfin) [Mpc]',size=txt)

ax1.set_xlim(-lim,lim)
ax1.set_ylim(-lim,lim)
ax4.set_xlim(-0.1,tlim)
ax4.set_ylim(-lim,lim)

# Set plot titles	
ax1.set_title('Initial Positions',size=1.2*txt)
ax4.set_title(str(N)+'-Body Simulation; Data from file: '+Dat+'.dat',size=1.5*txt)
#ax4.set_title('Force Softened - All +ve mass',size=1.5*txt)
ax7.set_title('Final Positions',size=1.2*txt)

# Set Plot Legend
pos_mass = plt.scatter([],[],s=s,color='b',marker='.',linewidths=s,label='+ve Mass')
neg_mass = plt.scatter([],[],s=s,color='r',marker='.',linewidths=s,label='-ve Mass')
ax1.legend(handles=[pos_mass,neg_mass],bbox_to_anchor=(-.28,1.03,1.,.12),loc=2,borderaxespad=0.,fontsize=0.7*txt,frameon=False)

# Format the figure margins
plt.subplots_adjust(left=0.05,right=0.98,top=0.95,bottom=0.05,wspace=0.4,hspace=0.2)

print '====Plots Generated==='
# Save the figure
if (save):
	plt.savefig('/home/echo/cbm7/Summer_Project/Graphs/'+Image)
	print '======Image Saved====='
else:
	# Show the Figure
	plt.show()
