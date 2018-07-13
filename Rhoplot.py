import matplotlib.pylab as plt
import matplotlib.ticker as mtick
import numpy as np
import sys
from matplotlib.pyplot import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib import rc
from decimal import Decimal

# Save the plot?
save=False
# Files To Open
File  = 'Plot'
Dat   = File+'_Data'
Signs = File+'_Signs'			# Signs of the Masses
Image = Dat +'_DensityPlot'	# Name of the Image to Save

# Customisation variables
s=1			# Marker size
lw=1		# Line width
txt=16 		# Plot label text size
lblpad=10		# Plot label padding
tick=12		# Tick text size
tickh=12		# Tick size
pad=1			# Tick padding

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
t=[0 for i in xrange(dp)]
delta=[[0 for i in xrange(dp)] for j in xrange(N)]

print '=====Reading Data====='
for i,line in zip(xrange(dp),lines):
	data = line.split()
	t[i]=float(data[-(N+2)])/tdyn
	for j in xrange(N):
		delta[j][i]=float(data[-(N-j)])
print '======Data Read======='
print '===Generating Plots==='

tmin  = min(t)
tmax  = max(t)+0.01

# Customisation
plt.rc('path',simplify=True)
plt.rc('savefig',dpi=300,format='png')
plt.rc('axes',linewidth=lw,labelpad=lblpad)
plt.rc('xtick',labelsize=tick,direction='inout')
plt.rc('xtick.major',size=tickh,width=lw,pad=pad)
plt.rc('ytick',labelsize=tick,direction='inout')
plt.rc('ytick.major',size=tickh,width=lw,pad=pad)

# Set up the plot
plt.figure(figsize=(24,12))

ax=plt.axes(xlim=(tmin,tmax))
for i in xrange(N):
	ax.plot(t,delta[i],color=colors[i],marker='o',ms=s,lw=0.5,ls='none')
	ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
	ax.set_yscale('log')

# Set axes labels
plt.xlabel('t/tdyn',size=txt)
plt.ylabel(r'$\log\delta$',size=1.3*txt)

# Set plot titles	
ax.set_title('Local Densities',size=1.5*txt)

# Set Plot Legend
pos_mass = plt.scatter([],[],s=3*s,color='b',marker='.',linewidths=s,label='+ve Mass')
neg_mass = plt.scatter([],[],s=3*s,color='r',marker='.',linewidths=s,label='-ve Mass')
ax.legend(handles=[pos_mass,neg_mass],bbox_to_anchor=(0.,1.,1.,.075),loc=2,borderaxespad=0.,fontsize=0.8*txt,frameon=False)

print '====Plots Generated==='
# Save the figure
if (save):
	plt.savefig('/home/echo/cbm7/Summer_Project/Graphs/'+Image)
	print '======Image Saved====='
else:
	# Show the Figure
	plt.show()

