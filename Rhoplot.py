import matplotlib.pylab as plt
import matplotlib.ticker as mtick
import sys
from matplotlib import rc

# Save the plot?
save=False
debug=False
# Files To Open
File  = '4cubAllPos'
Dat   = File+'_Data'
Image = Dat +'_DensityPlot'	# Name of the Image to Save

# Customisation variables
s=1.5			# Marker size
lw=1			# Line width
txt=16 		# Plot label text size
lblpad=11	# Plot label padding
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

# Read data and make plots
t=[0 for i in xrange(dp)]
deltap=[0 for i in xrange(dp)]
deltan=[0 for i in xrange(dp)]

print '=====Reading Data====='
for i,line in zip(xrange(dp),lines):
	data = line.split()
	t[i]=float(data[-4])/tdyn
	deltap[i]=float(data[-2])
	deltan[i]=float(data[-1])

print '======Data Read======='

if (debug):
	f = open('Debugging.txt','w')
	for i in xrange(dp):
		s = str([t[i],deltap[i],deltan[i]])+'\n'
		f.write(s)
	f.close()
	print 'Debug file written'
	sys.exit()	

print '===Generating Plots==='

tmin  = min(t)
tmax  = max(t)+0.01
dmin  = 0.5*10**(-1)
dmax  = 2*10**4


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

ax=plt.axes(xlim=(tmin,tmax))#,ylim=(dmin,dmax))
ax.plot(t,deltap,ls='solid',lw=1,color='b')
ax.plot(t,deltan,ls='solid',lw=1,color='r')

ax.hlines(0.,tmin,tmax,linestyles='dashed',lw=1.2)

#ax.set_xscale('log')
ax.set_yscale('log')

# Set axes labels
plt.xlabel('t/tdyn',size=txt)
plt.ylabel(r'$\delta$',size=1.3*txt)

# Set plot titles	
ax.set_title('Growth Rate - '+File,size=1.5*txt)

# Set Plot Legend
pos_mass = plt.scatter([],[],s=3*s,color='b',marker='.',linewidths=s,label='+ve Mass')
neg_mass = plt.scatter([],[],s=3*s,color='r',marker='.',linewidths=s,label='-ve Mass')
ax.legend(handles=[pos_mass,neg_mass],bbox_to_anchor=(0.,1.,1.,.075),loc=2,borderaxespad=0.,fontsize=0.8*txt,frameon=False)

print '===Plots Generated===='
# Save the figure
if (save):
	plt.savefig('/home/echo/cbm7/Summer_Project/Graphs/'+Image)
	print '=====Image Saved======'
else:
	# Show the Figure
	plt.show()

