import matplotlib.pylab as plt

File  = 'Plot'
Dat   = File+'_Data'				# The Plot Data
Signs = File+'_Signs'			# Signs of the Masses

f=open(Dat+'.dat','r')
line1=f.readline()
lines=f.readlines()
f.close()

paras=line1.split()
N=int(paras[0])
L=float(paras[1])/2.
lim=1.1*L
dp=len(lines)

x=[[0 for i in xrange(dp)] for j in xrange(N)]
y=[[0 for i in xrange(dp)] for j in xrange(N)]
z=[[0 for i in xrange(dp)] for j in xrange(N)]

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

for j,line in zip(xrange(dp),lines):
	data=line.split()
	for i in xrange(N):
		x[i][j]=float(data[3*i])
		y[i][j]=float(data[3*i+1])
		z[i][j]=float(data[3*i+2])

fig=plt.figure(figsize=(24,12))
nrows=3;ncols=1
ax1=plt.subplot2grid((nrows,ncols),(0,0),rowspan=1,colspan=1)
ax2=plt.subplot2grid((nrows,ncols),(1,0),rowspan=1,colspan=1,sharex=ax1,sharey=ax1)
ax3=plt.subplot2grid((nrows,ncols),(2,0),rowspan=1,colspan=1,sharex=ax1,sharey=ax1)

for i in xrange(N):
	ax1.scatter(i+1,x[i][0],color=colors[i],marker='.')
	ax2.scatter(i+1,y[i][0],color=colors[i],marker='.')
	ax3.scatter(i+1,z[i][0],color=colors[i],marker='.')

ax1.set_title(r'Initial Positions at t=0')
ax1.set_ylabel(r'x-Position')
ax2.set_ylabel(r'y-Position')
ax3.set_ylabel(r'z-Position')
ax3.set_xlabel(r'Particle')
ax1.set_xlim(0.9,N+0.1)
ax1.set_ylim(-lim,lim)
ax1.set_xticks(range(1,N+1,1))
ax1.grid()
ax2.grid()
ax3.grid()
plt.subplots_adjust(left=0.05,right=0.98,top=0.95,bottom=0.05,wspace=0.4,hspace=0.1)
plt.show()

