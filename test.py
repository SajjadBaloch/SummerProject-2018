import matplotlib.pylab as plt
File='Plot_Data.dat'
f=open(File,'r')
line1=f.readline()
lines=f.readlines()
f.close()
paras=line1.split()
N=int(paras[0])
L=float(paras[1])/2.
lim=1.1*L
dp=len(lines)
x=[[0 for i in xrange(dp)] for j in xrange(N)]
for j,line in zip(xrange(dp),lines):
	data=line.split()
	for i in xrange(N):
		x[i][j]=float(data[3*i])
fig=plt.figure(figsize=(24,12))
ax=plt.axes()
for i in xrange(N):
	ax.scatter(i+1,x[i][0],color='b',marker='.')
ax.set_xlabel(r'N')
ax.set_ylabel(r'x(t=0) [Mpc]')
ax.set_ylim(-lim,lim)
ax.set_xticks(range(1,N+1,1))
plt.grid()
plt.show()
