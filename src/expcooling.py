'''
Created on Mar 17, 2015

@author: hxie
'''
from numpy import genfromtxt
import numpy
import matplotlib.pyplot as plt
data = genfromtxt('cooling_exp.csv',delimiter=',')
#data = numpy.array(data)
print data
x = data[:,0]

fig, ax1 = plt.subplots()
ax1.plot(data[:,0],data[:,1],'ro',label='%.1f' % data[0,1])
ax1.plot(data[:,0],data[:,2],'o',label='-1')
ax1.plot(data[:,0],data[:,3],'o',label='-32')
ax1.plot(data[:,0],data[:,4],'o',label='-50')
ax1.plot(data[:,0],data[:,5],'o',label='-72')
ax1.plot(data[:,0],data[:,6],'o',label='-80')
ax1.plot(data[:,0],data[:,7],'o',label='-89.7')
ax1.plot(data[:,0],data[:,8],'o',label='-93.5')
ax1.plot(data[:,0],data[:,9],'o',label='-96.8')
ax1.plot(data[:,0],data[:,10],'o',label='-100.4')
ax1.plot(data[:,0],data[:,11],'o',label='-107')
ax1.plot(data[:,0],data[:,12],'o',label='-107')
legend = ax1.legend(loc='best')
ax1.set_xlabel('photon energy(eV)')
ax1.set_ylabel('QE')
ax1.set_yscale('log')
plt.show()
