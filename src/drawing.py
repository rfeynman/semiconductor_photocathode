'''
Created on Mar 30, 2015

@author: wange
'''
import random
from scipy.stats import maxwell,expon
from scipy.constants import *
import matplotlib.pyplot as plt
import numpy as np

import math
from numpy import genfromtxt
from scipy.interpolate import interp1d
from scipy import integrate
import time

simdata= genfromtxt('Cool down beam.csv',delimiter=',')
expdata=genfromtxt('cooling_exp.csv',delimiter=',')
fig,ax1=plt.subplots(figsize=(16,8))

ax1.plot(simdata[:,0],simdata[:,2],'b-',label='sim_%.1f' % simdata[0,1])
ax1.plot(simdata[:,4],simdata[:,6],'r-',label='sim_%.1f' % simdata[0,5])
ax1.plot(simdata[:,8],simdata[:,10],'g-',label='sim_%.1f' % simdata[0,9])
ax1.plot(simdata[:,12],simdata[:,14],'y-',label='sim_%.1f' % simdata[0,13])
ax1.plot(simdata[:,16],simdata[:,18],'k-',label='sim_%.1f' % simdata[0,17])

ax1.plot(expdata[:,0],expdata[:,1],'bo',label='exp_%.1f' % expdata[0,1])
ax1.plot(expdata[:,0],expdata[:,3],'ro',label='exp_%.1f' % expdata[0,3])
ax1.plot(expdata[:,0],expdata[:,5],'go',label='exp_%.1f' % expdata[0,5])
ax1.plot(expdata[:,0],expdata[:,8],'yo',label='exp_%.1f' % expdata[0,8])
ax1.plot(expdata[:,0],expdata[:,11],'ko',label='exp_%.1f' % expdata[0,11])

legend = ax1.legend(loc='best')
ax1.set_xlabel('photon energy[eV]')
   
ticks=np.arange(1.5,6.0,0.2)
ax1.set_xticks(ticks)
         
ax1.set_ylabel('rate[%]',color='b')
ax1.set_yscale('log')
ax1.margins(0.02)  


plt.show()