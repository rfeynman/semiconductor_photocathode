import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt
from scipy.interpolate import interp1d
from scipy import integrate
from scipy.stats import expon
import random
import math

me=9.11*10**-31  # electron mess,kg; maybe need to transfer to effective mass
bandgap=1.64

def DosDist(n,loct,e_photon,thick,data):# n, the photon numbers; loct, the postion(z) of the electrons;thick, unit, nm.
    f = interp1d(data[:,0], data[:,1])

    n1 = int((e_photon-1)/0.01) # change to n 
    energy = np.linspace(1.,e_photon,n1)
    norm,err= integrate.quad(lambda e: f(e-e_photon)*f(e), 1,e_photon,limit=100)
    
    data_ene = []
    num_energy = []
    i=0
    while i< n1:
        n3 =round(n*f(energy[i]-e_photon)*f(energy[i])*0.01/norm) #using n instead of n1
        num_energy.append(n3)
        ener_array = np.empty(n3)
        ener_array.fill(energy[i])
        data_ene.extend(ener_array)
        i+=1
    np.random.shuffle(data_ene)
    print len(data_ene)
    
    plt.subplot(211)
    plt.plot(data[:,0],data[:,1])
    plt.subplot(212)
    plt.hist(data_ene,bins=40)
    plt.show()
    
    p2D=[]
    wl=((19.82-27.95*e_photon+11.15*e_photon**2)*10**-3 )**-1
    pens = expon.rvs(loc=0, scale=wl,size=n)
    penss=filter(lambda x:x<=thick,pens)
    params_exp=expon.fit(pens,floc=0)    
    
    i=0
    for i in range(len(penss)):
        phi=random.uniform(0,2*math.pi) # initial angular
        poy=random.uniform(-1*10**6,1*10**6)  # initial y direction position, nm
        v = np.sqrt(2 * (data_ene[i]-bandgap)*1.6*(10**(-19))/me)
        p2D.append([penss[i],poy,v*math.cos(phi),v*math.sin(phi),v,data_ene[i]-bandgap])  #p2D: (z,y,vz,vy,v,ene)
        i+=1  
    p2D=np.array(p2D)       
    return p2D,penss,params_exp,data_ene

def main():
    data = genfromtxt('K2CsSb DOS.csv',delimiter=',')
    p2D,penss,parames_exp,data_ene=DosDist(10000,0,4,50,data)

    print p2D
    me=9.11*10**-31

   

if __name__=='__main__':
    main()    