'''
Created on Feb 5, 2015

@author: wange

|----------------------------|
|                            |
|                            |
|                            |    ^y
|                            |    |
|----------------------------|    |
------------------------------>z
(z,y,vz,vy,v,ene)
{
The probability density function for expon is:

expon.pdf(x) = lambda * exp(- lambda*x)
for x >= 0.

The scale parameter is equal to scale = 1.0 / lambda.
rvs(loc=0, scale=1, size=1)
}
{
A special case of a chi distribution, with df = 3, loc = 0.0, and given scale = a, where a is the parameter used in the Mathworld description [R240].

The probability density function for maxwell is:

maxwell.pdf(x) = sqrt(2/pi)x**2 * exp(-x**2/2)
for x > 0.
rvs(loc=0, scale=1, size=1)    Random variates.
pdf(x, loc=0, scale=1)    Probability density function.
}

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


me=9.11*10**-31  # electron mess,kg
bandgap=1.64
'''
def Max_boltE(ene,le):  ## Maxwell b distribution
    return 2(((ene/pi)**0.5)/le**1.5)*np.exp(-ene/le)
'''


ec=1.6*(10**(-19))

def MBdist(n,e_photon,thick): # n: particle number, loct: start point(x-x0), scale: sigma, wl: wavelength,thick: thickness of the cathode
    assert e_photon> bandgap
    if e_photon-bandgap-0.8<=0:
        scale=e_photon-bandgap
        loct=0
    else:
        scale=0.8
        loct=e_photon-bandgap-scale
    data = maxwell.rvs(loc=loct, scale=scale, size=n)
    data_ene=np.array(data)
    params = maxwell.fit(data, floc=0)
    data_v=np.sqrt(2*data_ene*ec/me)
    p2D=[]
    wl=((19.82-27.95*e_photon+11.15*e_photon**2)*10**-3 )**-1
    pens = expon.rvs(loc=0, scale=wl,size=n)
    penss=filter(lambda x:x<=thick,pens)
    params_exp=expon.fit(pens,floc=0)
    i=0
    
    for n in range(len(penss)):
        phi=random.uniform(0,2*math.pi) # initial angular
        poy=random.uniform(-1*10**6,1*10**6)  # initial y direction position
        p2D.append([penss[i],poy,data_v[i]*math.cos(phi),data_v[i]*math.sin(phi),data_v[i],data[i]])  #p2D: (z,y,vz,vy,v,ene)
        i+=1  
    p2D=np.array(p2D)     
    return params,p2D,penss,params_exp




def DosDist(n,e_photon,thick,data):# n, the photon numbers; loct, the postion(z) of the electrons;thick, unit, nm.
    f = interp1d(data[:,0], data[:,1])

    n1 = int((e_photon-1)/0.01) # change to n 
    energy = np.linspace(1.,e_photon,n1)
    norm,err= integrate.quad(lambda e: f(e-e_photon)*f(e), 1,e_photon,limit=10000)
    
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
    
    '''
    plt.subplot(211)
    plt.plot(data[:,0],data[:,1])
    plt.subplot(212)
    plt.hist(data_ene,bins=40)
    plt.show()
    '''
    p2D=[]
    wl=((19.82-27.95*e_photon+11.15*e_photon**2)*10**-3 )**-1
    pens = expon.rvs(loc=0, scale=wl,size=n)
    penss=filter(lambda x:x<=thick,pens)
    params_exp=expon.fit(pens,floc=0)    
    
    i=0
    for i in range(len(penss)):
        phi=random.uniform(0,2*math.pi) # initial angular
        poy=random.uniform(-1*10**6,1*10**6)  # initial y direction position, nm
        v = np.sqrt(2 * np.abs((data_ene[i]-bandgap))*1.6*(10**(-19))/me)
        p2D.append([penss[i],poy,v*math.cos(phi),v*math.sin(phi),v,np.abs(data_ene[i]-bandgap)])  #p2D: (z,y,vz,vy,v,ene)
        i+=1  
    p2D=np.array(p2D) 
   
   
    return data_ene,p2D,penss,params_exp
    
def diff(stept,endT,p2Di,mfp_ps,bounde,bounda,types):
    
    
    p2D_emission=[]
    p2D_trap=[]
    p2D_back=[]

    
        
    if types==1:
        tmatix_ns=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[endT,0,1,0,0,0],[0,endT,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])# non scattering transistion matrix
        p2De=np.dot(p2Di,tmatix_ns)
        #p2De[:,0]=np.clip(p2De[:,0],bounde,bounda) 
        se,st,sb,p2De=swap(bounde,bounda,0,p2De)
        p2D_emission.extend(se.tolist())
        p2D_trap.extend(st.tolist())
        p2D_back.extend(sb.tolist())
        p2De=p2De.tolist()
        
        p2De.extend(p2D_back)
        p2De.extend(p2D_trap)
        p2De.extend(p2D_emission)
        
        p2De=np.array(p2De)
        p2De[:,5]=np.maximum(p2De[:,5],0)
        p2De[:,0]=np.clip(p2De[:,0],bounde,bounda) 
        
    elif types==2:
        tmatix_diff=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[endT,0,1,0,0,0],[0,endT,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
        t=0
        for t in range(int(endT/stept)):
            q=random.uniform(0,2*math.pi)
            tmatix_difft=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[stept,0,0,0,0,0],[0,stept,0,0,0,0],[0,0,math.cos(q),math.sin(q),1,0],[0,0,0,0,0,1]])
            tmatix_diff=np.dot(tmatix_diff,tmatix_difft)
            t+=1
        p2De=np.dot(p2Di,tmatix_diff)
        #p2De[:,0]=np.clip(p2De[:,0],bounde,bounda) 
        se,st,sb,p2De=swap(bounde,bounda,0,p2De)
        p2D_emission.extend(se.tolist())
        p2D_trap.extend(st.tolist())
        p2D_back.extend(sb.tolist())

        p2De=p2De.tolist()
        p2De.extend(p2D_back)
        p2De.extend(p2D_trap)
        p2De.extend(p2D_emission)
   
   
        p2De=np.array(p2De)
        print p2De
        p2De[:,5]=np.maximum(p2De[:,5],0)
        
        print p2De,bounde,bounda
        
        p2De[:,0]=np.clip(p2De[:,0],bounde,bounda) 
    elif types==3:
        tmatix_i=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[endT,0,1,0,0,0],[0,endT,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
        t=0
        p2De=np.dot(p2Di,tmatix_i)

        while t < endT:
            q=random.uniform(0,2*math.pi)
            ee=np.abs((0.01/5)*np.random.randn()+0.01 )  # random energy loss  2.5 * np.random.randn(2, 4) + 3// sigma * np.random.randn(...) + mu
            mfp=np.abs((mfp_ps/5)*np.random.randn()+mfp_ps)   #random mfp
            
            p2De[:,4]=np.sqrt(2*p2De[:,5]*ec/me)
            stept=mfp/((np.sqrt(2*np.mean(p2De[:,5])*ec)/me)*10**(-9))
            
            
            
            
            
            p2De[:,5]=p2De[:,5]-ee
            tmatix_epsim=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[stept,0,1,0,0,0],[0,stept,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
            p2De=np.dot(p2De,tmatix_epsim)
            se,st,sb,p2De=swap(bounde,bounda,0,p2De,ee)
            p2De[:,2]=p2De[:,4]*math.cos(q)
            p2De[:,3]=p2De[:,4]*math.sin(q)
            
            p2D_emission.extend(se.tolist())
            p2D_trap.extend(st.tolist())
            p2D_back.extend(sb.tolist())
           
            t+=stept
            
            if len(p2De)==0:
                #print p2De,p2D_trap
                break
            
        p2De=p2De.tolist()
        p2De.extend(p2D_back)
        p2De.extend(p2D_trap)
        p2De.extend(p2D_emission)
        p2De=np.array(p2De)
        p2De[:,5]=np.maximum(p2De[:,5],0)
        p2De[:,0]=np.clip(p2De[:,0],bounde,bounda) 
            
    
    else:
        print "wrong type: 1 is non scattering; 2 is brown diffusion 3 is simplified scattering"
            
   
    '''
    p2De[:,0]=map(lambda x:bounde if (x < bounde) else x,p2De[:,0])
    p2De[:,0]=map(lambda x:bounda if (x > bounda) else x,p2De[:,0])
    '''
    p2D_emission= np.array(p2D_emission)
    p2D_emission[:,0]=bounde
    p2D_back= np.array(p2D_back)
    p2D_back[:,0]=bounda
    p2D_trap = np.array(p2D_trap)
    if len(p2D_trap)!=0:
        p2D_trap[:,2:]=0
    return p2D_emission,p2D_trap,p2D_back,p2De

def swap(be,ba,bb,mat,el):
    assert ba > be
    emissind=mat[:,0]<=be
    backind=mat[:,0]>=ba
    loweind=mat[:,5]<=bb
    
    swape=mat[emissind,:]
    swape[:,5]=swape[:,5]+el
    swapt=mat[(~emissind) & (~backind)& loweind,:]
    swapb=mat[backind,:]
    swapr=mat[(~emissind) & (~loweind) & (~backind),:]
    #print swapt.shape
    return swape, swapt, swapb, swapr

def emission(mat,ea,sctk):
    
    matind=mat[:,5]>=(ea+sctk)
    mat2=mat[matind,:]
    
    mat2[:,5]=mat2[:,5]-ea-sctk
    mat2[:,4]=np.sqrt(2*mat2[:,5]*ec/me)
    mat2ind=np.abs(mat2[:,4])> np.abs(mat2[:,3])
    #print mat2ind
    emission=mat2[mat2ind,:]
    
    emission[:,2]=np.sqrt(emission[:,4]**2-emission[:,3]**2)
    surface_trap=len(mat)-len(emission)
    '''
    gamma=(np.mean(emission[:,5])+511000.0)/511000.0
    beta=np.sqrt(1-gamma**-2)
    Temittance= gamma*beta*np.sqrt(np.mean(emission[:,1]**2)*(10**-9)*np.mean((emission[:,3]/emission[:,2])**2)-np.mean(emission[:,1]*(10**-9)*(emission[:,3]/emission[:,2]))**2)  #sqrt(<x**2><x'**2>-<xx'>**2)
    '''
    Temittance=np.sqrt(np.mean((emission[:,1]*10**-9)**2))*np.sqrt(np.mean(emission[:,3]**2))/(3*10**8)  #sigma_x*sqrt(<v_x**2>)/c
    #print np.mean(emission[:,1]**2),np.mean(emission[:,1]**2),'\n',np.mean((emission[:,3]/emission[:,2])**2),'\n',np.mean(emission[:,1]*(emission[:,3]/emission[:,2]))**2
    #print np.sqrt(np.mean((emission[:,1]*10**-9)**2)),np.sqrt(np.mean(emission[:,3]**2))/(3*10**8)
    return emission, Temittance,surface_trap

def plot(p2D,pens,p2De,para,params_exp,thickness):
    
    plt.subplot(411)
    
    plt.hist(p2D[:,5], bins=10, normed=True)
    x = np.linspace(0, 5, 80)  
    
    plt.plot(x, maxwell.pdf(x, *para),'r',x,maxwell.cdf(x, *para), 'g')
    
    plt.subplot(412)
    plt.hist(pens, bins=20, normed=True)
    z=np.linspace(0,500,200)
    plt.xlim(0,2*thickness)  
    plt.plot(z, expon.pdf(z, *params_exp),'r')#plt.plot(z, expon.pdf(z, *params_exp),'r',z,expon.cdf(z, *params_exp), 'g')
    
    plt.subplot(413)
    plt.xlim(0,thickness+1)  
    plt.plot(p2D[:,0],p2D[:,1],'b.')
    
    plt.subplot(414)
    plt.ylim(-1*10**6,1*10**6)
    plt.xlim(0,thickness+1)
    plt.plot(p2De[:,0],p2De[:,1],'b.')
    plt.show()
    return

def scanplot(wlscan):
    
    fig,ax1=plt.subplots()
    ax2=ax1.twinx()
    
    ax1.plot(wlscan[:,0],wlscan[:,1]*10**6,'g-')
    ax2.plot(wlscan[:,0],wlscan[:,2],'b-')
    ax1.set_xlabel('photon energy[eV]')
    #plt.xlim(1.6,6.5)
    
    ax1.set_ylabel('emittance [mrad]',color='g')
   
    ax2.set_ylabel('QE[%]',color='b')
    ax2.set_yscale('log')
    ax2.margins(0.02)        
    plt.show()

def main(op):
    
    dosdata = genfromtxt('K2CsSb DOS.csv',delimiter=',')
    thickness=40
    n=10000
    pestart=3
    peend=6.2
    pestep=0.1
    petest=2.3
    electron_affinity=0.3
    schottky=0
    tstep=0.00001
    tend=0.01
    mfp=4
    emisposition=0
    
    if op == 0:
        para,p2D,pens,params_exp=MBdist(n,petest,thickness)# partcle #, start energy, sigma, apsorbtion length, thichness
       #p2D: initial energy distribution, pens: initial depth distribution, thickness: sample thickness, para and params_exp are fitting parameters
        esurf,trap,back,p2De=diff(tstep,tend,p2D,mfp,emisposition,thickness,3) # diff(stept,endT,p2Di,mfp_ps,bounde,bounda,types):
        # time step, time end, inital beam, mfp, emission surface, thickness, calculation type
        # esurf: reach to surface, trap: trap in body, back: get to back, p2De: all the particulates information
        emiss,emittance,surtrap=emission(esurf,electron_affinity,schottky)  # beam reach to surface, electron affinity, schottky
        #emiss: emit out beam, emittance:rms nor thermal emittance, surtrap: trap at surface
        #print emittance, float(len(emiss)/n)
        wlscan=[]
        for ep in np.arange(pestart,peend,pestep):
            emisst,emittancet, surtrapt=emission(diff(tstep,tend,MBdist(n,ep,thickness)[1],mfp,emisposition,thickness,3)[0],electron_affinity,schottky)
            qe=len(emisst)*100/float(n)
            #print ep, '    ',emittancet,'   ', qe,'%' 
            wlscan.append([ep, emittancet,qe])
        wlscan=np.array(wlscan)
        scanplot(wlscan)
        plot(p2D,pens,p2De,para,params_exp,thickness)



    elif op==1:
        data_ene,p2D,pens,params_exp=DosDist(n,petest,thickness,dosdata)
        #p2D: initial energy distribution, pens: initial depth distribution, thickness: sample thickness, para and params_exp are fitting parameters
        esurf,trap,back,p2De=diff(tstep,tend,p2D,mfp,emisposition,thickness,3) # diff(stept,endT,p2Di,mfp_ps,bounde,bounda,types):
        # time step, time end, inital beam, mfp, emission surface, thickness, calculation type
        # esurf: reach to surface, trap: trap in body, back: get to back, p2De: all the particulates information
        emiss,emittance,surtrap=emission(esurf,electron_affinity,schottky)  # beam reach to surface, electron affinity, schottky
        #emiss: emit out beam, emittance:rms nor thermal emittance, surtrap: trap at surface
        #print emittance, float(len(emiss)/n)
        print emiss
        
        wlscan=[]
        for ep in np.arange(pestart,peend,pestep):
            emisst,emittancet, surtrapt=emission(diff(tstep,tend,DosDist(n,ep,thickness,dosdata)[1],mfp,emisposition,thickness,3)[0],electron_affinity,schottky)
            qe=len(emisst)*100/float(n)
            #print ep, '    ',emittancet,'   ', qe,'%' 
            wlscan.append([ep, emittancet,qe])
            print ep
        wlscan=np.array(wlscan)
        scanplot(wlscan)
        plot(p2D,pens,p2De,None,params_exp,thickness)

       

if __name__ == '__main__':
    
    main(1)