'''
Created on Feb 12, 2015

@author: wange
'''
################################################
# 3D Lennard-Jones (Reduced Units) 
# NVE Molecular Dyamics Python:
# Based on R. LeSars Computational
# Materials Science textbook
# Code intended for learning purposes.

# This code is written in away that should
# make it easy to tranport to other languages
# especially C and FORTRAN as well as make it
# into seperate modules.

# Origional Author: Stefan Bringuier
# Affiliation: University of Arizona
# Department: Materials Science and Engineering
# Date: August 13, 2013

# If you find errors please report them to
# stefanb at email dot arizona dot edu

# No license or warrenty is provided.
# You are free to modify, add ,
# or redistribute in anyway.
##################################################

import numpy as np
from math import sqrt

# LJ 3D simulation
# Reduced Units

def LAMMPSdump(n,t,l,p,v,f):
    ''' Output dump file in lammps custom format
    n - number of atoms
    t - timestep
    l - box length (same in all directions)
    p - position array n x 3
    v - velocity array n x 3
    f - force array n x 3
    output:
    dump.LJ-3D-MD - LAMMPS dump file
    '''

    fil = open('dump.LJ-3D-MD','a')

    #Header 
    fil.write('ITEM: TIMESTEP \n')
    fil.write('%i \n' %(t))
    fil.write('ITEM: NUMBER OF ATOMS \n')
    fil.write('%i \n' %(n))
    fil.write('ITEM: BOX BOUNDS pp pp pp \n')
    fil.write('0.00 %f \n' %(l))
    fil.write('0.00 %f \n' %(l))
    fil.write('0.00 %f \n' %(l))
    fil.write('ITEM: ATOMS id type x y z vx vy vz fx fy fz \n')

    for i in xrange(n):
        ids = i + 1
        #Unscale coordinates
        x,y,z = p[i,0]*l,p[i,1]*l,p[i,2]*l
        vx,vy,vz = v[i,0],v[i,1],v[i,2]
        fx,fy,fz = f[i,0],f[i,1],f[i,2]
        fil.write('%i 1 %f %f %f %f %f %f %f %f %f \n' 
                  %(ids,x,y,z,vx,vy,vz,fx,fy,fz))

    fil.close()

def FCC(nc):
    ''' Faced-centerd cubic lattice
    input:
    nc - number of cells
    output:
    poss - scaled postions in n x dim array
    '''
    dim = 3 #3D
    nbasis = 4 # Number of atoms
    basis = np.array([[0.0,0.0,0.0],
                      [0.5,0.5,0.0],
                      [0.5,0.0,0.5],
                      [0.0,0.5,0.5]])

    natoms = (nc*nc*nc)*nbasis

    #Scaled coordinates (vector arrays)
    poss = np.zeros((natoms,dim))

    #TODO - NumPy-iz
    index = 0
    for i in xrange(nc):
        for j in xrange(nc):
            for k in xrange(nc):
                for l in xrange(nbasis):
                    poss[index][0] = (basis[l,0] + i)/nc
                    poss[index][1] = (basis[l,1] + j)/nc
                    poss[index][2] = (basis[l,2] + k)/nc
                    index += 1

    return natoms,poss

def MBdist(n,t):
    '''
    Maxwell-Boltzmann distribution 
    input:
    n - number atoms
    t - target temperature
    ouput:
    vel - velocity distribution array
    '''

    dim = 3 #dimensions

    vel = np.zeros((n,dim))
    momentum = np.zeros((dim,1))
    
    for d in xrange(dim):
        r = np.random.rand(n)
        veldist = np.sqrt(-2.0 * np.log(r)) * np.cos(2.0 * np.pi * r)        
        vel[:,d] = veldist
        momentum[d] = np.sum(vel[:,d])

    #Subtract net momentum/ Sum kinetic
    permom = momentum/n
    kengr = 0.00
    for d in xrange(dim):
        vel[:,d] -= permom[d]
        kengr += np.sum(np.square(vel[:,d]))
                
    #Rescale Kinetic energy to temp
    kengr *= 0.5 
    ktarget = 3.0/2.0 * ( n * t)
    rescale = np.sqrt(ktarget / kengr)

    vel *= rescale

    return vel

def LJ(rij):
    ''' Lennard-Jones Function and Derivative
    input:
    rij - seperation distance between i and j
    output:
    phi - energy
    dphi - force magnitude
    
    Note: function not used
    '''

    phi = 4 * ( (1/rij**12) - (1/rij**6))
    dphi = 24/rij**2 * (2 * (1/rij**12) - (1/rij**6))
    return phi,dphi

def ForceCalc(l,n,rc,p):
    ''' calculation routine driver
    input:
    l - simulation box length (cubic only)
    n - number of atoms
    rc - cutoff distance
    p - position array n x 3
    ouput:
    pote - potential energy accum
    vir - virial term
    forc - forc array

    Notes: This is fairly slow for one because of the double
    for loop in python and also no use of neighbor list.
    '''

    #TODO - speed up/ numpize
    dim = 3
    vshift = 1.0/(rc**12) - 1.0/(rc**6)
    forc = np.zeros((n,dim))
    
    virial = 0.00
    pote = 0.00

    for i in xrange(0,n-1):
        ftx = 0.00
        fty = 0.00
        ftz = 0.00
        for j in xrange(i+1,n):
            dx = p[j,0] - p[i,0]
            dy = p[j,1] - p[i,1]
            dz = p[j,2] - p[i,2]

            #Min. Image Conv. 
            #Recall positions are scaled
            dx = dx - round(dx)
            dy = dy - round(dy)
            dz = dz - round(dz)

            dist = l * sqrt((dx**2) + (dy**2) + (dz**2))
            if dist <= rc:
                phi =  1.0/(dist**12) - 1.0/(dist**6)
                dphi = 2.0/(dist**12) - 1.0/(dist**6)
                
                pote = pote + phi - vshift
                virial = virial + dphi

                ffx = (dphi * l * dx) / dist**2
                ffy = (dphi * l * dy) / dist**2
                ffz = (dphi * l * dz) / dist**2
                
                ftx = ftx + ffx
                fty = fty + ffy
                ftz = ftz + ffz
                
                #Newtons 3d law
                forc[j,0] = forc[j,0] - ffx
                forc[j,1] = forc[j,1] - ffy
                forc[j,2] = forc[j,2] - ffz
        
        #Sum forces on atom i
        forc[i,0] = forc[i,0] + ftx
        forc[i,1] = forc[i,1] + fty
        forc[i,2] = forc[i,2] + ftz


    #Add factor 4epsilon (normalized)
    pote = 4.0*pote 
    virial = 24.0*virial
    forc *= -24.0
    
    return pote,virial,forc

def PreIntegrate(n,dt,l,p,po,v,f):
    ''' Update positions after distribution
    input:
    n - number of aoms
    dt - timestep
    l - length of box
    p - positions n x 3 array
    po - old positions n x 3 array
    v - velocity n x 3 array
    f - forces n x 3 array
    '''
    dim = 3
    dtDl = dt/l
    hdt2l = 0.50 * (dt**2 / l)
    for i in xrange(n):
        po[i,0] = p[i,0] - v[i,0]*dtDl + hdt2l * f[i,0]
        po[i,1] = p[i,1] - v[i,1]*dtDl + hdt2l * f[i,1]
        po[i,2] = p[i,2] - v[i,2]*dtDl + hdt2l * f[i,2]


def Integrate(n,dt,l,p,po,v,f):
    ''' Verlet Integration of equations of motion
    input:
    n - number of atoms
    dt - timestep
    p - position array n x 3
    po - old positions
    v - velocity array n x 3
    f - force array n x 3
    output
    updated p,v,f
    '''
    dim = 3

    p_new = np.zeros((n,dim))
    ##p_old = np.zeros((n,dim))

    dt2l = (dt**2) / l
    l2dt = l / (2*dt)

    
    for i in xrange(n):
        p_new[i,0] = 2.0*p[i,0] - po[i,0] + f[i,0]*dt2l
        p_new[i,1] = 2.0*p[i,1] - po[i,1] + f[i,1]*dt2l
        p_new[i,2] = 2.0*p[i,2] - po[i,2] + f[i,2]*dt2l

    
        v[i,0] = (p_new[i,0] - po[i,0]) * l2dt 
        v[i,1] = (p_new[i,1] - po[i,1]) * l2dt 
        v[i,2] = (p_new[i,2] - po[i,2]) * l2dt 
        
        po[i,:] = np.copy(p[i,:])
        p[i,:] = np.copy(p_new[i,:])

    
def Update(t,n,dn,v,vi,vl,pe,c=False):
    '''
    input:
    t - timestep
    n - number of atoms
    dn - density
    v - velocity array n x 3
    vi - virial pressure term
    v. - volume
    pe - potential energy
    c - conversion parameters
    '''
    if c == False:
        if t == 0:
            print ("NVE MD Run: %i atoms " %(n))
            print "LJ units!"
            print "Timestep    Energy[LJ]   Temperature[LJ]    Pressure[LJ]"
            print "--------------------------------------------------------"

        dim = 3
        ke = 0.0
        for d in xrange(dim):
            ke += np.sum(np.square(v[:,d]))

        ke *= 0.5
        temp = (2.0/3.0) * ( ke / n)
        energy = ke + pe
        press = dn*temp + vi/(3*vl)
        
        print ("%i       %f        %f        %f " %(t,energy,temp,press))

    else:
        if t == 0:
            print ("NVE MD Run: %i atoms " %(n))
            print "None LJ units!"
            print "Timestep    Energy[eV]   Temperature[K]    Pressure[GPa]"
            print "--------------------------------------------------------"

        dim = 3
        ke = 0.0
        for d in xrange(dim):
            ke += np.sum(np.square(v[:,d]))

        ke *= 0.5
        temp = (2.0/3.0) * ( ke / n) 
        energy = (ke + pe) * c['e']
        press = (dn*temp + vi/(3*vl)) * (c['e'] / c['s']**3)
        press *= 160.2176487 # eV/A^3 to GPa
        temp = (2.0/3.0) * ( ke / n) * c['tp']

        print ("%i       %f        %f        %f " %(t,energy,temp,press))

    return None

if __name__ == "__main__":

    #RUN MD

    #Ar Params in eV and Angstroms
    epsilon = 0.0108319 
    sigma = 3.345
    mass = 39.950

    time_factor = sqrt( mass*sigma*sigma/epsilon)
    vol_factor = sigma*sigma*sigma
    temp_factor = epsilon / 8.6173324e-5
    
    #If you want output units in Angstrom,eV,etc.
    convert = {'s':sigma,'e':epsilon,'m':mass,
                  't':time_factor,'v':vol_factor,'tp':temp_factor}
    #Or in LJ
    #convert = False
    
    #Runtime Parameters
    output = 250
    timestep = 0.0001 
    runtime = 10000
    ncells = 3
    temperature = 1.00
    density = 0.8975

    
    #Initial Coordinates/Velocities
    natoms,positions = FCC(ncells)
    velocity = MBdist(natoms,temperature)
    
   
    #Only Cubic system
    volume = natoms/density
    length = volume**(1./3.)
    rcut = length/2.0

    #Get Initial Forces
    potengr,virial,forces = ForceCalc(length,natoms,
                                      rcut,positions)
 

    #Preemptive Update
    positions_old = np.zeros((natoms,3))
    PreIntegrate(natoms,timestep,length,positions,
                 positions_old,velocity,forces)

    #Time Loop
    print "--------------------------------------------------------"
    print ("Timestep size: %f" %(timestep))
    
    t = 0
    while t <= runtime:
        Integrate(natoms,timestep,length,positions,
                  positions_old,velocity,forces)

        #Data Output
        if (t % output) == 0:
            Update(t,natoms,density,velocity,virial,volume,potengr,convert)
            LAMMPSdump(natoms,t,length,positions,velocity,forces)

        #Force Calculation
        potengr,virial,forces = ForceCalc(length,natoms,rcut,positions)

        t += 1
        
    
    print "--------------------------------------------------------"
    print "NVE MD Run Complete "





if __name__ == '__main__':
    pass