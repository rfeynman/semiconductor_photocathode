#
# Code for testing the implementation for the calculation of
# the emission probability using the two approximate models
# developed by Kevin Jensen.
#
from emitProb import *
from mProcData import *
from mplot1D import *
import sys
import cmath
from string import *

back_end = 'pdf'
back_end = 'eps'

#
# getTkTM is designed to calculate the transition probability
# using the transfer matrix method with the potential discretized
# in a step-wise way.
#
# input:
#   U [eV] - list of values for the potential at the emission surface
#            varying from 0 to N-1, the potential is considered
#            discretized at constant intervals delL
#   m [kg] - mass of the electron at the same spatial points as for
#            the potential for indeces 1 to N, the value of m[0] is
#            for the mass in the material from which electrons
#            are emitted (the effective mass)
#   dL [m] - length of the discretization interval for the potential
#   kPerp [1/m] - magnutude of the electron wavevector in the emission
#                 plane at the initial position (index 0)
#   kLong [1/m] - magnitude of the electron longitudinal wavevector at
#                 initial position (index 0). It is cosidered that this
#                 position is at the emission surface and the electron
#                 is moving towards the emission surface.
#   elecE [eV] - energy of the electron with which it attemps emission,
#                it is relative to the minimum of the conduction valley
#                the electron is on.
#   elecEm [eV] - the energy of the valley's minimum the electron is on
#                 relative to the origin from which electron energies on
#                 conduction bands are measured.
def getTkTM(U, m, dL, kPerp, kLong, elecE, elecEm = 0.0) :
    nU = len(U)
    if ((nU+1) != len(m)) :
        print "getTkTM(U, m, ...): the provided arguments U and m do not"
        print "  satisfy the required condition: len(m) == (len(U)+1)."
        print "  Exiting to the system ..."
        sys.exit(-1)
    kl = []
    kl.append(cmath.sqrt(kLong*kLong))#-1.32)
    kperpSq = kPerp*kPerp
    kFac = 2.0*eV_to_J/(h_bar*h_bar)
    totElecE = elecE + elecEm
    #print "totElecE = %13.5E" % (totElecE)
    #totElecE = kLong*kLong*h_bar*h_bar*0.5/(m[0]*eV_to_J)
    #print "totElecE = %13.5E" % (totElecE)
    for i in range(nU) :
        klt = cmath.sqrt((kFac*m[i+1]*(totElecE - U[i]))-kperpSq)
        kl.append(klt)
        #print "dL*kl.real = %15.5E, dL*k.imag = %15.5E" % (dL*klt.real, dL*klt.imag)

    # list to represent the total transfer matrix: 
    # trP[0] -> P_{11}, trP[1] -> P_{12}
    # trP[2] -> P_{21}, trP[3] -> P_{22}
    #print "kl = ", kl
    #print "m = ", m
    trP = [1.0, 0.0, 0.0, 1.0]
    #print U
    #print "kPerp = %13.5E" % (kPerp)
    #print "nU = %3d" % nU
    #print "U = ", U
    #print "m = ", m
    #print "kl = ", kl
    #print "dL = ", dL
    for i in range(nU) :
        kmf = m[i]*kl[i+1]/(m[i+1]*kl[i])
        kdL = kl[i]*dL
        #print "i = %3d: kmf.real = %15.5E, kmf.imag = %15.5E, kdL.real = %15.5E, kdL.imag = %15.5E" \
        #    % (i, kmf.real, kmf.imag, kdL.real, kdL.imag)
        ekLm = cmath.exp(-1j*kdL)
        ekLp = cmath.exp(1j*kdL)
        if (i == 0) :
            ekLm = cmath.exp(-1j*0.0)
            ekLp = cmath.exp(1j*0.0)
        p11 = 0.5*(1.0+kmf)*ekLm
        p12 = 0.5*(1.0-kmf)*ekLm
        p21 = 0.5*(1.0-kmf)*ekLp
        p22 = 0.5*(1.0+kmf)*ekLp
        # print "0.5*(1.0+kmf) = ", 0.5*(1.0+kmf)
        # print " 0.5*(1.0-kmf) = ",  0.5*(1.0-kmf)
        # print "ekLm, ekLp, 1/ekLp = ", ekLm, ekLp, 1.0/ekLp
        # print i, " p = ", [p11, p12, p21, p22]
        t11 = trP[0]*p11 + trP[1]*p21
        t12 = trP[0]*p12 + trP[1]*p22
        t21 = trP[2]*p11 + trP[3]*p21
        t22 = trP[2]*p12 + trP[3]*p22
        #print trP
        trP[0] = t11
        trP[1] = t12
        trP[2] = t21
        trP[3] = t22
        # print trP
        #absp11 = abs(trP[0])
        #print "P = %15.5E" % (1./(absp11*absp11))
    absP11 = abs(trP[0])
    T = 0.0
    if (kl[nU].real > 0.0) : 
        T = (1.0/(absP11*absP11))*(kl[nU]*m[0]/(kl[0]*m[nU]))
    Tr = abs(T)
    #print "pe = %13.5E" % pe
    #print trP
    #print "abs(trP[0]) = %13.5E " % abs(trP[0])
    #print "abs(trP[2]) = %13.5E " % abs(trP[2])
    #
    # a = 0.25*abs((1.0+kl[1]/kl[0])*(1.0+kl[2]/kl[1]))
    # b = 0.25*abs((1.0-kl[1]/kl[0])*(1.0-kl[2]/kl[1]))
    # Ta = abs(kl[2]/kl[0])/(a*a+b*b+2.0*a*b*cos(2.0*abs(kl[1]*dL)))
    # c = 0.25*abs((1.0-kl[1]/kl[0])*(1.0+kl[2]/kl[1]))
    # d = 0.25*abs((1.0+kl[1]/kl[0])*(1.0-kl[2]/kl[1]))
    # Ra = (c*c+d*d+2.0*c*d*cos(2.0*abs(kl[1]*dL)))/(a*a+b*b+2.0*a*b*cos(2.0*abs(kl[1]*dL)))
    #
    rP = abs(trP[2]/trP[0])
    Rt = rP*rP
    #print "Rt = %13.5E, Ra = %13.5E)" % (Rt, Ra)
    #print "1-Rt = %13.5E, 1-Ra = %13.5E)" % (1.0-Rt, 1.0-Ra)
    pe = 1.0 - rP*rP
    #print "1-R = %13.5E, T = %13.5E, Ta = %13.5E, 1 - Ra = %13.5E" % (1.0-rP*rP, abs(T), Ta, 1.0-Ra)
    print "1-R = %13.5E, T = %13.5E" % (pe, Tr)
    return Tr
    # if (pe > 1.0) :
    #     return 1.0
    # else : 
    #     return pe

def getVMKlongInRange(aV, aM, V0, Fm, Ks, totE, kperp, dL, minX = 0.0, rangeFlag = 0) :
    if (rangeFlag == 0) :
        return (aV, aM)
    Q = getQ(1.0, Ks)
    F = Fm*1e-10 #  (particle charge)*(Electric field) in eV/Angstrom
    longE = totE - (h_bar*h_bar*kperp*kperp*0.5/(0.36*m_e*eV_to_J))
    V_max = V0 - sqrt(4.0*F*Q)
    #print "Elong = %15.5E eV, V_max = %15.5E eV" % (longE, V_max)

    if (longE < V_max) :
        tv = V0 - longE
        maxX = ((tv+sqrt(tv*tv - 4.0*F*Q))*0.5/F) + 2.0*dL
    else :
        maxX = sqrt(Q/F) + 2.0*dL
    #print "maxX = %15.5E A" % maxX
    
    nV = []
    nM = [aM[0]]
    x = minX -dL
    for i in range(len(aV)) : 
        x += dL
        if ( x < maxX ) :
            nV.append(aV[i])
            nM.append(aM[i+1])
        else :
            break
    #print "len(V) = %5d, len(M) = %5d, len(nV) = %d, len(nM) %d" % (len(aV), len(aM), len(nV), len(nM))
    return (nV, nM)
        

def getTotEDiamond(klong, kperp) :
    m_l = 1.4*m_e
    m_t = 0.36*m_e
    ef_l = h_bar*h_bar*0.5/(m_l*eV_to_J)
    ef_t = h_bar*h_bar*0.5/(m_t*eV_to_J)
    return ((ef_l*klong*klong)+(ef_t*(kperp)*(kperp)))

def getECosSq(klong, kperp) :
    m_l = 1.4*m_e
    m_t = 0.36*m_e
    ef_l = h_bar*h_bar*0.5/(m_l*eV_to_J)
    ef_t = h_bar*h_bar*0.5/(m_t*eV_to_J)
    klSq = klong*klong
    kpSq = kperp*kperp
    totE = ((ef_l*klSq)+(ef_t*kpSq))
    #cosSq = klSq/(klSq + (m_l*m_l/(m_t*m_t))*kpSq)
    cosSq = klSq/(klSq + kpSq)
    return (totE*cosSq)

def getELongDiamond(klong_min, klong_max, kperp, nx) :
    m_l = 1.4*m_e
    m_t = 0.36*m_e
    ef_l = h_bar*h_bar*0.5/(m_l*eV_to_J)
    ef_t = h_bar*h_bar*0.5/(m_t*eV_to_J)
    en_L = []
    dx = (klong_max - klong_min)/(float(nx))
    xc = klong_min-dx
    for n in range(nx) :
        xc += dx
        enl = ef_l*xc*xc
        en_L.append(enl)
    return en_L

def getECosSqDiamond(klong_min, klong_max, kperp, nx) :
    m_l = 1.4*m_e
    m_t = 0.36*m_e
    ef_l = h_bar*h_bar*0.5/(m_l*eV_to_J)
    ef_t = h_bar*h_bar*0.5/(m_t*eV_to_J)
    kpsq = kperp*kperp
    e_t = ef_t*kpsq
    k0 = 2.0*pi*0.75E10/3.57
    en_L = []
    dx = (klong_max - klong_min)/(float(nx))
    xc = klong_min-dx
    for n in range(nx) :
        xc += dx
        klsq = xc*xc
        #klTot = xc+k0
        #klTotSq =  klTot*klTot
        cosSq = klsq/(klsq+kpsq)
        #cosSq = klTotSq/(klTotSq+kpsq)
        etot = ef_l*xc*xc + e_t
        en_L.append(etot*cosSq)
    return en_L

def compareTk(V0, Fm, numV, effMassFlag, rmEff, kLongMin, kLongMax, numKLong, kPerp, plotPotentialFlag = False) :
    tkFileSuffix = '_'
    #strFm = (split('%2.0f' % (Fm*1E-6)))[0]
    #plot_title = '$\chi = %3.0f$ meV, $F = %s$ MV/m' % (V0*1000, strFm)
    plot_title = 'exp stair step'
    #potential_file_name = 'potential_SS' % (V0*1000, strFm)
    #
    F = Fm*1e-10
    maxX = ((V0+sqrt(V0*V0 - 4.0*F*Q_diamond))*0.5/F)
    xMax = 1.0*maxX*1E-10
    print "xMax = ", xMax
    dL = xMax/numV
    energyConstVal = 0.120
    (xp, ssU, ye) = getPotential(V0, Fm, k_s, energyConstVal, numV, xMax = xMax*1E10)
    dump1DDataToFile(xp, ssU, "pot10p43.dat")
    #ssU = ssU[0:13]
    #print "dL = %13.5E" % (dL)
    #print "ssU = ", ssU
    #sys.exit(0)
    
    #
    ###dL = 8.0e-10 # in m
    ###ssU = [0.53, 0.25, 0.1] # in eV
    melec = [m_e for n in range(1+len(ssU))]
    if (effMassFlag) :
        melec[0] = 1.4*m_e # 0.067*m_e
        tkFileSuffix += 'MassJumpVsConst'
    # else :
    #     melec = [m_e for n in range(1+len(ssU))]
    #     tkFileSuffix += 'MassConst'
    kLongElec = 1.668E+10
    kPerpElec = 1.387E+09
    elecE = 6.002E-01
    emitP = getTkTM(ssU, melec, dL, kPerpElec, kLongElec, elecE, elecEm = 0.0)
    print "kPerp = %13.3E kLong = %13.3E elecE = %13.3E emitP = %13.3E" % (kPerpElec, kLongElec, elecE, emitP)
    #     
    kLongElec = 1.432E+10
    kPerpElec = 1.669E+09
    elecE = 3.531E-01
    emitP = getTkTM(ssU, melec, dL, kPerpElec, kLongElec, elecE, elecEm = 0.0)
    print "kPerp = %13.3E kLong = %13.3E elecE = %13.3E emitP = %13.3E" % (kPerpElec, kLongElec, elecE, emitP)
    #     
    kLongElec = 1.670E+10
    kPerpElec = 6.121E+08
    elecE = 4.404E-01
    emitP = getTkTM(ssU, melec, dL, kPerpElec, kLongElec, elecE, elecEm = 0.0)
    print "kPerp = %13.3E kLong = %13.3E elecE = %13.3E emitP = %13.3E" % (kPerpElec, kLongElec, elecE, emitP)
    #sys.exit(0)
    
    #
    dkLong = (kLongMax - kLongMin)/numKLong
    kLong = [(kLongMin + dkLong*nk) for nk in range(numKLong)]
    #
    kPerp = 0.0
    ###totE = [get_E(sqrt(kl*kl+kPerp*kPerp), eff_m = m_star, alpha = alpha_per_J)/eV_to_J for kl in kLong]
    totE0p0 = [getTotEDiamond(kl, kPerp) for kl in kLong]
    #totE = [get_E(sqrt(kl*kl+kPerp*kPerp), eff_m = m_e, alpha = 0.0)/eV_to_J for kl in kLong]
    #n = len(kLong) - 10
    #getTkTM(ssU, melec, dL, kPerp, kLong[n], totE[n], elecEm = 0.0)
    pEmitFX = [getTkTM(ssU, melec, dL, kPerp, kLong[n], totE0p0[n], elecEm = 0.0) for n in range(len(kLong))]
    #sys.exit(0)
    #
    kPerp = 0.1E10 # [1/m]
    #totE = [get_E(sqrt(kl*kl+kPerp*kPerp), eff_m = m_star, alpha = alpha_per_J)/eV_to_J for kl in kLong]
    totE0p1 = [getTotEDiamond(kl, kPerp) for kl in kLong]
    pEmitFX_kPerp0p1 = [getTkTM(ssU, melec, dL, kPerp, kLong[n], totE0p1[n], elecEm = 0.0) for n in range(len(kLong))]
    #
    kPerp = 0.15E10 # [1/m]
    #totE = [get_E(sqrt(kl*kl+kPerp*kPerp), eff_m = m_star, alpha = alpha_per_J)/eV_to_J for kl in kLong]
    totE0p15 = [getTotEDiamond(kl, kPerp) for kl in kLong]
    pEmitFX_kPerp0p15 = [getTkTM(ssU, melec, dL, kPerp, kLong[n], totE0p15[n], elecEm = 0.0) for n in range(len(kLong))]
    #
    kPerp = 0.25E10 # [1/m]
    #totE = [get_E(sqrt(kl*kl+kPerp*kPerp), eff_m = m_star, alpha = alpha_per_J)/eV_to_J for kl in kLong]
    totE0p25 = [getTotEDiamond(kl, kPerp) for kl in kLong]
    pEmitFX_kPerp0p3 = [getTkTM(ssU, melec, dL, kPerp, kLong[n], totE0p25[n], elecEm = 0.0) for n in range(len(kLong))]

    #
    # repeat the calculation for the const mass case.
    #
    melec[0] = melec[1]
    #melec = [1.4*m_e for n in range(1+len(ssU))]
    kPerp = 0.0
    totE = [getTotEDiamond(kl, kPerp) for kl in kLong]
    pEmitFX_Const = [getTkTM(ssU, melec, dL, kPerp, kLong[n], totE[n], elecEm = 0.0) for n in range(len(kLong))]
    #
    kPerp = 0.1E10 # [1/m]
    totE = [getTotEDiamond(kl, kPerp) for kl in kLong]
    pEmitFX_kPerp0p1_Const = [getTkTM(ssU, melec, dL, kPerp, kLong[n], totE[n], elecEm = 0.0) for n in range(len(kLong))]
    #
    kPerp = 0.15E10 # [1/m]
    totE = [getTotEDiamond(kl, kPerp) for kl in kLong]
    pEmitFX_kPerp0p15_Const = [getTkTM(ssU, melec, dL, kPerp, kLong[n], totE[n], elecEm = 0.0) for n in range(len(kLong))]
    #
    kPerp = 0.25E10 # [1/m]
    totE = [getTotEDiamond(kl, kPerp) for kl in kLong]
    pEmitFX_kPerp0p3_Const = [getTkTM(ssU, melec, dL, kPerp, kLong[n], totE[n], elecEm = 0.0) for n in range(len(kLong))]
    #
    if (kPerp > 0.0) : 
        plot_title = r'\rm{electrons in [100] valleys}' # % (kPerp*1E-10)
        tkVsKLongFileName = 'TkVsKLong' #'TkVsKLong_kP%2.0f' % (V0*1000, strFm, kPerp*1E-8)
    else :
        plot_title = r'\rm{$k_\perp = 0$ \AA$^{-1}$}' 
        tkVsKLongFileName = 'TkVsKLong_kP0' 
    tkVsKLongFileName += tkFileSuffix
    mplot1D(totE0p25, pEmitFX_kPerp0p3, scale_x=1e3, scale_y = 1, markershape = 'b-',
            root_file_name = '', open_flag = False)
    mplot1D(totE0p15, pEmitFX_kPerp0p15, scale_x=1e3, scale_y = 1, markershape = 'y-',
            root_file_name = '', open_flag = False)
    mplot1D(totE0p1, pEmitFX_kPerp0p1, scale_x=1e3, scale_y = 1, markershape = 'g-',
            root_file_name = '', open_flag = False)
    mplot1D(totE0p0, pEmitFX, scale_x=1e3, scale_y = 1, markershape = 'r-',
            root_file_name = '', open_flag = False)
    #
    mplot1D(totE0p25, pEmitFX_kPerp0p3_Const, scale_x=1e3, scale_y = 1, markershape = 'b--',
            root_file_name = '', open_flag = False)
    mplot1D(totE0p15, pEmitFX_kPerp0p15_Const, scale_x=1e3, scale_y = 1, markershape = 'y--',
            root_file_name = '', open_flag = False)
    mplot1D(totE0p1, pEmitFX_kPerp0p1_Const, scale_x=1e3, scale_y = 1, markershape = 'g--',
            root_file_name = '', open_flag = False)
    mplot1D(totE0p0, pEmitFX_Const, scale_x=1e3, scale_y = 1, markershape = 'r--',
            x_label = r'\rm{$E$ (meV)}',
            y_label = r'\rm{$T\left(E\left(k_\parallel, k_\perp\right)\right)$}',
            xlimt = (100.0, 1000.0), # slogy=True, 
            ylimt = (0.0, 1.0), 
            plot_legend = ('$k_\perp = 0.25$ \AA$^{-1}$', '$k_\perp = 0.15$ \AA$^{-1}$', '$k_\perp = 0.1$ \AA$^{-1}$', '$k_\perp = 0$ \AA$^{-1}$',),
            plot_location = 4, plot_title = plot_title,
            grid_flag = True, backend = back_end,
            root_file_name = tkVsKLongFileName, open_flag = True)
    #
    return

k_s = 5.7 # 12.9
Q_diamond = getQ(1.0, k_s)
print "%15s = %13.5E" % ("Q_diam eV*Angstrom", Q_diamond)
V0 = 0.3 # 0.55 # 0.28 # eV
Fm = 15e6 # 10e6 # 5e6
#V0 = 0.28 # 0.3 # 0.55 # 0.28 # eV
#Fm = 10.43e6 # 15.43e6 # 10e6 # 5e6
rmEff = 0.96 # 0.4259577 # 0.57

kLongMin = 0.0005E10
kLongMax = 0.5E10
numKLong = 50 # 200

#
# This is the case for direct aggreement with KJ's model when cos^2(theta) = 1 since kperp = 0
# The same effective mass value is used for both inside and outside of the material.
#  
kPerp = 0.0
effMassFlag = True # False
compareTk(V0, Fm, 500, effMassFlag, rmEff, kLongMin, kLongMax, numKLong, kPerp, plotPotentialFlag = True)
sys.exit(0)

#
# This is the case for direct aggreement with KJ's model when cos^2(theta) = 1 since kperp = 0
# However, the effective mass value jumps accross the emission surface interface - these plots
# are to consider the effect of the mass jumping accross the interface from its m_l = 1.4*m_e inside
# diamond to the m_e value in vaccuum. 
# 
kPerp = 0.0
effMassFlag = True
compareTk(V0, Fm, 500, effMassFlag, rmEff, kLongMin, kLongMax, numKLong, kPerp)

effMassFlag = True
kPerp = 0.15E10
compareTk(V0, Fm, 500, effMassFlag, rmEff, kLongMin, kLongMax, numKLong, kPerp)

effMassFlag = False
kPerp = 0.15E10
compareTk(V0, Fm, 500, effMassFlag, rmEff, kLongMin, kLongMax, numKLong, kPerp)
