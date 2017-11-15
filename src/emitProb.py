#
# Code for prototyping the calculation of the emission probability.
#
from math import *
from GaAsUtils import *

#
# imgageChargeR(s) -> calculates the R(s) function from Eq.(27)
#   in Jensen, Ch. 3 (2001).
#   - lower limit is set to 0
#   - upper limit is set to pi/2
#   - ng = 5 # of points used for the Gaussian quadrature algorithm
#   - x[i] - abscissas for the algorithm
#   - w[i] - weigths 
def imageChargeR(s) :
    ng = 5
    piOver4 = pi*0.25
    x = [0.148874338981631, 0.433395394129247, 0.679409568299024,
         0.865063366688985, 0.973906528517172]
    w = [0.295524224714753, 0.269266719309992, 0.219086362515982,
         0.149451349150581, 0.066671344308684]
    if (s > 60.) :
        return (pi*sqrt(s)/(16.0*s+4.0))
    else :
        Rs = 0.0
        for i in range(len(x)) :
            xi = x[i]
            wi = w[i]
            xp = piOver4*(xi+1.0)
            xm = piOver4*(-xi+1.0)
            cxp = cos(xp)
            sxp = sin(xp)
            sxpSq = sxp*sxp
            fxp = (cxp*cxp*sxpSq)/(sqrt(s+sxpSq))
            cxm = cos(xm)
            sxm = sin(xm)
            sxmSq = sxm*sxm
            fxm = (cxm*cxm*sxmSq)/(sqrt(s+sxmSq))
            Rs += wi*(fxp+fxm)
        return Rs*piOver4

#
# getTheta(V_max, E, F, Q, rmEff) : calculates theta(E) for the WKB theory
#   for the potential V(x) = V_max -F*x - Q/x from
#   Eq. (26) in Jensen, Ch. 3 (2001).
#
#   V_max - max of the potential barrier in eV
#   E - energy of the particle impinging on the barrier in eV
#   F - (particle charge)*(Electric field) in eV/Angstrom 
#   Q - (particle charge)^2/(16*pi*epsilon_0)*(Ks-1)/(Ks+1)
#         == alpha*hbar*c*(Ks-1)/(Ks+1)/4 in eV*Angstrom
#           with Ks the static dielectric constant
#   rmEff - reduced effective mass
#   Units: fs, eV, Angstrom, fundamental charge q = 1 =>
#          epsilon_0 = q^2/(4*pi*alpha*hbar*c)
#          mass in eV*fs^2/Angstrom^2
#
def getTheta(V_max, E, F, Q, rmEff) :
     h = V_max - E
     aa = 4.0*Q*F
     if ( (h - sqrt(aa)) > 0.0 ) :
          mass = m_e_eV_x_fs2_o_A2 * rmEff
          L_0 = sqrt(h*h - aa)/F
          x_0 = 2.0*Q/(h + F*L_0)
          ac = (2.0/h_bar_eVfs)*sqrt(2.0*F*mass*L_0*L_0*L_0)
          return (ac*imageChargeR(x_0/L_0))
     else : 
          return 0.0

#
#  getQ(rq, Ks) : returns (q^2/(16*pi*epsilon_0))*(Ks-1)/(Ks+1) in eV*Angstrom
#    rq -> particle charge in fundamental charge units
#    Ks -> Ks - static dielectric constant
#
def getQ(rq, Ks) :
     alpha = 1.0/137.036
     c_light_A_o_fs = c_light * 1e-5
     cof = rq*rq*alpha*h_bar_eVfs*c_light_A_o_fs/4
     if (Ks > 1e3) :
          return cof
     else :
          return (cof*(Ks-1.0)/(Ks+1.0))
#
# getTkApprox(V0, E, Fm, Q, rmEff) -> calculates the transmission propbability 
#   in the modified Airy approach with the approximation given 
#   by Eq. (10) in K. L. Jensen, J. Vac. Sci. technol. B21, 1528 (2003).
#   V(x) = V0 -F*x - Q/x and WKB-related quantities.
#   V0 - initial max of the potential barrier in eV
#   E - energy of the particle impinging on the barrier in eV
#   Fm - (particle charge)*(Electric field) in eV/m
#   Q - (particle charge)^2/(16*pi*epsilon_0)*(Ks-1)/(Ks+1)
#         == alpha*hbar*c*(Ks-1)/(Ks+1)/4 in eV*Angstrom
#           with Ks the static dielectric constant
#       for Q=0, the results is for a triangular potential.
#   rmEff - reduced effective mass (in units of the electron rest mass)
#   Units: fs, eV, Angstrom, fundamental charge q = 1 =>
#          epsilon_0 = q^2/(4*pi*alpha*hbar*c)
#          mass in eV*fs^2/Angstrom^2
#
def getTkApprox(V0, E, Fm, Q, rmEff) :
     F = Fm*1e-10 #  (particle charge)*(Electric field) in eV/Angstrom
     FQ4 = 4.0*F*Q
     sqrtFQ4 = sqrt(FQ4)
     L0 = sqrt(V0*V0-FQ4)/F
     x_min = 2.0*Q/(V0+F*L0)
     effF =  (V0 - sqrtFQ4)/(L0 - x_min)
     #effF = F
     mass = m_e_eV_x_fs2_o_A2 * rmEff
     rh_bar_Sq = h_bar_eV_x_fs*h_bar_eV_x_fs
     rf = pow((2.0*mass*effF)/(rh_bar_Sq), 1.0/3.0)
     kSq = 2.0*mass*E/rh_bar_Sq
     k0Sq = 2.0*mass*(V0-sqrtFQ4)/rh_bar_Sq
     k = sqrt(kSq)
     omega = (k0Sq - kSq)/(rf*rf)
     peta = 0.510697183837523
     neta = rf*pow((omega*omega + peta*peta), 0.25)
     if ((V0-sqrtFQ4 - E) > 0.0) :
         thetaE = getTheta(V0, E, F, Q, rmEff)
         #thetaEKJ = 2.*pow(omega, 1.5)/3.0
         #thetaELL = (2.0/3.0)*(sqrt(2.0*mass)/(effF*h_bar_eV_x_fs))*pow(V0-E-sqrtFQ4, 1.5)
         #print "t(E) = %13.5E, tKJ(E) = %13.5E, tLL(E) = %13.5E" % (thetaE, thetaEKJ, thetaELL)
         #thetaE = thetaELL
         tF = exp(2.0*thetaE)-0.25*(1.0-exp(-2.0*thetaE))
         Tk = 4.0*neta*k/(2.0*neta*k + (neta*neta+kSq)*tF)
     else :
         Tk = 4.0*neta*k/(2.0*neta*k + (neta*neta+kSq))
     # print "Tk = %13.5E, F = %13.5E, effF = %13.5E, sqrt(4FQ)/V0 = %13.5E" % (Tk, F, effF, sqrtFQ4/V0)
     return Tk
#
# getTkApproxWKB(V0, E, Fm, Q, rmEff) -> calculates the transmission 
#   propbability in the approximate WKB approach given by Eq. (18) in 
#   K. L. Jensen, J. Vac. Sci. technol. B21, 1528 (2003).
#   V(x) = V0 -F*x - Q/x and WKB-related quantities.
#   V0 - initial max of the potential barrier in eV
#   E - energy of the particle impinging on the barrier in eV
#   Fm - (particle charge)*(Electric field) in eV/m
#   Q - (particle charge)^2/(16*pi*epsilon_0)*(Ks-1)/(Ks+1)
#         == alpha*hbar*c*(Ks-1)/(Ks+1)/4 in eV*Angstrom
#           with Ks the static dielectric constant
#       for Q=0, the results is for a triangular potential.
#   rmEff - reduced effective mass (in units of the electron rest mass)
#   Units: fs, eV, Angstrom, fundamental charge q = 1 =>
#          epsilon_0 = q^2/(4*pi*alpha*hbar*c)
#          mass in eV*fs^2/Angstrom^2
#
def getTkApproxWKB(V0, E, Fm, Q, rmEff) :
    F = Fm*1e-10 #  (particle charge)*(Electric field) in eV/Angstrom
    FQ4 = 4.0*F*Q
    sqrtFQ4 = sqrt(FQ4)
    #
    mass = m_e_eV_x_fs2_o_A2 * rmEff
    rh_bar_Sq = h_bar_eV_x_fs*h_bar_eV_x_fs
    rf = 2.0*mass*F/rh_bar_Sq # different from the rf used in getTkApprox(.)
    kk = 2.0*mass*E/rh_bar_Sq
    k = sqrt(kk)
    xo = sqrt(Q/F)
    Eo = V0 - 1.25*sqrtFQ4
    aux = V0*V0-FQ4
    xmin = 2.0*Q
    if (aux > 0.0) :
        xmin /= (V0+sqrt(aux))
    else :
        xmin /= V0
    E3 = E*E*E
    aux = xo/xmin
    aux = (aux*aux)-1.0
    CofE = E3/(E3+(F*F*rh_bar_Sq/(128.0*mass))*(aux*aux))
    if ( E < Eo ) :
        thetaE = getTheta(V0, E, F, Q, rmEff)
        Tk = CofE/(1.0 + exp(2.0*thetaE))
    else :
        aux = sqrt(mass/rh_bar_Sq)
        theta_o = 1.1614992*aux*pow(Q*Q*Q/F, 0.25)
        dtheta = -2.4221121*aux*pow(Q/(F*F*F), 0.25)
        thetaE = theta_o + dtheta*(E - Eo)
        Tk = CofE/(1.0 + exp(2.0*thetaE))
    #
    return Tk

#
# getPotential(V0, Fm, Ks, E, nx) - calcualtes the potential:
#   V(x) = V0 -F*x - Q/x and WKB-related quantities.
#   V0 - initial max of the potential barrier in eV
#   E - energy of the particle impinging on the barrier in eV
#   Fm - (particle charge)*(Electric field) in eV/m
#   nx - number of points to use to discretize the interval [0, 1.01*x_max]
#   Q - (particle charge)^2/(16*pi*epsilon_0)*(Ks-1)/(Ks+1)
#         == alpha*hbar*c*(Ks-1)/(Ks+1)/4 in eV*Angstrom
#           with Ks the static dielectric constant
#   Units: fs, eV, Angstrom, fundamental charge q = 1 =>
#          epsilon_0 = q^2/(4*pi*alpha*hbar*c)
#          mass in eV*fs^2/Angstrom^2
#
def getPotential(V0, Fm, Ks, E, nx, xMax = 0) :
    Q = getQ(1.0, Ks)
    F = Fm*1e-10 #  (particle charge)*(Electric field) in eV/Angstrom
    V = V0 - E
    V_max = V0 - sqrt(4.0*F*Q)
    sq0 = sqrt(V0*V0 - 4.0*Q*F)
    L0 = sq0/F
    x_min = (V0 - sq0)/(2.0*F)
    x_max = x_min + L0
    sq = sqrt(V*V - 4.0*Q*F)
    L = sq/F
    xm = (V - sq)/(2.0*F)
    xp = xm + L
    print "%15s = %13.6E [%5s]" % ("Ks", Ks, "-")
    print "%15s = %13.6E [%5s]" % ("Q", Q, "eV*A")
    print "%15s = %13.6E [%5s]" % ("F", Fm, "eV/m")
    print "%15s = %13.6E [%5s]" % ("F", F, "eV/A")
    print "%15s = %13.6E [%5s]" % ("V0", V0, "eV")
    print "%15s = %13.6E [%5s]" % ("Vmax", V_max, "eV")
    print "%15s = %13.6E [%5s]" % ("V0-Vmax", V0-V_max, "eV")
    print "%15s = %13.6E [%5s]" % ("sqrt(4FQ)", sqrt(4.0*F*Q), "eV")
    print "%15s = %13.6E [%5s]" % ("sqrt(V0*V0-4FQ)", sqrt(V0*V0-4.0*F*Q), "eV")
    print "%15s = %13.6E [%5s]" % ("x_min", x_min, "A")
    print "%15s = %13.6E [%5s]" % ("x_max", x_max, "A")
    print "%15s = %13.6E [%5s]" % ("arg(V_max)", sqrt(Q/F), "A")
    print "%15s = %13.6E [%5s]" % ("x_max-x_min", L0, "A")
    print "%15s = %13.6E [%5s]" % ("xm", xm, "A")
    print "%15s = %13.6E [%5s]" % ("xp", xp, "A")
    print "%15s = %13.6E [%5s]" % ("xp-xm", L, "A")
    a0 = 0.0
    a1 = x_max
    if (xMax != 0) :
        a1 = xMax
    da = (a1-a0)/(float(nx))
    xc = []
    ycV = []
    ycE = [E for i in range(nx)]
    x = -da
    for c in range(nx) :
        x += da
        if (x < x_min ) :
            yv = 0.0
        else :
            yv = V0 - F*x - Q/x
            xc.append(x)
            ycV.append(yv)
    return (xc, ycV, ycE)

#####################################################################
