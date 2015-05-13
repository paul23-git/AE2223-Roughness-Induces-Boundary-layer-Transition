import numpy as np

beta = np.deg2rad(11.3)      #Oblique shock angle from table
theta = np.deg2rad(5.)
gamma = 1.4  
T01 = 579.
Tw = 290. 
cp = 1013. 
R = 287.05
S = 120.
Pr = 0.7

def stanton_variables(phot, Rex, q):
    pass

def st_Mn(Mhot):
    return Mhot * np.sin(beta)
def st_minf(Mhot):
    Mn = st_Mn(Mhot)
    return (1./np.sin(beta-theta))*np.sqrt((1.+(gamma-1.)/2.*Mn**2.)/(gamma*Mn**2.-(gamma-1.)/2.))

def st_rhototal_upstream(Phot):
    return Phot/R/T01

def st_rhostatic_upstream(rho_total_upstream, Mhot):
    return rho_total_upstream/(1.+(gamma-1.)/2.*Mhot**2)**(1./(gamma-1.))

def st_rhoinf(rho_total_upstream, Mhot):
    rho1 = st_rhostatic_upstream(rho_total_upstream, Mhot)
    Mn = st_Mn(Mhot)
    return rho1*((gamma+1.)*Mn**2.)/((gamma-1.)*Mn**2.+2.)

def st_Pstatic_upstream(Phot, Mhot):
    M1 = Mhot
    return Phot/((1.+(gamma-1.)/2.*M1**2)**(gamma/(gamma-1.)))

def st_Pinf(Phot, Mhot):
    p01 = Phot
    M1 = Mhot
    p1 = st_Pstatic_upstream(Phot, Mhot)
    Mn = st_Mn(Mhot)
    return p1*(1.+(2.*gamma/(gamma+1.))*(Mn**2.-1.))

def st_Tstatic_upstream(T01, Mhot):
    return T01/(1.+(gamma-1.)/2.*Mhot**2)

def st_Tinf(T01, Phot, Mhot, rho_total_upstream):    
    T1 = st_Tstatic_upstream(T01, Mhot)
    pinf = st_Pinf(Phot, Mhot)
    p1 = st_Pstatic_upstream(Phot, Mhot)
    rho1 = st_rhostatic_upstream(rho_total_upstream, Mhot)
    rhoinf = st_rhoinf(rho_total_upstream, Mhot)
    return T1*(pinf/p1)*(rho1/rhoinf)

def st_Ttotal_downstream(Tinf, Minf):
    return Tinf*(1+(gamma-1)/2*Minf**2)

def st_Muinf(mu0, Tinf, Ttot_down):
    return mu0*((Tinf/Ttot_down)**(3/2.))*(Ttot_down+S)/(Tinf+S) 


def stanton_experiment(phot, Rex,q, Mhot, Phot):
    Minf = st_minf(Mhot)
    
    
    rho_total_upstream = st_rhototal_upstream(Phot)
    rhoinf = st_rhoinf(rho_total_upstream, Mhot)
    
    Tinf = st_Tinf(T01, Phot, Mhot, rho_total_upstream)
    np.sqrt(gamma*R*Tinf)
    c = np.sqrt(gamma*R*Tinf)
    Vinf = Minf*c
    
    T02 = st_Ttotal_downstream(Tinf, Minf)
    C = (cp * rhoinf * Vinf * (T02-Tw) )
    sexp = q/C
    return sexp


def stanton_laminar(measurement, Rex):
    m = measurement
    
    mew0 = m.mu0
    Mhot = m.M
    Phot = m.pressure_pascal()
    rho_total_upstream = m.rho_total_upstream()
    Tw = m.Tw
    gamma = m.gamma
    R = m.R
    
    Tinf = st_Tinf(m.T0, Phot, Mhot, rho_total_upstream)
    Minf = st_minf(Mhot)
    pinf = st_Pinf(Phot, Mhot)
    rhoinf = st_rhoinf(rho_total_upstream, Mhot)
    T02 = st_Ttotal_downstream(Tinf, Minf)
    mewinf = st_Muinf(mew0, Tinf, T02)
    c = np.sqrt(gamma*R*Tinf)
    Vinf = Minf*c
    
    
    
    #Theory laminar values
    rlam=np.sqrt(Pr)   
    
    Treflam = Tinf*(1.+0.032*Minf**2.+0.58*(Tw/Tinf - 1.))            #Laminar reference Temperature
    mewreflam = mew0*((Treflam/T02)**(3/2.))*(T02 + S)/(Treflam+S)      #Viscosity
    rhoreflam = pinf/R/Treflam                                        #Reference density
    
    #Compressible laminar friction coefficient
    cflam = np.sqrt(rhoreflam*mewreflam/rhoinf/mewinf)*0.664/np.sqrt(Rex)
    
    #Theoretical laminar Stanton nubmer
    chlam = cflam*0.5*(Pr**(-2/3.))
    
    #Laminar values used for normalization of Stanton number
    Tawlam = Tinf*(1+rlam*((gamma-1.)/2.)*Minf**2.)
    
    qlam = chlam*cp*rhoinf*Vinf*(Tawlam-Tw)
    return qlam/cp/rhoinf/Vinf/(T02-Tw)
    return qlam


def stanton_turbulent(measurement, Rex):
    m = measurement
    
    mew0 = m.mu0
    Mhot = m.M
    Phot = m.pressure_pascal()
    rho_total_upstream = m.rho_total_upstream()
    Tw = m.Tw
    gamma = m.gamma
    R = m.R
    
    Tinf = st_Tinf(m.T0, Phot, Mhot, rho_total_upstream)
    Minf = st_minf(Mhot)
    pinf = st_Pinf(Phot, Mhot)
    rhoinf = st_rhoinf(rho_total_upstream, Mhot)
    T02 = st_Ttotal_downstream(Tinf, Minf)
    mewinf = st_Muinf(mew0, Tinf, T02)
    c = np.sqrt(gamma*R*Tinf)
    Vinf = Minf*c
    
    #Theory turbulent values
    rturb=Pr**(1/3.)  
    
    Trefturb = Tinf*(1.+0.035*(Minf**2.)+0.45*(Tw/Tinf - 1.))           #Turbulent reference Temperature
    mewrefturb = mew0*((Trefturb/T02)**(3/2.))*(T02 + S)/(Trefturb+S)     #Viscosity
    rhorefturb = pinf/R/Trefturb                                        #Reference density
    
    
    #Compressible turbulent friction coefficient
    cfturb = (rhorefturb/rhoinf)**(4/5.)*(mewrefturb/mewinf)**(1/5.)*0.0592/((Rex)**(1/5.))
    
    #Theoretical turbulent Stanton number
    chturb = cfturb*0.5*Pr**(-2/3.)
    
    #Turbulent values used for normalization
    Tawturb = Tinf*(1.+rturb*((gamma-1.)/2.)*Minf**2.)
    
    qturb = chturb*cp*rhoinf*Vinf*(Tawturb-Tw)
    return qturb/cp/rhoinf/Vinf/(T02-Tw)
    return qturb
      


def stanton(phot,Rex,q):

    #Input known values for given conditions
    S = 120.                     #Sutherlands constant for air
    R = 287.05                   #Ideal gas constant
    mew0 = 18.7*(10**(-6))       #Standard for air
    gamma = 1.4                  #ratio of specific heat
    Pr = 0.7                     #Prandtl number for air
    cp = 1013.                   #Standard for air
    beta = np.deg2rad(11.3)      #Oblique shock angle from table
    theta = np.deg2rad(5.)

    #Input from standard test conditions
    Tw = 290.                    #Given
    Mhot = 7.5                   #Given
    
    #####----- General calculations -----#####
    
    #Total conditions upstream oblique shock wave
    p01 = phot
    T01 = 579.
    rho01 = p01/R/T01
    
    #Upstream Mach number
    M1 = Mhot
    
    #Values for r
    rlam=np.sqrt(Pr)
    rturb=Pr**(1/3.)
    
    #Normal Mach number, upstream of the shockwave
    Mn = M1*np.sin(beta)
    
    #Static conditions upstream of the shockwave
    p1 = p01/((1.+(gamma-1.)/2.*M1**2)**(gamma/(gamma-1.)))  #Pressure
    
    rho1 = rho01/(1.+(gamma-1.)/2.*M1**2)**(1./(gamma-1.))   #Density
    
    T1 = T01/(1.+(gamma-1.)/2.*M1**2)                        #Temperature
    
    #Static conditions downstream of the shockwave
    Minf = (1./np.sin(beta-theta))*np.sqrt((1.+(gamma-1.)/2.*Mn**2.)/(gamma*Mn**2.-(gamma-1.)/2.))
    
    pinf = p1*(1.+(2.*gamma/(gamma+1.))*(Mn**2.-1.))
    
    rhoinf = rho1*((gamma+1.)*Mn**2.)/((gamma-1.)*Mn**2.+2.)
    
    Tinf = T1*(pinf/p1)*(rho1/rhoinf)
    
    #Total temperature downstream of the shockwave
    T02 = Tinf*(1+(gamma-1)/2*Minf**2)
      
    #speed of sound downstream
    c = np.sqrt(gamma*R*Tinf)
    
    #Flow velocity downstream
    Vinf = Minf*c
    
    #Viscosity infinity
    mewinf = mew0*((Tinf/T02)**(3/2.))*(T02+S)/(Tinf+S) 
    
    #-------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------
    
    #Theory laminar values
    Treflam = Tinf*(1.+0.032*Minf**2.+0.58*(Tw/Tinf - 1.))            #Laminar reference Temperature
    mewreflam = mew0*((Treflam/T02)**(3/2.))*(T02 + S)/(Treflam+S)      #Viscosity
    rhoreflam = pinf/R/Treflam                                        #Reference density
    
    #Compressible laminar friction coefficient
    cflam = np.sqrt(rhoreflam*mewreflam/rhoinf/mewinf)*0.664/np.sqrt(Rex)
    
    #Theoretical laminar Stanton nubmer
    chlam = cflam*0.5*(Pr**(-2/3.))
    
    #Laminar values used for normalization of Stanton number
    Tawlam = Tinf*(1+rlam*((gamma-1.)/2.)*Minf**2.)
    
    qlam = chlam*cp*rhoinf*Vinf*(Tawlam-Tw)
    
    #-------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------
    
    #Theory turbulent values
    Trefturb = Tinf*(1.+0.035*(Minf**2.)+0.45*(Tw/Tinf - 1.))           #Turbulent reference Temperature
    mewrefturb = mew0*((Trefturb/T02)**(3/2.))*(T02 + S)/(Trefturb+S)     #Viscosity
    rhorefturb = pinf/R/Trefturb                                        #Reference density
    
    
    #Compressible turbulent friction coefficient
    cfturb = (rhorefturb/rhoinf)**(4/5.)*(mewrefturb/mewinf)**(1/5.)*0.0592/((Rex)**(1/5.))
    
    #Theoretical turbulent Stanton number
    chturb = cfturb*0.5*Pr**(-2/3.)
    
    #Turbulent values used for normalization
    Tawturb = Tinf*(1.+rturb*((gamma-1.)/2.)*Minf**2.)
    
    qturb = chturb*cp*rhoinf*Vinf*(Tawturb-Tw)
      
    #-------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------
    
    #Normalised stanton numbers
    sexp = q/cp/rhoinf/Vinf/(T02-Tw)
    slam = qlam/cp/rhoinf/Vinf/(T02-Tw)
    sturb = qturb/cp/rhoinf/Vinf/(T02-Tw)
    
    #------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------

    #K-Value
    K = 100*(q-qlam)/(qturb-qlam)

    #-------------------------------------------------------------------------------------
    
    return sexp, slam, sturb, K


if (__name__ == "__main__"):
    pass