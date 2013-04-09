# Python-based Guderley imploding shock calculations
# This is the post-processor
# A. Zylstra 2012/05/24

from Guderley import *
from Fusion import *
from Plasma import *
from Constants import *
import csv
import math
from scipy.integrate import quad

class Postproc:
    """A post-processor for Guderley calculations."""
    #Class variables
    g = 0
    dr = 2e-4
    rmin = 0
    rmax = 0.02
    dt = 5e-12
    tmin = -0.2e-9
    tmax = 0.5e-9
    name = 'Postproc'
    # arrays of physical variables for smoothing
    sni = numpy.zeros((1,1),dtype=numpy.float)
    sTi = numpy.zeros((1,1),dtype=numpy.float)
    #smoothing?
    smooth = 0

    def __init__(self,name,g):
        self.name = name
        d = os.path.dirname(name)
        if not os.path.exists(name):
            os.makedirs(name)
        self.g = g
        self.tmin += g.tc
        self.tmax = g.tFF()
        self.rmax = g.r0
        #get time steps
        self.dr = pow(10,-4)*float(input("Radial step (um): "))
        self.dt = pow(10,-12)*float(input("Time step (ps): "))
        #set up postproc data arrays
        nr = self.g.r0 / self.dr
        nt = (self.tmax - self.tmin) / self.dt
        self.sni = numpy.zeros((nr,nt),dtype=numpy.float)
        self.sTi = numpy.zeros((nr,nt),dtype=numpy.float)
        #ask for smoothing factor
        self.smooth = float(input("Smoothing factor (0-1): "))
        #read data in from g
        self.read_data()

    def read_data(self):
        """Read data from Guderley solution into the postprocessor, doing smoothing if necessary."""
        nr = self.g.r0 / self.dr
        nt = (self.tmax - self.tmin) / self.dt
        tempni = numpy.zeros((nr,nt),dtype=numpy.float)
        tempTi = numpy.zeros((nr,nt),dtype=numpy.float)
        for i in range(int(self.g.r0 / self.dr)):
            for j in range(int((self.tmax-self.tmin)/self.dt)):
                if (j*self.dt+self.tmin) == self.g.tc:
                    tempni[i,j] = 0.
                    tempTi[i,j] = 0.
                else:
                    tempni[i,j] = self.g.ni(i*self.dr, j*self.dt+self.tmin)
                    tempTi[i,j] = self.g.Ti(i*self.dr, j*self.dt+self.tmin)

        if self.smooth > 0:
            for i in range(int(self.g.r0 / self.dr)):
                for j in range(int((self.tmax-self.tmin)/self.dt)):
                    if (j*self.dt+self.tmin) == self.g.tc:
                        tempni[i,j] += 0.
                        tempTi[i,j] += 0.
                    else:
                        ni = self.g.ni(i*self.dr, j*self.dt+self.tmin)
                        Ti = self.g.Ti(i*self.dr, j*self.dt+self.tmin)
                        sigmar = self.IonMFP(i*self.dr, j*self.dt+self.tmin)
                        sigmat = self.Tauii(i*self.dr, j*self.dt+self.tmin)
                        kern = self.getKernel(i,j,sigmar,sigmat)
                        tempni += kern*ni*self.dr*self.dt
                        tempTi += kern*Ti*self.dr*self.dt
        self.sni = tempni
        self.sTi = tempTi

    def getKernel(self,indexr,indext,sigmar,sigmat):
        """Get a kernel for smoothing calculations."""
        nr = self.g.r0 / self.dr
        nt = (self.tmax - self.tmin) / self.dt
        kern = numpy.zeros((nr,nt),dtype=numpy.float)
        mur = indexr*self.dr
        mut = indext*self.dt+self.tmin
        #go 3 sigma out
        di = int(3*sigmar/self.dr)
        dj = int(3*sigmat/self.dt)
        for i in range(-di,di):
            for j in range(-dj,dj):
                if int(math.fabs(indexr + i)) < nr-1 and (j+indext) < nt-1 and (j+indext) >= 0 and int(math.fabs(indexr + i)) >= 0:
                    kern[int(math.fabs(indexr + i))][indext + j] += self.gauss(int(math.fabs(indexr + i))*self.dr,(indext+j)*self.dt+self.tmin, mur, sigmar, mut, sigmat)

        #Need to do a correction for the effective volume
        InitVol = 4*math.pi*pow(mur,2)*self.dr
        FinalVol = 0
        for i in range(-di,di):
            for j in range(-dj,dj):
                if int(math.fabs(indexr + i)) < nr-1 and (j+indext) < nt-1 and (j+indext) >= 0 and int(math.fabs(indexr + i)) >= 0:
                    FinalVol += kern[int(math.fabs(indexr + i))][indext + j]*4*3.1415*pow(int(math.fabs(indexr + i))*self.dr,2)*self.dr
        if FinalVol > 0:
            kern = kern * pow(InitVol/FinalVol,2) #one factor of Vi/Vf for normalization of kernel, the other factor for "ideal gas" expansion cooling & number conservation
        return kern

    def gauss(self, r,t,mur,sigmar,mut,sigmat):
        """Gaussian in r and t."""
        return (1/(2*math.pi*sigmar*sigmat))*math.exp(-pow(r-mur,2)/(2*pow(sigmar,2)))*math.exp(-pow(t-mut,2)/(2*pow(sigmat,2)))
                

    def write(self, t):
        """Write out at time t (s)."""
        output = csv.writer(open(os.path.join(self.name,'Snapshot.csv'),'w'))
        output.writerow( ["r (cm)", "u (cm/s)", "cs (cm/s)", "rho (g/cc)", "T (keV)", "P (GBar)", "Ion MFP (cm)", "tau_ii (s)", "LogL", "Debye Length (cm)", "Smoothed ni (1/cc)", "Smoothed Ti (keV)"] )
        for i in list(numpy.arange(self.rmin, self.rmax, self.dr)):
            output.writerow( [i, self.g.u(i,t), self.g.c(i,t), self.g.rho(i,t), self.g.T(i,t), self.g.P(i,t), self.IonMFP(i,t), self.Tauii(i,t), self.LogL(i,t), self.LambdaD(i,t), self.sni[i/self.dr,(t-self.tmin)/self.dt], self.sTi[i/self.dr,(t-self.tmin)/self.dt] ])

    
    def run(self):
        """Calculate the total yield (between tmin and tmax)."""
        YDD = 0
        YD3He = 0
        YHeHe = 0
        YDDMolvig = 0
        YD3HeMolvig = 0
        YHeHeMolvig = 0
        rateFile = csv.writer(open(os.path.join(self.name,'BurnRate.csv'),'w'))
        rateFile.writerow( ["t (s)", "DD (1/s)", "D3He (1/s)", "3He3He (1/s)"] )
        yieldFile = csv.writer(open(os.path.join(self.name,'Yield.csv'),'w'))
        MolvigYieldFile = csv.writer(open(os.path.join(self.name,'MolvigYield.csv'),'w'))
        TiDD = 0
        TiD3He = 0
        TiHeHe = 0
        TiFile = csv.writer(open(os.path.join(self.name,'Ti.csv'),'w'))
        
        #iterate over all time:
        for t in list(numpy.arange(self.tmin, self.tmax, self.dt)):
            #DD
            [dDD, dDDMolvig] = self.DDrate(t)
            YDD += dDD*self.dt
            YDDMolvig += dDDMolvig*self.dt
            TiDD += self.DDTirate(t)*self.dt
            #D3He
            [dD3He, dD3HeMolvig] = self.D3Herate(t)
            YD3He += dD3He*self.dt
            YD3HeMolvig += dD3HeMolvig*self.dt
            TiD3He += self.D3HeTirate(t)*self.dt
            #3He3He
            [d3He3He, d3He3HeMolvig] = self.HeHerate(t)
            YHeHe += d3He3He*self.dt
            YHeHeMolvig += d3He3HeMolvig*self.dt
            TiHeHe += self.HeHeTirate(t)*self.dt
            #output
            rateFile.writerow( [t, dDD, dD3He, d3He3He] )

        #If there is DD:
        if YDD > 0:
            TiDD = TiDD / YDD
            print("DD yield = " + '{:.2e}'.format(YDD))
            print("DD Molvig yield = " + '{:.2e}'.format(YDDMolvig))
            print("DD Ti = " + '{:.2f}'.format(TiDD))
            yieldFile.writerow( ["DD",YDD] )
            MolvigYieldFile.writerow( ["DD",YDDMolvig] )
            TiFile.writerow( ["DD",TiDD] )
        if YD3He > 0:
            TiD3He = TiD3He / YD3He
            print("D3He yield = " + '{:.2e}'.format(YD3He))
            print("D3He Molvig yield = " + '{:.2e}'.format(YD3HeMolvig))
            print("D3He Ti = " + '{:.2f}'.format(TiD3He))
            yieldFile.writerow( ["D3He",YD3He] )
            MolvigYieldFile.writerow( ["D3He",YD3HeMolvig] )
            TiFile.writerow( ["D3He",TiD3He] )
        if YHeHe > 0:
            TiHeHe = TiHeHe / YHeHe
            print("3He3He yield = " + '{:.2e}'.format(YHeHe))
            print("3He3He Molvig yield = " + '{:.2e}'.format(YHeHeMolvig))
            print("3He3He Ti = " + '{:.2f}'.format(TiHeHe))
            yieldFile.writerow( ["3He3He",YHeHe] )
            MolvigYieldFile.writerow( ["3He3He",YHeHeMolvig] )
            TiFile.writerow( ["3He3He",TiHeHe] )


    # ------------------------------------
    # Burn rate calculators
    # ------------------------------------
    def DDrate(self, t):
        """Calculate the DD burn rate at time t (s). Returns [normal, Molvig] rate."""
        fD = self.g.f1 # D fraction (atomic)
        if fD == 0:
            return [0,0]
        
        ret1 = 0
        ret2 = 0
        it = (t-self.tmin) / self.dt
        for r in range(int(self.g.rFF(t) / self.dr)):
            ret1 += DD(self.sTi[r,it])*pow(self.sni[r,it],2)*(fD*fD/2)*4*math.pi*pow(r*self.dr,2)*self.dr
            ret2 += self.svDDMolvig(r*self.dr,self.tmin+self.dt*it)*pow(self.sni[r,it],2)*(fD*fD/2)*4*math.pi*pow(r*self.dr,2)*self.dr
        return [ret1, ret2]
    def D3Herate(self, t):
        """Calculate the D3He burn rate at time t (s). Returns [normal, Molvig] rate."""
        fD = self.g.f1 # D fraction (atomic)
        f3He = self.g.f2 # 3He fraction (atomic)
        if fD*f3He == 0:
            return [0,0]
        
        ret1 = 0
        ret2 = 0
        it = (t-self.tmin) / self.dt
        for r in range(int(self.g.rFF(t) / self.dr)):
            ret1 += D3He(self.sTi[r,it])*pow(self.sni[r,it],2)*(fD*f3He)*4*math.pi*pow(r*self.dr,2)*self.dr
            ret2 += self.svD3HeMolvig(r*self.dr,self.tmin+self.dt*it)*pow(self.sni[r,it],2)*(fD*f3He)*4*math.pi*pow(r*self.dr,2)*self.dr
        return [ret1, ret2]
    def HeHerate(self, t):
        """Calculate the 3He3He burn rate at time t (s). Returns [normal, Molvig] rate."""
        f3He = self.g.f2 # 3He fraction (atomic)
        if f3He == 0:
            return [0,0]
        
        ret1 = 0
        ret2 = 0
        it = (t-self.tmin) / self.dt
        for r in range(int(self.g.rFF(t) / self.dr)):
            ret1 += HeHe(self.sTi[r,it])*pow(self.sni[r,it],2)*(f3He*f3He/2)*4*math.pi*pow(r*self.dr,2)*self.dr
            ret2 += self.svHeHeMolvig(r*self.dr,self.tmin+self.dt*it)*pow(self.sni[r,it],2)*(f3He*f3He/2)*4*math.pi*pow(r*self.dr,2)*self.dr
        return [ret1, ret2]

    # ------------------------------------
    # Burn-averaged Ti "rate" calculators
    # ------------------------------------
    def DDTirate(self, t):
        """Calculate the DD Ti 'burn rate' at time t (s)."""
        fD = self.g.f1 # D fraction (atomic)
        if fD == 0:
            return 0
        
        ret = 0
        it = (t-self.tmin) / self.dt
        for r in range(int(self.g.rFF(t) / self.dr)):
            ret += self.sTi[r,it]*DD(self.sTi[r,it])*pow(self.sni[r,it],2)*(fD*fD/2)*4*math.pi*pow(r*self.dr,2)*self.dr
        return ret
    def D3HeTirate(self, t):
        """Calculate the D3He burn rate at time t (s)."""
        fD = self.g.f1 # D fraction (atomic)
        f3He = self.g.f2 # 3He fraction (atomic)
        if fD*f3He == 0:
            return 0
        
        ret = 0
        it = (t-self.tmin) / self.dt
        for r in range(int(self.g.rFF(t) / self.dr)):
            ret += self.sTi[r,it]*D3He(self.sTi[r,it])*pow(self.sni[r,it],2)*(fD*f3He)*4*math.pi*pow(r*self.dr,2)*self.dr
        return ret
    def HeHeTirate(self, t):
        """Calculate the DD Ti 'burn rate' at time t (s)."""
        f3He = self.g.f2 # D fraction (atomic)
        if f3He == 0:
            return 0
        
        ret = 0
        it = (t-self.tmin) / self.dt
        for r in range(int(self.g.rFF(t) / self.dr)):
            ret += self.sTi[r,it]*HeHe(self.sTi[r,it])*pow(self.sni[r,it],2)*(f3He*f3He/2)*4*math.pi*pow(r*self.dr,2)*self.dr
        return ret

    # ------------------------------------
    # Print the trajectories
    # ------------------------------------
    def PrintTrajectories(self):
        """Print relevant trajectories in the problem."""
        t0 = self.g.tc-self.g.t0
        File = csv.writer(open(os.path.join(self.name,'Trajectories.csv'),'w'))
        File.writerow( ["t (s)", "rs (cm)", "us (cm/s)", "rFF (cm)", "rShell (cm)"] )
        for t in list(numpy.arange(t0, self.tmax, self.dt)):
            File.writerow( [t, self.g.rs(t), self.g.us(t), self.g.rFF(t), self.g.rShell(t)] )

    # ------------------------------------
    # Calculators for various plasma parameters
    # See formulary for reference
    # ------------------------------------        
    def IonMFP(self,r,t):
        """Calculate the ion mean free path (cm)."""
        return IonMFP( self.g.ni(r,t) , self.g.FuelZ , self.g.FuelA )
    def Tauii(self,r,t):
        """Calculate the ion-ion collision time (s)."""
        return Tauii( self.g.ni(r,t) , self.g.Ti(r,t) , self.g.FuelZ, self.g.FuelA )
    def LambdaD(self,r,t):
        """Calculate the Debye length (cm)."""
        return LambdaD( self.g.ne(r,t) , self.g.Te(r,t) )
    def uTherm(self,r,t):
        """Calculate the ion thermal velocity (cm/s)."""
        return uTherm( self.g.Ti(r,t) , self.g.FuelZ )
    def LogL(self,r,t):
        """Calculate the Coulomb logarithm."""
        return LogL( self.g.ni(r,t) , self.g.Ti(r,t) , self.g.FuelZ , self.g.FuelA )
    def PFermi(self,r,t):
        """Calculate the Fermi pressure at radius r (cm), time t (ns)."""
        return PFermi( self.g.ne(r,t) )

    # ------------------------------------
    # Calculate energy, entropy in the gas
    # ------------------------------------
    def ThermalEnergy(self,t):
        """Calculate the total thermal energy in the gas at time t (s). Returns Joules."""
        Energy = 0.
        it = (t-self.tmin) / self.dt
        for r in range(int(self.g.rShell(t) / self.dr)-1):
            Vol = 4*math.pi*pow(r*self.dr,2)*self.dr
            #ion thermal energy
            Energy += 1.5*self.sni[r,it]*Vol*kB*self.sTi[r,it]*11600.*1000.
            #electron thermal energy
            Energy += 1.5*self.g.ne(r*self.dr,it*self.dt)*Vol*kB*self.g.Te(r*self.dr,it*self.dt)*11600.*1000.
        return Energy*1e-7
    def alpha(self,t):
        """Calculate mass-weighted ratio of hydro pressure to Fermi pressure at time t."""
        it = (t-self.tmin) / self.dt #time index
        TotalMass = 1e-12
        alpha = 0.
        for r in range(int(self.g.rShell(t) / self.dr)-1):
            Vol = 4*math.pi*pow(r*self.dr,2)*self.dr
            Mass = Vol*self.sni[r,it]*self.g.FuelZ*mp
            Pf = max(self.PFermi(r*self.dr,t) , 1e-9)
            alpha += (self.g.P(r*self.dr,t) / Pf ) * Mass
            TotalMass += Mass
        alpha = alpha / TotalMass
        return alpha
    def EnergyEntropy(self):
        """Print the thermal energy versus time, and calculates entropy."""
        t0 = self.g.tc-self.g.t0
        File = csv.writer(open(os.path.join(self.name,'Energy.csv'),'w'))
        maxE = 0. #max energy dumped into gas after shock collapse
        maxA = 0. #max alpha in gas after shock collapse
        File.writerow( ["t (s)", "E (J)", "alpha"] )
        for t in list(numpy.arange(t0, self.tmax-self.dt, self.dt)):
            tempE = self.ThermalEnergy(t)
            tempA = self.alpha(t)
            if t > self.g.tc and tempE > maxE:
                maxE = tempE
            if t > self.g.tc and tempA > maxA:
                maxA = tempA
            File.writerow( [t, tempE, tempA ] )
        
        File = csv.writer(open(os.path.join(self.name,'SE_Summary.csv'),'w'))
        File.writerow( ["Thermal energy in gas after shock = " + str(maxE) + "J"] )
        File.writerow( ["Max gas P / Pf = " + str(maxA)] )
        
    
    # ------------------------------------
    # Make Lagrange plots
    # ------------------------------------
    def LagrangePlots(self):
        """Make Lagrange plots of initial gas material."""
        File = csv.writer(open(os.path.join(self.name,'Lagrange.csv'),'w'))
        dr = 30e-4 #30um spacing
        for r in list(numpy.arange(0,self.g.r0, dr)):
            self.g.rLagrangeInit(r)
            for t in list(numpy.arange(self.g.tc - self.g.t0, self.g.tFF()-self.dt, self.dt)):
                File.writerow( [ t , self.g.rLagrange(t) ] )
            File.writerow([])
            
        
    # -----------------------------------------------
    #      Molvig-Knudsen reactivities
    # -----------------------------------------------
    #Helper functions:
    def phi(self, theta, r, R):
        """Helper function for Molvig L calculations."""
        return ( math.pi/2.0 - theta - math.asin(r*math.sin(theta)/R) )
    def L(self, theta, r, R):
        """spherical geometry L helper function."""
        return math.sqrt( math.pow(R*math.sin(self.phi(theta,r,R)),2) + math.pow(r-R*math.cos(self.phi(theta,r,R)),2) )
    def integrand1(self, theta,r,R):
        """integrand for Leff calculation"""
        return ( math.sin(theta) / math.pow(self.L(theta,r,R),2) )
    def Leff(self, r,t):
        """Effective scale length in spherical geometry."""
        R = self.g.rShell(t) #shell radius
        if r >= R:
            return 0
        return math.sqrt( 2. / quad(self.integrand1, 0, math.pi, args=(r,R))[0] )
    # Knudsen number:
    def Nk(self, r, t, Z1, Z2):
        """Knudsen number for fusion with ion charges Z1,Z2."""
        return (math.sqrt(0.26) * pow(1000*11600*kB*self.g.Ti(r,t),2) / (math.pi*self.g.ni(r,t)*pow(Z1*Z2*e*e,2)*self.LogL(r,t)*self.Leff(r,t)) )
    # Knudsen distribution function:
    def fK(self, En,Nk,r,t):
        """Knudsen distribution function."""
        #normalized energy
        eps = En / self.g.Ti(r,t)
        #evaluated as prefactor and exponential
        p1 = 2 / math.sqrt( math.pi + Nk*math.pow(eps,1.5) )
        p2 = math.exp( -1.0*(eps+0.8*Nk*math.pow(eps,2.5)+0.32*math.pow(Nk*eps*eps,2))/(1 + 0.8*Nk*math.pow(eps,1.5)) )
        return p1*p2
    # DD reactivity:
    def integrand2(self, En, r, t):
        """Integrand for Molvig-Knudsen DD reactivity."""
        NkDD = self.Nk(r, t, 1, 1)
        return sigmaDDn(En)*c*math.sqrt(En/(1000*938))*self.fK(En,NkDD,r,t)*math.sqrt(En)
    def svDDMolvig(self, r,t):
        """Calculate DD reactivity (Molvig-Knudsen)."""
        R = self.g.rShell(t)
        if r > R or self.g.Ti(r,t) < 0.5:
            return 0
        return quad(self.integrand2, 0, 1000, args=(r,t))[0]
    # D3He reactivity:
    def integrand3(self, En, r, t):
        """Integrand for Molvig-Knudsen D3He reactivity."""
        NkDHe = self.Nk(r, t, 1, 2)
        return sigmaD3He(En)*c*math.sqrt(En/(1250*938))*self.fK(En,NkDHe,r,t)*math.sqrt(En)
    def svD3HeMolvig(self, r,t):
        """Calculate D3He reactivity (Molvig-Knudsen)."""
        R = self.g.rShell(t)
        if r > R or self.g.Ti(r,t) < 0.5:
            return 0
        return quad(self.integrand3, 0, 1000, args=(r,t))[0]
    # 3He3He reactivity:
    def integrand4(self, En, r, t):
        """Integrand for Molvig-Knudsen 3He3He reactivity."""
        NkHeHe = self.Nk(r, t, 2, 2)
        return sigmaHeHe(En)*c*math.sqrt(En/(1500*938))*self.fK(En,NkHeHe,r,t)*math.sqrt(En)
    def svHeHeMolvig(self, r,t):
        """Calculate 3He3He reactivity (Molvig-Knudsen)."""
        R = self.g.rShell(t)
        if r > R or self.g.Ti(r,t) < 0.5:
            return 0
        return quad(self.integrand4, 0, 1000, args=(r,t))[0]
    