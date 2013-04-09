# Python-based Guderley imploding shock calculations
# This is the post-processor
# A. Zylstra 2012/02/08

from Guderley import *
from Reactivity import *
import csv
import math

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
        InitVol = 4*3.1415*pow(mur,2)*self.dr
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
        return (1/(2*3.1415926*sigmar*sigmat))*math.exp(-pow(r-mur,2)/(2*pow(sigmar,2)))*math.exp(-pow(t-mut,2)/(2*pow(sigmat,2)))
                

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
        rateFile = csv.writer(open(os.path.join(self.name,'BurnRate.csv'),'w'))
        rateFile.writerow( ["t (s)", "DD (1/s)", "D3He (1/s)", "3He3He (1/s)"] )
        yieldFile = csv.writer(open(os.path.join(self.name,'Yield.csv'),'w'))
        TiDD = 0
        TiD3He = 0
        TiHeHe = 0
        TiFile = csv.writer(open(os.path.join(self.name,'Ti.csv'),'w'))
        
        #iterate over all time:
        for t in list(numpy.arange(self.tmin, self.tmax, self.dt)):
            #DD
            dDD = self.DDrate(t)
            YDD += dDD*self.dt
            TiDD += self.DDTirate(t)*self.dt
            #D3He
            dD3He = self.D3Herate(t)
            YD3He += dD3He*self.dt
            TiD3He += self.D3HeTirate(t)*self.dt
            #3He3He
            d3He3He = self.HeHerate(t)
            YHeHe += d3He3He*self.dt
            TiHeHe += self.HeHeTirate(t)*self.dt
            #output
            rateFile.writerow( [t, dDD, dD3He, d3He3He] )

        #If there is DD:
        if YDD > 0:
            TiDD = TiDD / YDD
            print("DD yield = " + '{:.2e}'.format(YDD))
            print("DD Ti = " + '{:.2f}'.format(TiDD))
            yieldFile.writerow( ["DD",YDD] )
            TiFile.writerow( ["DD",TiDD] )
        if YD3He > 0:
            TiD3He = TiD3He / YD3He
            print("D3He yield = " + '{:.2e}'.format(YD3He))
            print("D3He Ti = " + '{:.2f}'.format(TiD3He))
            yieldFile.writerow( ["D3He",YD3He] )
            TiFile.writerow( ["D3He",TiD3He] )
        if YHeHe > 0:
            TiHeHe = TiHeHe / YHeHe
            print("3He3He yield = " + '{:.2e}'.format(YHeHe))
            print("3He3He Ti = " + '{:.2f}'.format(TiHeHe))
            yieldFile.writerow( ["3He3He",YHeHe] )
            TiFile.writerow( ["3He3He",TiHeHe] )


    # ------------------------------------
    # Burn rate calculators
    # ------------------------------------
    def DDrate(self, t):
        """Calculate the DD burn rate at time t (s)."""
        fD = self.g.f1 # D fraction (atomic)
        if fD == 0:
            return 0
        
        ret = 0
        it = (t-self.tmin) / self.dt
        for r in range(int(self.g.rFF(t) / self.dr)):
            ret += DD(self.sTi[r,it])*pow(self.sni[r,it],2)*(fD*fD/2)*4*math.pi*pow(r*self.dr,2)*self.dr
        return ret
    def D3Herate(self, t):
        """Calculate the D3He burn rate at time t (s)."""
        fD = self.g.f1 # D fraction (atomic)
        f3He = self.g.f2 # 3He fraction (atomic)
        if fD*f3He == 0:
            return 0
        
        ret = 0
        it = (t-self.tmin) / self.dt
        for r in range(int(self.g.rFF(t) / self.dr)):
            ret += D3He(self.sTi[r,it])*pow(self.sni[r,it],2)*(fD*f3He)*4*math.pi*pow(r*self.dr,2)*self.dr
        return ret
    def HeHerate(self, t):
        """Calculate the 3He3He burn rate at time t (s)."""
        f3He = self.g.f2 # 3He fraction (atomic)
        if f3He == 0:
            return 0
        
        ret = 0
        it = (t-self.tmin) / self.dt
        for r in range(int(self.g.rFF(t) / self.dr)):
            ret += HeHe(self.sTi[r,it])*pow(self.sni[r,it],2)*(f3He*f3He/2)*4*math.pi*pow(r*self.dr,2)*self.dr
        return ret

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
        if self.g.ni(r,t) == 0:
            return 0
        return (2.28e7)*(1./self.g.FuelZ)*math.sqrt( self.g.FuelA/self.g.ni(r,t) )
    def Tauii(self,r,t):
        """Calculate the ion-ion collision time (s)."""
        if self.g.ni(r,t) == 0:
            return 0
        nu = (4.80e-8)*pow(self.g.FuelZ,4)*math.sqrt(1/(self.g.FuelA*pow(1000*self.g.Ti(r,t),3)))*self.g.ni(r,t)*self.LogL(r,t)
        return (1./nu)
    def LambdaD(self,r,t):
        """Calculate the Debye length (cm)."""
        if self.g.ni(r,t) == 0:
            return 0
        #kB=1.381e-16
        #e=4.803e-10
        #return math.sqrt(kB*self.g.Ti(r,t)*1000*11600 / (4*3.1415*self.g.ni(r,t)*e*e))
        return (7.43e2)*math.sqrt( 1000*self.g.Ti(r,t) / ((self.g.FuelZ+1)*self.g.ni(r,t)) )
    def uTherm(self,r,t):
        """Calculate the ion thermal velocity (cm/s)."""
        ret = (9.79e5)*math.sqrt( 1000*self.g.Ti(r,t)/self.g.FuelZ )
        return max(1,ret) #causes problems if u = 0
    def LogL(self,r,t): #See C.K. Li PRL 1993
        """Calculate the Coulomb logarithm."""
        if self.g.ni(r,t) == 0:
            return 0
        mp = 1.6726e-24 #g
        mr = mp*self.g.FuelA/2.
        e = 4.803e-10
        pperp = pow(e*self.g.FuelZ,2) / (mr * pow(self.uTherm(r,t),2) )
        hbar = 1.0546e-27 #cgs
        pmin = math.sqrt( pow(pperp,2) + pow(hbar/(2*mr*self.uTherm(r,t)),2) )
        return math.log( self.LambdaD(r,t) / pmin )
    def PFermi(self,r,t):
        """Calculate the Fermi pressure at radius r (cm), time t (ns)."""
        h = 6.626e-27 #erg*s
        me = 9.109e-28 #g
        ne = self.g.ne(r,t)
        Pf = ( pow(h,2) / (20*me) )*pow(3/3.1415,2/3)*pow(ne,5/3) #dyn
        return Pf *1e-6 * 1e-9 #Gbar

    # ------------------------------------
    # Calculate energy, entropy in the gas
    # ------------------------------------
    def ThermalEnergy(self,t):
        """Calculate the total thermal energy in the gas at time t (s). Returns Joules."""
        Energy = 0.
        it = (t-self.tmin) / self.dt
        kB=1.381e-16 #cgs
        for r in range(int(self.g.rShell(t) / self.dr)-1):
            Vol = 4*3.1415*pow(r*self.dr,2)*self.dr
            #ion thermal energy
            Energy += 1.5*self.sni[r,it]*Vol*kB*self.sTi[r,it]*11600.*1000.
            #electron thermal energy
            Energy += 1.5*self.g.ne(r*self.dr,it*self.dt)*Vol*kB*self.g.Te(r*self.dr,it*self.dt)*11600.*1000.
        return Energy*1e-7
    def alpha(self,t):
        """Calculate mass-weighted ratio of hydro pressure to Fermi pressure at time t."""
        it = (t-self.tmin) / self.dt #time index
        mp = 1.672e-24 #g
        TotalMass = 1e-12
        alpha = 0.
        for r in range(int(self.g.rShell(t) / self.dr)-1):
            Vol = 4*3.1415*pow(r*self.dr,2)*self.dr
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
