# Python-based Guderley imploding shock calculations
# This is the post-processor
# A. Zylstra 2012/01/27

from Guderley import *
from Reactivity import *
import csv
import math

class Postproc:
    """A post-processor for Guderley calculations."""
    #Class variables
    g = 0
    dr = 1e-4
    rmin = 0
    rmax = 0.02
    dt = 2e-12
    tmin = -0.2e-9
    tmax = 0.5e-9
    name = 'Postproc'

    def __init__(self,name,g):
        self.name = name
        d = os.path.dirname(name)
        if not os.path.exists(name):
            os.makedirs(name)
        self.g = g
        self.tmin += g.tc
        self.tmax += g.tc

    def write(self, t):
        """Write out at time t (s)."""
        output = csv.writer(open(os.path.join(self.name,'Snapshot.csv'),'w'))
        output.writerow( ["r (cm)", "u (cm/s)", "cs (cm/s)", "rho (g/cc)", "T (keV)", "P (GBar)"] )
        for i in list(numpy.arange(self.rmin, self.rmax, self.dr)):
            output.writerow( [i, self.g.u(i,t), self.g.c(i,t), self.g.rho(i,t), self.g.T(i,t), self.g.P(i,t)] )

    
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
        for r in list(numpy.arange(self.rmin, self.rmax, self.dr)):
            ret += DD(self.g.Ti(r,t))*pow(self.g.ni(r,t),2)*(fD*fD/2)*4*math.pi*pow(r,2)*self.dr
        return ret
    def D3Herate(self, t):
        """Calculate the D3He burn rate at time t (s)."""
        fD = self.g.f1 # D fraction (atomic)
        f3He = self.g.f2 # 3He fraction (atomic)
        if fD*f3He == 0:
            return 0
        
        ret = 0
        for r in list(numpy.arange(self.rmin, self.rmax, self.dr)):
            ret += D3He(self.g.Ti(r,t))*pow(self.g.ni(r,t),2)*(fD*f3He)*4*math.pi*pow(r,2)*self.dr
        return ret
    def HeHerate(self, t):
        """Calculate the 3He3He burn rate at time t (s)."""
        f3He = self.g.f2 # 3He fraction (atomic)
        if f3He == 0:
            return 0
        
        ret = 0
        for r in list(numpy.arange(self.rmin, self.rmax, self.dr)):
            ret += HeHe(self.g.Ti(r,t))*pow(self.g.ni(r,t),2)*(f3He*f3He/2)*4*math.pi*pow(r,2)*self.dr
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
        for r in list(numpy.arange(self.rmin, self.rmax, self.dr)):
            ret += self.g.Ti(r,t)*DD(self.g.Ti(r,t))*pow(self.g.ni(r,t),2)*(fD*fD/2)*4*math.pi*pow(r,2)*self.dr
        return ret
    def D3HeTirate(self, t):
        """Calculate the D3He burn rate at time t (s)."""
        fD = self.g.f1 # D fraction (atomic)
        f3He = self.g.f2 # 3He fraction (atomic)
        if fD*f3He == 0:
            return 0
        
        ret = 0
        for r in list(numpy.arange(self.rmin, self.rmax, self.dr)):
            ret += self.g.Ti(r,t)*D3He(self.g.Ti(r,t))*pow(self.g.ni(r,t),2)*(fD*f3He)*4*math.pi*pow(r,2)*self.dr
        return ret
    def HeHeTirate(self, t):
        """Calculate the 3He3He burn rate at time t (s)."""
        f3He = self.g.f2 # 3He fraction (atomic)
        if f3He == 0:
            return 0
        
        ret = 0
        for r in list(numpy.arange(self.rmin, self.rmax, self.dr)):
            ret += self.g.Ti(r,t)*HeHe(self.g.Ti(r,t))*pow(self.g.ni(r,t),2)*(f3He*f3He/2)*4*math.pi*pow(r,2)*self.dr
        return ret

    # ------------------------------------
    # Print the trajectories
    # ------------------------------------
    def PrintTrajectories(self):
        """Print relevant trajectories in the problem."""
        t0 = self.g.tc-self.g.t0
        File = csv.writer(open(os.path.join(self.name,'Trajectories.csv'),'w'))
        File.writerow( ["t (s)", "rs (cm)", "us (cm/s)", "rFF (1/s)"] )
        for t in list(numpy.arange(t0, self.tmax, self.dt)):
            File.writerow( [t, self.g.rs(t), self.g.us(t)] )
        
