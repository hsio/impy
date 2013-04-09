# Guderley calculation
# This code uses CGS units
#
# Author: Alex Zylstra
# Date: 2012/02/07

import math
import numpy
import scipy.optimize
import scipy.integrate
import scipy.interpolate
import csv
import os

class Guderley:
    """A wrapper class for a Guderley simulation."""
    # -----------------------------------------------------------
    # Global variables within this class
    # -----------------------------------------------------------
    # user specified parameters (default values here)
    r0 = 0.08
    rho0 = 0.005
    tc = pow(10,-9)
    xish = 400
    xirsh = xish*.74
    f1 = 0.5
    f2 = 0.5
    eiCoup = 0.
    name = "Guderley"

    # constants
    gamma = 5./3.
    n = 3.
    alpha = 0.68838
    Ion1 = "D"
    Ion1A = 2.
    Ion1Z = 1.
    Ion2 = "3He"
    Ion2A = 3.
    Ion2Z = 2.
    t0 = 0.
    FuelA = 0
    FuelZ = 0

    # Some physical constants
    Na = 6.022*pow(10,23)
    e = 1.602*pow(10,-19)
    AMU = 1.66053873*pow(10,-24)
    Erg2eV = 1/(6.2415*pow(10,11))
    P2GBAR = pow(10,-10)*1.0/(1.01325*pow(10,5))

    # Some Guderley constants
    gKappa = 0.
    gLambda = 1./alpha - 1.
    gEpsilon = gKappa*(gamma-1)+2*gLambda
    gMu = 2/(gamma-1)
    gBeta = n - gMu*gLambda
    gNu = n*gamma + gKappa - 2*gLambda

    #Singular points in UC plane
    P1=[0,1]
    P2=[0,0]
    P3=[0,0]
    P4=[0,0]
    P5=[0,0]
    P6=[0,0]
    #Other points
    A = [math.sqrt(2*gamma*(gamma-1))/(gamma+1) , 2/(gamma+1)]
    S1=[10,0]
    S2=[0,0]

    #Trajectories in UC plane
    UincList = []
    UincInt = 0
    UouterInt = 0
    UcentralInt = 0

    #Shell trajectory
    rShellInt = 0

    #List of points for C vs xi
    CvsXiShock = []
    CvsXiShockInt = 0
    CvsXiCentral = []
    CvsXiCentralInt = 0
    CvsXiOuter = []
    CvsXiOuterInt = 0
    GA = (gamma+1)/(gamma-1)
    K3 = GA*pow(alpha*A[0],-1.*gMu*(n+gKappa)/gBeta)*pow(1-A[1],-1.*((gKappa+gMu*gLambda)/gBeta))
    
    # ------------------------------------
    # Initialization
    # ------------------------------------
    def __init__(self, name):
        self.name = name
        #create a directory named name if it doesn't exist
        d = os.path.dirname(name)
        if not os.path.exists(name):
            os.makedirs(name)
        
        self.runGuderley()
    
    # ------------------------------------
    # Calculate the Guderley solution
    # ------------------------------------
    def runGuderley(self):
        """Top-level method to do the Guderley calculation."""
        self.doConfig()
        self.calcSingPoints()
        self.UCTrajectories()
        self.CalcUCG()
        self.rShellInit()
        return    

    # ------------------------------------
    # Set up the problem
    # ------------------------------------
    def doConfig(self):
        """Configure the top-level parameters."""
        self.r0 = pow(10,-4.)*float(input("Radius (um): "))
        #self.rho0 = pow(10,-3)*float(input("Density (mg/cm3): "))
        self.tc = pow(10,-9)*float(input("Collapse time (ns): "))
        self.xish = pow(10,-4)*pow(10,9*self.alpha)*float(input("Shock strength (um/ns^a): "))
        self.xirsh = 0.740 * self.xish
        #fuel info
        self.f1 = float(input("D2 fill (atm): "))
        self.rho0 = self.f1*0.08988*2/1000
        self.f2 = float(input("3He fill (atm): "))
        self.rho0 += self.f2*0.1786*(3/4)/1000
        temp = 2*self.f1 + self.f2 #in atm
        self.f1 = self.f1*2 / temp
        self.f2 = self.f2 / temp
        self.FuelA = self.f1*self.Ion1A + self.f2*self.Ion2A
        self.FuelZ = self.f1*self.Ion1Z + self.f2*self.Ion2Z
        self.t0 = pow(self.r0/self.xish,1./self.alpha)
        self.eiCoup = float(input("Rygg coupling parameter: "))
        
        return   

    # ------------------------------------
    # Calculate shock trajectories and velocities
    # ------------------------------------
    def rs(self,t):
        """Calculate the shock position at time t (s). Returns rs in cm."""
        strength = self.xish
        if t > self.tc:
            strength = self.xirsh
        return strength * pow( math.fabs(t-self.tc) , self.alpha )
    def us(self,t):
        """Calculate the shock velocity at time t (s). Returns us in cm/s."""
        strength = self.xish
        if t > self.tc:
            strength = self.xirsh
        return strength * self.alpha * pow( math.fabs(t-self.tc) , self.alpha-1 )
    def rFF(self,t):
        """Calculate the 'free fall' radius at time t (s). Returns rFF in cm."""
        uFF = self.u(self.r0 , (self.tc - self.t0 + 1e-12) )
        ret = max( 0. , (self.r0 - (t-(self.tc-self.t0))*uFF) )
        return ret
    def tFF(self):
        """Calculate when the 'free fall' mass hits the origin. Returns tFF in s."""
        uFF = self.u(self.r0 , (self.tc - self.t0 + 1e-12) )
        return ( self.r0/uFF + (self.tc-self.t0) )
    def rShellInit(self):
        """Do initial numerical integration for shell position vs time calculations."""
        times = list(numpy.arange(self.tc-self.t0,self.tFF(),1e-12))
        rList = scipy.integrate.odeint(self.uLab, self.r0, times)
        #fix list formatting
        rList2 = []
        for i in rList:
            rList2.append(i[0])
        #interpolate
        self.rShellInt = scipy.interpolate.interp1d(times, rList2)
    def rShell(self,t):
        """Shell position at time t (s). Returns rShell in cm."""
        return self.rShellInt(t)[()]

    # ------------------------------------
    # Self-similarity coordinate
    # ------------------------------------
    def xi(self, r, t):
        """Calculate the self-similarity coordinate as a function of r (cm) and t (s)."""
        return (r/self.r0) * pow( math.fabs( (t-self.tc)/self.t0 ) , -1.*self.alpha )
        
    # ------------------------------------
    # Calculate the U-C plane singular points
    # ------------------------------------
    def U23(self,x):
        return (self.n-1)*self.gamma*pow(x,2)+(self.gKappa-2*self.gLambda-self.gamma*(self.n-1-self.gLambda))*x-(self.gKappa-2*self.gLambda)

    def calcSingPoints(self):
        """Calculate the singular points in the UC plane."""
        U2 = scipy.optimize.fsolve(self.U23,0.4)[0]
        U3 = scipy.optimize.fsolve(self.U23,0.6)[0]
        self.P2 = [1-U2,U2]
        self.P3 = [1-U3,U3]
        self.P5 = [(math.sqrt(self.n)/(self.alpha*(self.n+self.gMu)))/math.sqrt(1.+(self.gKappa+self.gMu*self.gLambda)*(self.n+self.gMu)/((self.gMu+2)*(self.n-self.gMu*self.gLambda))) , self.gMu/(self.alpha*(self.n+self.gMu))]
        self.P6 = [10. , -(self.gKappa-2*self.gLambda)/(self.n*self.gamma)]
        
        #write singular points to file for reference
        singWriter = csv.writer(open(os.path.join(self.name,'SingularPoints.csv'), 'w'))
        singWriter.writerow(self.P1)
        singWriter.writerow(self.P2)
        singWriter.writerow(self.P3)
        singWriter.writerow(self.P4)
        singWriter.writerow(self.P5)
        singWriter.writerow(self.P6)
        
        return    

    # ------------------------------------
    # Some helper functions
    # ------------------------------------
    def Delta(self, x1, x2):
        return pow(x2,2.) - pow(1.-x1,2.)
    def Delta1(self, x1, x2):
        return x1*(1.-x1)*(1./self.alpha-x1)-pow(x2,2.)*(self.n*x1+(self.gKappa-2*self.gLambda)/self.gamma)
    def Delta2(self, x1, x2):
        return x2*((1.-x1)*(1/self.alpha-x1)+(x1/self.gMu)*(self.gLambda+(self.n-1)*(x1-1))-pow(x2,2.)+(self.gEpsilon/(2.*self.gamma))*pow(x2,2.)/(x1-1.))
    def dUdC(self, x1, x2):
        return self.Delta1(x1,x2) / self.Delta2(x1,x2)
    def ShockJump(self, x):
        """Shock jump in the UC plane. x is a 2D point."""
        ret1 = 2*self.gamma*(self.gamma-1)*pow((1-x[1]),2)/pow(self.gamma+1,2)+pow(x[0],2)*(1-2*pow(((self.gamma-1)/(self.gamma+1)),2)-2*(self.gamma-1)*pow((x[0]/(1-x[1])),2)/pow(self.gamma+1,2))
        ret2 = 1-(1-x[1])*((self.gamma-1)/(self.gamma+1)+2*pow((x[0]/(1-x[1])),2)/(self.gamma+1))
        return [ret1 , ret2]
    def ShockJumpErr1(self, x):
        return math.fabs(self.ShockJump([float(x),float(self.UCcentral(x))])[1] - self.UCouter(x))
    def ShockJumpErr2(self, x):
        return math.fabs(self.ShockJump([float(x),float(self.UCouter(x))])[1] - self.UCcentral(x))

    # ------------------------------------
    # Calculate trajectories in the UC plane
    # ------------------------------------
    def UCTrajectories(self):
        """Solve the Guderley solution for UC trajectories"""
        dInt = 0.0005
        
        #Incoming shock
        #Find the UC trajectory for shock between A and P3
        intRange = list(numpy.arange( self.A[0], self.P3[0], -1.0*dInt))
        temp = scipy.integrate.odeint( self.dUdC , self.A[1] , intRange )
        Uinc1 = []
        for i in range(0,len(temp)):
            Uinc1.append([ intRange[i] , temp[i,0] ])
        Uinc1.reverse() #ascending order
        
        #Find the UC trajectory for the shock between P4 and P3
        intRange = list(numpy.arange(self.P4[0]+dInt, self.P3[0]-0.02, dInt))
        temp = scipy.integrate.odeint( self.dUdC , self.P4[1]+dInt , intRange )
        Uinc2 = []
        for i in range(0,len(temp)):
            Uinc2.append([ intRange[i] , temp[i,0] ])
        
        #Stitch the whole trajectory together
        self.UincList = [self.P4] + Uinc2 + [self.P3] + Uinc1
        # set up interpolation
        UincListx = []
        UincListy = []
        for i in self.UincList:
            UincListx.append(i[0])
            UincListy.append(i[1])
        self.UincInt = scipy.interpolate.interp1d(UincListx, UincListy)
            
        #Outer flow region
        #numerical integration
        intRange = list(numpy.arange( dInt , 1.1 , dInt ))
        temp = scipy.integrate.odeint( self.dUdC , self.P4[0]-dInt , intRange )
        Uouterflow1 = []
        for i in range(0,len(temp)):
            Uouterflow1.append([ intRange[i] , temp[i,0] ])
        # set up interpolation
        x = []
        y = []
        for i in Uouterflow1:
            x.append(i[0])
            y.append(i[1])
        self.UouterInt = scipy.interpolate.interp1d(x, y)

        #Central flow region
        intRange = list(numpy.arange( self.P6[0] , 1-dInt , -1.*dInt ))
        temp = scipy.integrate.odeint( self.dUdC , self.P6[1] , intRange )
        Ucentralflow1 = []
        for i in range(0,len(temp)):
            Ucentralflow1.append([ intRange[i] , temp[i,0] ])
        # set up interpolation
        x = []
        y = []
        for i in Ucentralflow1:
            x.append(i[0])
            y.append(i[1])
        x.reverse()
        y.reverse()
        self.UcentralInt = scipy.interpolate.interp1d(x, y)
        

        #Find the shock jump points
        self.S1 = [10,0]
        self.S2 = [0,0]
        S2C = scipy.optimize.fmin(self.ShockJumpErr1 , 1.089, disp=0)[0]
        S1C = self.ShockJump([S2C,self.UCcentral(S2C)])[0]
        self.S1 = [S1C, self.UCouter(S1C)]
        self.S2 = [S2C, self.UCcentral(S2C)]

        #Output the trajectories to a file
        trajWriter = csv.writer(open(os.path.join(self.name,'UC_Trajectories_Inc.csv'), 'w'))
        for i in list(numpy.arange(dInt, self.A[0], dInt)):
            trajWriter.writerow([i,self.UCinc(i)])
        trajWriter = csv.writer(open(os.path.join(self.name,'UC_Trajectories_Outer.csv'), 'w'))
        for i in list(numpy.arange(dInt, self.S1[0], dInt)):
            trajWriter.writerow([i,self.UCouter(i)])
        trajWriter = csv.writer(open(os.path.join(self.name,'UC_Trajectories_Central.csv'), 'w'))
        for i in list(numpy.arange(self.S2[0]+dInt, self.P6[0], dInt)):
            trajWriter.writerow([i,self.UCcentral(i)])
        
        return

    # ------------------------------------
    # Interpolating functions for UC trajectories
    # ------------------------------------
    def UCinc(self, x):
        """U(C) for the incoming shock."""
        if type(x) is list:
            x = x[0]
        return self.UincInt(x)[()]
    def UCouter(self, x):
        """U(C) for the outer flow."""
        if type(x) is list:
            x = x[0]
        if x < self.S1[0] and x > 0 and x < 1.1:
            return self.UouterInt(x)[()]
        return 0
    def UCcentral(self, x):
        """U(C) for the central flow."""
        if type(x) is list:
            x = x[0]
        if x > self.S2[0] and x < self.P6[0] and x > 1:
            return self.UcentralInt(x)[()]
        return 0

    # ------------------------------------
    # Initial calculation of U,C,G vs self-similarity index
    # ------------------------------------
    def CalcUCG(self):
        """Calculate U,C,G vs xi."""
        #Incoming shock
        dXi=0.005
        XiMin=1
        XiMax=100
        self.CvsXiShock.append([XiMin,self.A[0]])
        for i in range(1, int(1 + (XiMax-XiMin)/dXi)):
            Xii = self.CvsXiShock[i-1][0]
            Ci = self.CvsXiShock[i-1][1]
            self.CvsXiShock.append([Xii+dXi,Ci+self.Delta2(self.UCinc(Ci),Ci)/self.Delta(self.UCinc(Ci),Ci)*(math.log(Xii+dXi) - math.log(Xii))])
        #interpolation
        x = []
        y = []
        for i in self.CvsXiShock:
            x.append(i[0])
            y.append(i[1])
        self.CvsXiShockInt = scipy.interpolate.interp1d(x, y)

        #Central flow region
        dXi = 0.001
        XiMin = dXi
        XiMax = 0.74
        self.CvsXiCentral.append([XiMax, self.S2[0]])
        for i in range(1, int(1+(XiMax-XiMin)/dXi)):
            Xii = self.CvsXiCentral[i-1][0]
            Ci = self.CvsXiCentral[i-1][1]
            self.CvsXiCentral.append([Xii-dXi,Ci-self.Delta2(self.UCcentral(Ci),Ci)/self.Delta(self.UCcentral(Ci),Ci)*(math.log(Xii+dXi) - math.log(Xii))])
        self.CvsXiCentral.reverse()
        #interpolation
        x = []
        y = []
        for i in self.CvsXiCentral:
            x.append(i[0])
            y.append(i[1])
        self.CvsXiCentralInt = scipy.interpolate.interp1d(x, y)

        #Outer flow region
        dXi = 0.005
        XiMin = 0.74
        XiMax = 100
        self.CvsXiOuter.append([XiMin,self.S1[0]])
        for i in range(1, int(1+(XiMax-XiMin)/dXi)):
            Xii = self.CvsXiOuter[i-1][0]
            Ci = self.CvsXiOuter[i-1][1]
            self.CvsXiOuter.append([Xii+dXi,Ci+self.Delta2(self.UCouter(Ci),Ci)/self.Delta(self.UCouter(Ci),Ci)*(math.log(Xii+dXi) - math.log(Xii))])
        #interpolation
        x = []
        y = []
        for i in self.CvsXiOuter:
            x.append(i[0])
            y.append(i[1])
        self.CvsXiOuterInt = scipy.interpolate.interp1d(x, y)
                          
        return
                               
    # ------------------------------------
    # U,C,G vs self-similarity index
    # ------------------------------------
    #incoming shock
    def Cshock(self, x):
        """C(xi) for the incoming shock."""
        if x > 1 and x < 100:
            return self.CvsXiShockInt(x)[()]
        return 0
    def Ushock(self, x):
        """U(xi) for the incoming shock."""
        if x > 1 and x < 100:
            return self.UCinc(self.CvsXiShockInt(x)[()])
        return 0
    def Gshock(self, x):
        """G(xi) for the incoming shock."""
        if x > 1 and x < 100:
            return self.K3*pow(self.alpha*pow(x,1/self.alpha)*self.Cshock(x),self.gMu*(self.n+self.gKappa)/self.gBeta)*pow(1-self.Ushock(x),(self.gKappa+self.gMu*self.gLambda)/self.gBeta)
        return 0
    #central flow region
    def Ccentral(self, x):
        """C(xi) for the central flow region."""
        if x > 0 and x <= 0.74:
            return self.CvsXiCentralInt(x)[()]
        return 0
    def Ucentral(self, x):
        """U(xi) for the central flow region."""
        if x > 0 and x <= 0.74:
            return self.UCcentral(self.CvsXiCentralInt(x)[()])
        return 0
    def Gcentral(self, x):
        """G(xi) for the central flow region."""
        if x > 0 and x <= 0.74:
            return self.K3*pow(self.alpha*pow(x,1/self.alpha)*self.Ccentral(x),self.gMu*(self.n+self.gKappa)/self.gBeta)*pow(1-self.Ucentral(x),(self.gKappa+self.gMu*self.gLambda)/self.gBeta)
        return 0
    #outer flow region
    def Couter(self, x):
        """C(xi) for the outer flow region."""
        if x > 0.74 and x < 100:
            return self.CvsXiOuterInt(x)[()]
        return 0
    def Uouter(self, x):
        """U(xi) for the outer flow region."""
        if x > 0.74 and x < 100:
            return self.UCouter(self.CvsXiOuterInt(x)[()])
        return 0
    def Gouter(self, x):
        """G(xi) for the outer flow region."""
        if x > 0.74 and x < 100:
            return self.K3*pow(self.alpha*pow(x,1/self.alpha)*self.Couter(x),self.gMu*(self.n+self.gKappa)/self.gBeta)*pow(1-self.Uouter(x),(self.gKappa+self.gMu*self.gLambda)/self.gBeta)
        return 0
    #Overall U,C,G vs r and t
    def C(self, r, t):
        """Calculate C(r,t) with r in cm and t in s."""
        if t < self.tc: #incoming shock
            return self.Cshock( self.xi(r,t) )
        if r < self.rs(t): #central flow
            return self.Ccentral( self.xi(r,t) )
        if r >= self.rs(t): #outer flow
            return self.Couter (self.xi(r,t) )
    def U(self, r, t):
        """Calculate U(r,t) with r in cm and t in s"""
        if t < self.tc: #incoming shock
            return self.Ushock( self.xi(r,t) )
        if r < self.rs(t): #central flow
            return self.Ucentral( self.xi(r,t) )
        if r >= self.rs(t): #outer flow
            return self.Uouter (self.xi(r,t) )
    def G(self, r, t):
        """Calculate G(r,t) with r in cm and t in s"""
        if t < self.tc: #incoming shock
            return self.Gshock( self.xi(r,t) )
        if r < self.rs(t): #central flow
            return self.Gcentral( self.xi(r,t) )
        if r >= self.rs(t): #outer flow
            return self.Gouter (self.xi(r,t) )

    # ------------------------------------
    # Hydro variables
    # ------------------------------------
    def u(self, r, t):
        """Fluid velocity u(r,t) with r in cm and t in s. Returns u in um/ns, relative to shock direction."""
        if t != self.tc:
            return self.alpha*r/math.fabs(t-self.tc)*self.U(r,t)
        return 0
    def uLab(self, r, t):
        """Fluid velocity u(r,t) with r in cm and t in s. Returns u in um/ns, in lab frame."""
        ret = self.u(r,t)
        if t != self.tc:
            if t < self.tc: #incoming shock
                return -1.*ret
            if r < self.rs(t): #central flow
                return ret
            if r >= self.rs(t): #outer flow
                return ret
        return 0
    def c(self, r, t):
        """Sound speed c(r,t) with r in cm and t in s Returns c in um/ns."""
        if t != self.tc:
            return self.alpha*r/math.fabs(t-self.tc)*self.C(r,t)
        return 0
    def rho(self, r, t):
        """Density rho(r,t) with r in cm and t in s Returns rho in g/cm3."""
        return self.rho0 * pow((r/self.r0), self.gKappa) * self.G(r,t)
    def T(self, r, t):
        """One-fluid hydro temperature T(r,t) with r in cm and t in s Returns T in keV."""
        return pow(10,-3) * self.FuelA * self.AMU * pow(self.c(r,t),2) / self.Erg2eV
    def P(self, r, t):
        """One-fluid hydro pressure P(r,t) with r in cm and t in s Returns P in GBar."""
        return (1 + self.FuelZ) * self.rho(r,t) * pow(self.c(r,t),2) * pow(10,-15)
    def Ti(self, r, t):
        """Rygg-style ion temperature with defined e-i coupling Ti(r,t) with r in cm and t in s Returns Ti in keV."""
        return (0.001 * self.FuelA * self.AMU * pow(self.c(r,t),2) / self.Erg2eV) * (1.0 - self.eiCoup * self.FuelZ / (1.+self.FuelZ))
    def ni(self, r, t):
        """Ion number density ni(r,t) with r in cm and t in s Returns ni in 1/cm^3."""
        return self.rho(r,t) / (self.FuelA * self.AMU)
    def ne(self, r, t):
        """Electron number density ne(r,t) with r in cm and t in s Returns ne in 1/cm^3."""
        return self.ni(r,t) * self.FuelZ
