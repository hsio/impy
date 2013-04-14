# Calculate burn-averaged plasma parameters
# A. Zylstra 2013/04/13

from Implosion import *
from Resources.IO import *
from Resources.Plasma import *
from Resources import Fusion
from numpy import arange
from numpy import zeros
import math
import csv
import os

# global implosion
impl = 0

# reactions to weight
reactions = []
# syntax: [ name , A1, Z1, A2, Z2, reactivity fn ]
reactions.append( [ "DD" , 2, 1, 2, 1, Fusion.DD ] )
reactions.append( [ "D3He" , 2, 1, 3, 2, Fusion.D3He ] )
reactions.append( [ "Vol" , 0, 0, 0, 0, lambda x: 1 ]) # volume weighted
reactions.append( [ "Mass" , 0, 0, 0, 0, lambda x: 1 ]) # volume weighted
    
# ------------------------------------
# Burn rate calculators
# ------------------------------------
def rate(it,WeightRxn):
    """Calculate the 'rate' at time index it. Returns [lD , ND, PlasmaF , MFP , tau_i , R_cs , R_vt , tau_ei, Te]"""
    # return values
    Yield = lD = ND = PlasmaF = MFP = tau_i = R_cs = R_vt = tau_ei = Terate = Tirate = nerate = nirate = 0

    # pull data from reaction:
    A1 = WeightRxn[1]
    Z1 = WeightRxn[2]
    A2 = WeightRxn[3]
    Z2 = WeightRxn[4]
    
    # iterate over radius
    for ir in range( impl.ir_min() , impl.ir_fuel() ):
        vol = impl.vol(ir,it)
        f1 = impl.f(ir,it,A1,Z1)
        f2 = impl.f(ir,it,A2,Z2)

        # for fusion reaction weighting:
        dblcount = 1
        if (A1 == A2) and (Z1 == Z2):
            dblcount = 2 #account for factor of 2 if reactants are identical

        # plasma conditions:
        ni = impl.ni(ir,it)
        Ti = impl.Ti(ir,it)
        ne = impl.ne(ir,it)
        Te = impl.Te(ir,it)
        A = impl.Abar(ir,it)
        Z = impl.Zbar(ir,it)
        
        # weight function
        w = WeightRxn[5](Ti)*pow(ni,2)*(f1*f2/dblcount)*vol
        # special case for volume weighted:
        if WeightRxn[0] == "Vol":
            w = vol
        # special case for mass weighted:
        if WeightRxn[0] == "Mass":
            w = ni*vol
        
        # add weighted incremental values
        Yield += w
        lD += w*LambdaD(ne,Te)
        ND += w*PlasmaParameter(ne,Te)
        PlasmaF += w*fpe(ne)
        MFP += w*IonMFP(ni,Ti,Z,A)
        tau_i += w*Taui(ni,Ti,Z,A)
        rfuel = impl.r( impl.ir_fuel() , it )
        R_cs += w * rfuel / impl.c(ir,it)
        R_vt += w * rfuel / uTherm(Ti,A)
        tau_ei += w * Tauei(ni,Ti,Z)
        Terate += w * Te
        Tirate += w * Ti
        nerate += w * ne
        nirate += w * ni
    
    # return array
    return [Yield, lD , ND, PlasmaF , MFP , tau_i , R_cs , R_vt , tau_ei , Terate , Tirate , nerate , nirate]
        
 
# ------------------------------------
# Main method
# ------------------------------------
def run(i):
    """Calculate total yield."""
    # input sanity check:
    if not isinstance(i,Implosion):
        print("WARNING: invalid input.")
        return
    global impl
    impl = i # global implosion variable for this module
        
    # things to calculate
    Yield = 0
    lD = 0 # debye length in cm
    ND = 0 # plasma parameter
    PlasmaF = 0 # plasma frequency
    MFP = 0 # ion MFP
    tau_i = 0 # ion-ion collision time
    R_vt = 0 # ion transit time R/vt
    R_cs = 0 # radius / sound speed
    tau_ei = 0 # ion-electron collision time
    Te = 0 # electron temp
    Ti = 0 # ion temp
    ne = 0 # electron density
    ni = 0 # ion denisty
    
    # output files
    file = csv.writer(open(os.path.join(OutputDir,'BurnAvgPlasma.csv'),'w'))
    file.writerow( ['Weight','Quantity','Value','Units'] )
    file_TR = csv.writer(open(os.path.join(OutputDir,'BurnAvgPlasma_TR.csv'),'w'))
    file_TR.writerow( ['Weight','Time','Debye Length','Plasma Parameter','Plasma Freq.','MFP','Tau_i','R/vt','R/cs','Tau_ei','Te','Ti','ne','ni'] )

    # iterate over each weighting function:
    for WeightRxn in reactions:
        Yield = lD = ND = PlasmaF = MFP = tau_i = R_vt = R_cs = tau_ei = Te = Ti = ne = ni = 0
        #iterate over all time:
        for it in range(impl.it_min(), impl.it_max()):
            dRate = rate(it,WeightRxn)
            Yield += dRate[0]
            lD += dRate[1]
            ND += dRate[2]
            PlasmaF += dRate[3]
            MFP += dRate[4]
            tau_i += dRate[5]
            R_vt += dRate[6]
            R_cs += dRate[7]
            tau_ei += dRate[8]
            Te += dRate[9]
            Ti += dRate[10]
            ne += dRate[11]
            ni += dRate[12]

            # get absolute time:
            t = impl.t(it)

            # write to time-resolved file:
            file_TR.writerow( [ WeightRxn[0], t, 
                dRate[1]/dRate[0],
                dRate[2]/dRate[0],
                dRate[3]/dRate[0],
                dRate[4]/dRate[0],
                dRate[5]/dRate[0],
                dRate[6]/dRate[0],
                dRate[7]/dRate[0],
                dRate[8]/dRate[0],
                dRate[9]/dRate[0],
                dRate[10]/dRate[0],
                dRate[11]/dRate[0],
                dRate[12]/dRate[0] ] )
            
        #normalization
        lD = lD/Yield
        ND = ND/Yield
        PlasmaF = PlasmaF / Yield
        MFP = MFP / Yield
        tau_i = tau_i / Yield
        R_vt = R_vt / Yield
        R_cs = R_cs / Yield
        tau_ei = tau_ei / Yield
        Te = Te / Yield
        Ti = Ti / Yield
        ne = ne / Yield
        ni = ni / Yield
        
        file.writerow([ WeightRxn[0] ,  " Debye length",lD,"cm"])
        file.writerow([ WeightRxn[0] , "Plasma parameter",ND])
        file.writerow([ WeightRxn[0] ,  "Plasma freq." , PlasmaF , "Hz" ])
        file.writerow([ WeightRxn[0] ,  "Mean Free Path" , MFP , "cm" ])
        file.writerow([ WeightRxn[0] ,  "Ion collision time" , tau_i , "s" ])
        file.writerow([ WeightRxn[0] ,  "R/vt" , R_vt , "s" ])
        file.writerow([ WeightRxn[0] ,  "R/cs" , R_cs , "s" ])
        file.writerow([ WeightRxn[0] ,  "tau_ei" , tau_ei , "s" ])
        file.writerow([ WeightRxn[0] ,  "Te" , Te , "keV" ])
        file.writerow([ WeightRxn[0] ,  "Ti" , Ti , "keV" ])
        file.writerow([ WeightRxn[0] ,  "ne" , ne , "1/cc" ])
        file.writerow([ WeightRxn[0] ,  "ni" , ni , "1/cc" ])

