from pysb import *
from pysb.macros import *
from numpy import linspace
import matplotlib.pyplot as plt
from pysb.bng import run_ssa
from pysb.integrate import odesolve
from scipy.constants import N_A
from sympy import sympify

Model() 

# Parameter('test1', 1.)
# Parameter('test2', 1.)
# Parameter('test3', 1.)
# Parameter('test4', 1.)
# Expression('TEST1', sympify("test1*test2"))
# Expression('TEST2', sympify("TEST1*test3"))
# Expression('TEST3', sympify("TEST2*test4"))

# Cell volume (1 picoliter is typical for eukaryotes)
Parameter('VOL', 1e-12) # L 

# Initial protein copy numbers (estimates from Prof. Tansey)
Parameter('WDR5_0',   16000) # molec
Parameter('Myc_0',    5000)  # molec
Parameter('RbBP5_0',  4500)  # molec
Parameter('KANSL2_0', 2000)  # molec

def set_volume(vol):
    WDR5_0.value   *= (vol/VOL.value)
    Myc_0.value    *= (vol/VOL.value)
    RbBP5_0.value  *= (vol/VOL.value)
    KANSL2_0.value *= (vol/VOL.value)
    VOL.value       = vol

def set_initial_conditions(wdr5=None, myc=None, rbbp5=None, kansl2=None, conc=False):
    for (X,y) in [(WDR5_0, wdr5), (Myc_0, myc), (RbBP5_0, rbbp5), (KANSL2_0, kansl2)]:
        if (y != None):
            if conc: 
                y *= N_A * VOL.value # M -> molec
            X.value   = y   # molec
    
# Dissocation constants (from Prof. Tansey)
Parameter('Kd_WDR5_Myc',    2.7*1e-6) # uM -> M
Parameter('Kd_WDR5_RbBP5',  5.6*1e-6) # uM -> M
Parameter('Kd_WDR5_KANSL2', 8.6*1e-6) # uM -> M

# Reverse binding rate constants
Parameter('kr_WDR5_Myc',    100) # /time
Parameter('kr_WDR5_RbBP5',  100) # /time
Parameter('kr_WDR5_KANSL2', 100) # /time

# Forward binding rate constants
Expression('kf_WDR5_Myc',    sympify("kr_WDR5_Myc    / Kd_WDR5_Myc    / %g / VOL" % N_A)) # /M-time -> /time
Expression('kf_WDR5_RbBP5',  sympify("kr_WDR5_RbBP5  / Kd_WDR5_RbBP5  / %g / VOL" % N_A)) # /M-time -> /time
Expression('kf_WDR5_KANSL2', sympify("kr_WDR5_KANSL2 / Kd_WDR5_KANSL2 / %g / VOL" % N_A)) # /M-time -> /time

# Binding domain names are taken from the summary model found in the 
# Dias, Nguyen, Georgiev et al., 2014 article

# Monomers
Monomer('WDR5',   ['s1'])
Monomer('Myc',    ['DVDL'])
Monomer('RbBP5',  ['DVDL'])
Monomer('KANSL2', ['DVDL'])

# Initial conditions
Initial(WDR5(s1=None),     WDR5_0)
Initial(Myc(DVDL=None),    Myc_0)
Initial(RbBP5(DVDL=None),  RbBP5_0)
Initial(KANSL2(DVDL=None), KANSL2_0)

# Reversible binding rules
Rule('WDR5_binds_Myc',    WDR5(s1=None) + Myc(DVDL=None)    <> WDR5(s1=1) % Myc(DVDL=1),    kf_WDR5_Myc,    kr_WDR5_Myc)
Rule('WDR5_binds_RbBP5',  WDR5(s1=None) + RbBP5(DVDL=None)  <> WDR5(s1=1) % RbBP5(DVDL=1),  kf_WDR5_RbBP5,  kr_WDR5_RbBP5)
Rule('WDR5_binds_KANSL2', WDR5(s1=None) + KANSL2(DVDL=None) <> WDR5(s1=1) % KANSL2(DVDL=1), kf_WDR5_KANSL2, kr_WDR5_KANSL2)

# Observables

# Free species
Observable('WDR5_free',   WDR5(s1=None))
Observable('Myc_free',    Myc(DVDL=None))
Observable('RbBP5_free',  RbBP5(DVDL=None))
Observable('KANSL2_free', KANSL2(DVDL=None))

# Complexes
Observable('WDR5_Myc',    WDR5(s1=1) % Myc(DVDL=1))
Observable('WDR5_RbBP5',  WDR5(s1=1) % RbBP5(DVDL=1))
Observable('WDR5_KANSL2', WDR5(s1=1) % KANSL2(DVDL=1))
   
######################################################################################################
# From: "Tansey, William P." <william.p.tansey@Vanderbilt.Edu>
# To: "Lopez, Carlos F" <c.lopez@Vanderbilt.Edu>
# Subject: Myc
# Date: October 22, 2014 at 9:16:35 AM CDT
#
# Hi Carlos,
# 
# Enjoyed our chat yesterday, as always. I think that we will be able to generate some juicy data 
# for modeling as we develop better assays for MYC transcription etc. moving forward. In the meantime, 
# here's a recap of the more immediate challenge we discussed yesterday. As yet, none of this is 
# published so please keep to yourself/group.
# 
# In a nutshell, we've found that MYC needs to interact with WDR5 to bind chromatin and work as a 
# transcription factor. WDR5 is present in MLL-type histone methyltransferases, where it binds to 
# many proteins, including RBBP5. WDR5 is also present in the NSL histone acetyltransferase complex, 
# where it binds many proteins, including KANSL2. The RBBP5-WDR5 and KANSL2-WDR5 interfaces are 
# virtually identical to each other, and this has been argued as an important mechanism that directs 
# WDR5 into one complex and not the other.
# 
# As mentioned yesterday, MYC binds WDR5 in exactly the same way as does RBBP5 and KANSL2 (picture 
# below and more detail in attached papers). 
# 
# Our model is that MYC displaces RBBP5 and/or KANSL2 from the MLL and/or NSL complexes by direct 
# competition for binding. What would be fun to determine, is whether this is possible, and what the 
# relative levels of MYC, WDR5, RBBP5, and KANSL2 would need to be to achieve a certain % of replacement.
# 
# The relevant Kds are
# 
# MYC:WDR5 2.7 uM
# RBBP5:WDR5 5.6 uM
# KANSL2:WDR5 8.6 uM.
# 
# Cellular concentrations are harder to pin down. For MYC, we can say that normal levels are 5,000 molecules 
# per cell. Typical overexpression is 75,000 molecules per cell. But the numbers can go as high as 500,000 
# molecules per cell, depending on experimental setup. For WDR5, the attached paper gives 16,000 molecules 
# per U2OS cell (I think that's on the low side), 4,500 molecules per cell for RBBP5, and they don't see KANSL2. 
# Bottom line is that these numbers probably all need to be measured experimentally in our system to be able to 
# get something we can believe.
# 
# Let me know if you need anything else!
# 
# All the best,
# 
# Bill
######################################################################################################


