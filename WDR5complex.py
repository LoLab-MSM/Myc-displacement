from pysb import *
from pysb.macros import *
from numpy import linspace
import matplotlib.pyplot as plt
from pysb.bng import run_ssa
from pysb.integrate import odesolve
import scipy.constants as sc

Model() 

# Binding domain names are taken from the summary model found in the 
# Dias, Nguyen, Georgiev et al., 2014 article

VOL = 1e-12 # L (1 picoliter is typical volume of eukaryotic cell)

Monomer('WDR5',   ['s1'])
Monomer('Myc',    ['DVDL'])
Monomer('RbBP5',  ['DVDL'])
Monomer('KANSL2', ['DVDL'])

# Dissocation constants (from experiment)
Kd_WDR5_Myc    = 2.7*1e-6*sc.N_A*VOL # uM -> molecules
Kd_WDR5_RbBP5  = 5.6*1e-6*sc.N_A*VOL # uM -> molecules
Kd_WDR5_KANSL2 = 8.6*1e-6*sc.N_A*VOL # uM -> molecules

# Forward binding rate constants in units of /s
Parameter('kf_WDR5_Myc',    1e-6)
Parameter('kf_WDR5_RbBP5',  1e-6)
Parameter('kf_WDR5_KANSL2', 1e-6)

# Reverse binding rates are in units of /s
Parameter('kr_WDR5_Myc',    kf_WDR5_Myc.value    * Kd_WDR5_Myc)
Parameter('kr_WDR5_RbBP5',  kf_WDR5_RbBP5.value  * Kd_WDR5_RbBP5)
Parameter('kr_WDR5_KANSL2', kf_WDR5_KANSL2.value * Kd_WDR5_KANSL2)

# The initial protein levels (uM -> molecules) are set to the estimates provided by Dr. Tansey. 
Parameter('WDR5_0',   round(2.657e-02*1e-6*sc.N_A*VOL)) # 16000
Parameter('Myc_0',    round(8.303e-03*1e-6*sc.N_A*VOL)) # 5000
Parameter('RbBP5_0',  round(7.473e-03*1e-6*sc.N_A*VOL)) # 4500
Parameter('KANSL2_0', round(3.321e-03*1e-6*sc.N_A*VOL)) # 2000

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


    


