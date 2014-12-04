from pysb import *
from pysb.macros import *
from numpy import linspace
import matplotlib.pyplot as plt
from pysb.bng import run_ssa
from pysb.integrate import odesolve

Model() 

#Binding domain names are taken from the summary model found in the 
#Dias, Nguyen, Georgiev et al., 2014 article

#Individual protein/protein complexes that are known to associate with WDR5 as 
#a trimeric complex, including WDR5
#Abbreviations: 'b'= bound, 'u'= unbound

Monomer('WDR5',['s1'])
#or Monomer('WDR5',['s1','s2'])
Monomer('Myc',['DVDL'])
Monomer('RbBP5',['DVDL'])
Monomer('KANSL2',['DVDL'])
#Monomer('MLL',['WIN'])
#Monomer('KANSL1',['WIN'])

kd_WDR5_Myc = 2.7e-2
kd_WDR5_RbBP5 = 5.6e-2
kd_WDR5_KANSL2 = 8.6e-2
kdimf = 5.0e-2
kdimr_WDR5_Myc = kdimf/kd_WDR5_Myc
kdimr_WDR5_RbBP5 = kdimf/kd_WDR5_RbBP5
kdimr_WDR5_KANSL2 = kdimf/kd_WDR5_KANSL2
print kdimr_WDR5_Myc, kdimr_WDR5_RbBP5, kdimr_WDR5_KANSL2

#Forward binding rates are set to 50 uM^-2 * s^-1.
Parameter('kf_WDR5_Myc',kdimf)
Parameter('kf_WDR5_RbBP5',kdimf)
Parameter('kf_WDR5_KANSL2',kdimf)
#Parameter('kf_WDR5_MLL',kdimf)
#Parameter('kf_WDR5_KANSL1',kdimf)
#?Parameter('kf_RbBP5_WDR5_Myc',kdimf)
#?Parameter('kf_RbBP5_WDR5_MLL',kdimf)
#?Parameter('kf_KANSL2_WDR5_MLL',kdimf)
#?Parameter('kf_KANSL2_WDR5_KANSL1',kdimf)
#?Parameter('kf_RbBP5_WDR5_KANSL1',kdimf)
#?Parameter('kf_KANSL2_WDR5_Myc',kdimf)
#?Parameter('kf_KANSL2_WDR5_MLL',kdimf)
#?Parameter('kf_MLL_WDR5_RbBP5',kdimf)
#?Parameter('kf_Myc_WDR5_RbBP5',kdimf)

#Reverse binding rates are in units uM^-1 * s^-1
Parameter('kr_WDR5_Myc',kdimf/kd_WDR5_Myc)
Parameter('kr_WDR5_RbBP5',kdimf/kd_WDR5_RbBP5)
Parameter('kr_WDR5_KANSL2',kdimf/kd_WDR5_KANSL2)

#Reverse binding rates that have yet to be determined
#Parameter('kr_WDR5_MLL',kdimf/kd_WDR5_MLL)
#Parameter('kr_WDR5_KANSL1',kdimf/kd_WDR5_KANSL1)
#?Parameter('kr_RbBP5_WDR5_Myc',kdimf/kd_RbBP5_WDR5_Myc)
#?Parameter('kr_RbBP5_WDR5_MLL',kdimf/kd_RbBP5_WDR5_MLL)
#?Parameter('kr_KANSL2_WDR5_MLL',kdimf/kd_KANSL2_WDR5_MLL)
#?Parameter('kr_KANSL2_WDR5_KANSL1',kdimf/kd_KANSL2_WDR5_KANSL1)
#?Parameter('kr_RbBP5_WDR5_KANSL1',kdimf/kd_RbBP5_WDR5_KANSL12)
#?Parameter('kr_KANSL2_WDR5_Myc',kdimf/kd_KANSL2_WDR5_Myc)
#?Parameter('kr_KANSL2_WDR5_MLL',kdimf/kd_KANSL2_WDR5_MLL)
#?Parameter('kr_MLL_WDR5_RbBP5',kdimf/kd_MLL_WDR5_RbBP5)
#?Parameter('kr_Myc_WDR5_RbBP5',kdimf/kd_Myc_WDR5_RbBP5)

#The parameters are set to the molecular estimates provided by Dr. Tansey. 
#Will need to determine quantity of the number of molecules by 
#experimentation.
Parameter('WDR5_0',16000)
Parameter('Myc_0',5000)
Parameter('RbBP5_0',4500)
Parameter('KANSL2_0',2000)
#Parameter('MLL_0',2e3)
#Parameter('KANSL1_0',2e3)


#All monomers are initially described in the unbound state at all possible
#binding sites.
Initial(WDR5(s1=None), WDR5_0)
Initial(Myc(DVDL=None), Myc_0)
Initial(RbBP5(DVDL=None), RbBP5_0)
Initial(KANSL2(DVDL=None), KANSL2_0)
#Initial(MLL(WIN=None), MLL_0)
#Initial(KANSL1(WIN=None), KANSL1_0)


#Rules for Myc binding WDR5
Rule('WDR5_binds_Myc', WDR5(s1=None) + Myc(DVDL=None) <> WDR5(s1=1) % Myc(DVDL=1), kf_WDR5_Myc, kr_WDR5_Myc)
#Rule('WDR5_binds_Myc', WDR5(s1=None,s2=None) + Myc(WIN=None) <> WDR5(s1=1,s2=None) % Myc(WIN=1), kf_WDR5_Myc, kr_WDR5_Myc)
#Rule('WDR5_RbBP5_binds_Myc', WDR5(s1=None,s2=1) % RbBP5(DVDL=1) + Myc(WIN=None) <> WDR5(s1=2,s2=1) % RbBP5(DVDL=1) % Myc(WIN=2), kf_RbBP5_WDR5_Myc, kr_RbBP5_WDR5_Myc)
#Rule('WDR5_KANSL2_binds_Myc', WDR5(s1=None,s2=1) % KANSL2(DVDL=1) + Myc(WIN=None) <> WDR5(s1=2,s2=1) % KANSL2(DVDL=1) % Myc(WIN=2), kf_KANSL2_WDR5_Myc, kr_KANSL2_WDR5_Myc)

#Rules for RbBP5 binding WDR5
Rule('WDR5_binds_RbBP5', WDR5(s1=None) + RbBP5(DVDL=None) <> WDR5(s1=1) % RbBP5(DVDL=1), kf_WDR5_RbBP5, kr_WDR5_RbBP5)
#Rule('WDR5_binds_RbBP5', WDR5(s1=None,s2=None) + RbBP5(DVDL=None) <> WDR5(s1=None,s2=1) % RbBP5(DVDL=1), kf_WDR5_RbBP5, kr_WDR5_RbBP5)
#Rule('WDR5_MLL_binds_RbBP5', WDR5(s1=1,s2=None) % MLL(WIN=1) + RbBP5(DVDL=None) <> WDR5(s1=1,s2=2) % MLL(WIN=1) + RbBP5(DVDL=2), kf_MLL_WDR5_RbBP5, kr_MLL_WDR5_RbBP5)
#Rule('WDR5_Myc_binds_RbBP5', WDR5(s1=1,s2=None) % Myc(WIN=1) + RbBP5(DVDL=None) <> WDR5(s1=1,s2=2) % Myc(WIN=1) + RbBP5(DVDL=2), kf_Myc_WDR5_RbBP5, kr_Myc_WDR5_RbBP5)

#Rules for KANSL2 binding WDR5
Rule('WDR5_binds_KANSL2', WDR5(s1=None) + KANSL2(DVDL=None) <> WDR5(s1=1) % KANSL2(DVDL=1), kf_WDR5_KANSL2, kr_WDR5_KANSL2)
#Rule('WDR5_binds_KANSL2', WDR5(s1=None,s2=None) + KANSL2(DVDL=None) <> WDR5(s1=None,s2=1) % RbBP5(DVDL=1), kf_WDR5_KANSL2, kr_WDR5_KANSL2)
#Rule('WDR5_MLL_binds_KANSL2', WDR5(s1=1,s2=None) % MLL(WIN=1) + KANSL2(DVDL=None) <> WDR5(s1=1,s2=2) % MLL(WIN=1) % KANSL2(DVDL=2), kf_MLL_WDR5_KANSL2, kr_MLL_WDR5_KANSL2)
#Rule('WDR5_Myc_binds_KANSL2', WDR5(s1=1,s2=None) % Myc(WIN=1) + KANSL2(DVDL=None) <> WDR5(s1=1,s2=2) % Myc(WIN=1) % KANSL2(DVDL=2), kf_Myc_WDR5_KANSL2, kr_Myc_WDR5_KANSL2)
#Rule('WDR5_KANSL1_binds_KANSL2', WDR5(s1=1,s2=None) % KANSL1(WIN=1) + KANSL2(DVDL=None) <> WDR5(s1=1,s2=2) % KANSL1(WIN=1) % KANSL2(DVDL=2), kf_KANSL1_WDR5_KANSL2, kr_KANSL1_WDR5_KANSL2)

#Rules for MLL binding WDR5
#Rule('WDR5_binds_MLL', WDR5(s1=None) + MLL(WIN=None) <> WDR5(s1=1) % MLL(WIN=1), kf_WDR5_MLL, kr_WDR5_MLL)
#Rule('WDR5_binds_MLL', WDR5(s1=None,s2=None) + MLL(WIN=None) <> WDR5(s1=1,s2=None) % Myc(WIN=1), kf_WDR5_MLL, kr_WDR5_MLL)
#Rule('WDR5_RbBP5_binds_MLL', WDR5(s1=None,s2=1) % RbBP5(DVDL=1) + MLL(WIN=None) <> WDR5(s1=2,s2=1) % RbBP5(DVDL=1) % MLL(WIN=2), kf_RbBP5_WDR5_MLL, kr_RbBP5_WDR5_MLL)
#Rule('WDR5_KANSL2_binds_MLL', WDR5(s1=None,s2=1) % KANSL2(DVDL=1) + MLL(WIN=None) <> WDR5(s1=2,s2=1) % KANSL2(DVDL=1) % MLL(WIN=2), kf_KANSL2_WDR5_MLL, kr_KANSL2_WDR5_MLL)

#Rules for KANSL1 binding WDR5
#Rule('WDR5_binds_KANSL1', WDR5(s1=None) + KANSL1(WIN=None) <> WDR5(s1=1) % KANSL1(WIN=1), kf_WDR5_KANSL1, kr_WDR5_KANSL1)
#Rule('WDR5_binds_KANSL1', WDR5(s1=None,s2=None) + KANSL1(WIN=None) <> WDR5(s1=1,s2=None) % KANSL1(WIN=1), kf_WDR5_KANSL1, kr_WDR5_KANSL1)
#Rule('WDR5_RbBP5_binds_KANSL1', WDR5(s1=None,s2=1) % RbBP5(DVDL=1) + KANSL1(WIN=None) <> WDR5(s1=2,s2=1) % RbBP5(DVDL=1) % KANSL1(WIN=2), kf_RbBP5_WDR5_KANSL1, kr_RbBP5_WDR5_KANSL1)
#Rule('WDR5_KANSL2_binds_KANSL1', WDR5(s1=None,s2=1) % KANSL2(DVDL=1) + KANSL1(WIN=None) <> WDR5(s1=2,s2=1) % KANSL2(DVDL=1) % KANSL1(WIN=2), kf_KANSL2_WDR5_KANSL1, kr_KANSL2_WDR5_KANSL1)

#Observable species
Observable('WDR5_free', WDR5(s1=None))
Observable('WDR5_Myc', WDR5(s1=1) % Myc(DVDL=1))
Observable('WDR5_RbBP5', WDR5(s1=1) % RbBP5(DVDL=1))
Observable('WDR5_KANSL2', WDR5(s1=1) % KANSL2(DVDL=1))
Observable('Myc_free', Myc(DVDL=None))

    


