# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 13:11:44 2014

@author: marion
"""
from WDR5complex import model
from numpy import linspace
import matplotlib.pyplot as plt
from pysb.bng import run_ssa
from pysb.integrate import odesolve
import numpy as np
import re

t = linspace(0,5,100)
ref = np.array([p.value for p in model.parameters])

# Run simulations
x = odesolve(model,t,verbose=True)
y = run_ssa(model,t,seed=100,verbose=True)

# Plot observables (each in a different figure)
for obs in model.observables:
    plt.figure(str(obs))
    plt.plot(t[:len(x[str(obs)])],x[str(obs)],lw=3,label=str(obs)+" (ODE)")
    plt.plot(t[:len(y[str(obs)])],y[str(obs)],lw=2,label=str(obs)+" (SSA)")
    plt.legend(loc=0)
    plt.xlim(xmax=t[-1])
    plt.xlabel("Time")
    plt.ylabel("Copy number")
    plt.ticklabel_format(useOffset=False)
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
        
plt.show()
