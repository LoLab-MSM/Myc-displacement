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
t=linspace(0,1,50)
plt.figure()

#x=odesolve(model,t,verbose=True)
#for obs in model.observables:
#    plt.plot(t,x[str(obs)], label=str(obs) + '_ode')


#plt.savefig('WDR5_equalMycRbBP5KANSL2.png')
#plt.show()
colors = ('r','k','b','g','0.5')
for i in range(100):
    print i
    x=run_ssa(model,t,verbose=False)
    for j in range(len(model.observables)):
        obs = model.observables[j]
        if i == 0:
            plt.plot(t,x[str(obs)],color=colors[j], label=str(obs))
        else:    
            plt.plot(t,x[str(obs)],color=colors[j])
   
plt.legend(loc=0)
#plt.savefig('WDR5_equalMycRbBP5KANSL2.png')
 
ref = np.array([p.value for p in model.parameters])
for i in linspace(0.01,1,100):
    print i
    for j in range(len(model.parameters)):
        if re.search('_0$', model.parameters[j].name):
            model.parameters[j].value *= i
    x=run_ssa(model,t,verbose=False)
    for j in range(len(model.observables)):
        obs = model.observables[j]
        plt.figure(str(obs))
        #plt.yscale('log')
        plt.legend(loc=0)
        if i == 0.1:
            plt.plot(t,x[str(obs)]/i,color=colors[j], label=str(obs))
            print x[str(obs)]
        else:    
            plt.plot(t,x[str(obs)]/i,color=colors[j])
    for j in range(len(model.parameters)):
        model.parameters[j].value = ref[j]


    
plt.show()   