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


colors = ('r','k','b','g','0.5')
points = linspace(0.01,1,100) 
ref = np.array([p.value for p in model.parameters])

for i in points:  
    for j in range(len(model.parameters)):
        if re.search('_0$', model.parameters[j].name):
            model.parameters[j].value *= i
    x=run_ssa(model,t,verbose=False)
    for j in range(len(model.observables)):
        obs = model.observables[j]
        plt.figure(str(obs))
        plt.legend(loc=0)
        if i == 0.01:
            plt.plot(t,x[str(obs)]/np.max(x[str(obs)]),color=colors[j], label=str(obs))
        else:    
            plt.plot(t,x[str(obs)]/np.max(x[str(obs)]),color=colors[j])
    for j in range(len(model.parameters)):
        model.parameters[j].value = ref[j]
n_conc = np.zeros((len(model.observables),np.shape(points)[0]))
counter = 0
for i in points:
    for j in range(len(model.parameters)):
        if re.search('_0$', model.parameters[j].name):
            model.parameters[j].value *= i
    x=odesolve(model,t,verbose=False)
    for j in range(len(model.observables)):
        obs = model.observables[j]
        plt.figure(str(obs))
        plt.legend(loc=0)
        plt.title(str(model.observables[j].name))
        n_conc[j,counter] = x[str(obs)][-1]#/np.max(x[str(obs)])
        if i == 0.01:
            plt.plot(t,x[str(obs)]/np.max(x[str(obs)]),color=colors[j], label=str(obs))
            print x[str(obs)]
        else:    
            plt.plot(t,x[str(obs)]/np.max(x[str(obs)]),color=colors[j])
    for j in range(len(model.parameters)):
        model.parameters[j].value = ref[j]
    counter +=1


for j in range(len(model.observables)):
    plt.title(str(model.observables[j].name))
    for i in range(len(points)):
        plt.plot(points[i], n_conc[j,i], '.',color=colors[j])
    plt.show()