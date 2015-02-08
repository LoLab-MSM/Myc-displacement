# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 13:11:44 2014

@author: marion
"""
from WDR5complex import model
from WDR5complex import set_volume
from WDR5complex import set_initial_conditions
from numpy import linspace
import matplotlib.pyplot as plt
from pysb.bng import run_ssa
from pysb.integrate import odesolve
import numpy as np
import re

run_ode = True

# set_volume(10.*model.parameters['VOL'].value)
# set_initial_conditions(myc=100.*model.parameters['Myc_0'].value)

t = linspace(0,5,101)
# ref = np.array([p.value for p in model.parameters])

# Run simulations
x = run_ssa(model,t,seed=100,verbose=True)
if run_ode:
    y = odesolve(model,t,verbose=True)

print
for ic in model.initial_conditions:
    print "%s: %g" % (ic[1].name, ic[1].value)
print
for p in model.parameters_rules():
    print "%s: %g" % (p.name, p.value)
print
for e in model.expressions:
    subs = [(p.name, p.value) for p in model.parameters]
    print "%s: %s " % (e.name, e.expand_expr(model).subs(subs))

# Plot observables (each in a different figure)

#Free species
fig_free, axs_free = plt.subplots(4, 1, sharex=True)
fig_free.canvas.set_window_title('Free species')
n=0

# Complexes
fig_cpx, axs_cpx = plt.subplots(3, 1, sharex=True)
fig_cpx.canvas.set_window_title('Complexes')
m=0

def change_yticks(mult, old_ticks):
    new_ticks = [old_ticks[0]]
    while new_ticks[-1] < old_ticks[-1]:
        new_ticks.append( new_ticks[-1] + mult*(old_ticks[1]-old_ticks[0]) )
    return new_ticks

def plot_fig(axs):
    axs.plot(t[:len(x[str(obs)])], x[str(obs)], lw=3, label=str(obs)+" (SSA)", color='blue')
    if run_ode:
        axs.plot(t[:len(y[str(obs)])], y[str(obs)], lw=2, label=str(obs)+" (ODE)", color='red')
    axs.legend(loc=0)
    axs.set_xlim(xmax=t[-1])
    axs.set_xlabel("Time")
    axs.set_ylabel("Copy number")
    axs.ticklabel_format(useOffset=False)
    axs.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    # Double the distance between the y-ticks
    axs.yaxis.set_ticks(change_yticks(2., axs.yaxis.get_ticklocs()))

for obs in model.observables:
    if re.search(r'_free$', str(obs)):
        plot_fig(axs_free[n])
        n += 1
    else:
        plot_fig(axs_cpx[m])
        m += 1

plt.show()
