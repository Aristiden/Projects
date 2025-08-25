# -*- coding: utf-8 -*-
"""
Created on Mar 12, 2024

by Robert Wolle
Phys512 Midterm Problem 1c

This program plots contours of equal energy in k-space for input parameters of the tight binding approximation. 
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200 # setting the dpi of outputted graphs.

t = 1
a = 1
E_0 = 0
z = [-t, 0, t]

steps = 100



def Fermi_energy(k_x,k_y, t, a, E_0):
    Fermi_energy = E_0 + 2*t*(np.cos(2*np.pi*k_x) + np.cos(2*np.pi*k_y))
    return Fermi_energy


k_x = np.linspace(0,1,steps)
k_y = np.linspace(0,1,steps)

    
E = np.zeros([steps,steps])
for i in range(steps):
    for j in range(steps):
        E[i][j] = Fermi_energy(k_x[i],k_y[j], t, a, E_0)



contours = 5

plt.contour(k_x,k_y,E, levels=z)
plt.colorbar()
plt.show()

"""
def Fermi_energy(k_x,k_y, t, a, E_0):
    Fermi_energy = E_0 + 2*t*(np.cos(2*np.pi*k_x) + np.cos(2*np.pi*k_y))
    return Fermi_energy

k_x = np.linspace(0,1,steps)
k_y = np.linspace(0,1,steps)
for t in t_list:
    
    E = np.zeros([steps,steps])
    for i in range(steps):
        for j in range(steps):
            E[i][j] = Fermi_energy(k_x[i],k_y[j], t, a, E_0)
    
    
    
    contours = 5
    
    plt.contourf(k_x,k_y,E, levels=contours)
    plt.colorbar()
    plt.show()

"""
