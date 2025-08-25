# -*- coding: utf-8 -*-
"""
Created on Mar 12, 2024

by Robert Wolle
Phys512 Midterm Problem 1c

This program calculates and plots energy bands in the Brillouin Zone for graphene
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200 # setting the dpi of outputted graphs.

t = 1
a = 1

# sets number of points between unit distances in k-space:
steps = 50

def ff_star(k_x,k_y, a):
    ff_star = 1 + 4*np.cos(np.pi*k_x - np.pi*k_y)*np.cos(-np.pi*k_x - np.pi*k_y) + 4*(np.cos(-np.pi*k_x - np.pi*k_y)**2)
    return ff_star

def E_plus(t, ff_star):
    E_plus = t*np.sqrt(ff_star)
    return E_plus


def range_k(start, end, steps):
    x_start = start[0]
    y_start = start[1]
    x_end = end[0]
    y_end = end[1]
    k_len = np.sqrt((x_end-x_start)**2 + (y_end-y_start)**2)
    num = round(steps*k_len)
    k_x = np.linspace(x_start,x_end,num=num)
    k_y = np.linspace(y_start,y_end,num=num)
    k_dist = np.linspace(0, k_len,num)
    return k_x, k_y, k_dist

def plot_contour_name(current, contour_vline, Gamma, M, K):
    name = 'unknown'
    if current[0] == Gamma[0] and current[1] == Gamma[1]:
        name = 'Gamma'
    if current[0] == M[0] and current[1] == M[1]:
        name = 'M'
    if current[0] == K[0] and current[1] == K[1]:
        name = 'K'    
    offsetx = -0.06
    offsety = -0.5
    plt.text(contour_vline+offsetx,offsety,name,rotation=90)


# Contour points:
Gamma = np.array([0,0])
M = np.array([1/2,0])
K = np.array([1/3,1/3])

# Contour path:
contour = [Gamma, M, K, Gamma]


k_dist_prev = 0 # first plotted contour
contour_vlines = np.zeros(len(contour))
plt.axvline(contour_vlines[0],color='k')
plot_contour_name(contour[0], contour_vlines[0], Gamma, M, K)
for c in range(1,len(contour)):
    start_point = contour[c-1]
    end_point = contour[c]
    k_x, k_y, k_dist = range_k(start_point, end_point, steps)
    E = np.zeros([len(k_dist)])
    for i in range(len(k_dist)):
        E[i] = E_plus(t, ff_star(k_x[i], k_y[i], a))
    
    plt.plot(k_dist+k_dist_prev,E,'b')
    plt.plot(k_dist+k_dist_prev,-E,'r')
    
    
    contour_vlines[c] = k_dist[len(k_dist)-1] + contour_vlines[c-1]
    plot_contour_name(contour[c], contour_vlines[c], Gamma, M, K)
    plt.axvline(contour_vlines[c],color='k')
    
    k_dist_prev = k_dist[len(k_dist)-1]+k_dist_prev

#plt.gca().set_xticks([])

plt.show()
