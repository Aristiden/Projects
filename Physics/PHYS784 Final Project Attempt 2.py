# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 15:44:08 2025

Robert Wolle

PHYS784 Final Project


"""

import time
import math
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200 # setting the dpi of outputted graphs.
from numba import njit
from scipy import sparse


#Initialize time:
start = time.time()

"""
-----------------------------------------------------------
                        FUNCTIONS
-----------------------------------------------------------

    atom_potential(atom_r, test_r, atom_q, epsilon_0, charge)
      - Takes the position and charge of an atom and calculates the
          electric potential due to that atom at a test position
          
    calculate_V(test_r, lattice_locations_list, epsilon_0, charge, total_charge)
      - Takes in a list of atom locations and charges, 
    

"""

@njit
def atom_potential(atom_r, test_r, atom_q, epsilon_0, charge):
    dist = np.sqrt((atom_r[0]-test_r[0])**2 + (atom_r[1]-test_r[1])**2 + (atom_r[2]-test_r[2])**2)
    V = 1/(4*np.pi*epsilon_0)*charge*atom_q/dist # eVm/e^2 * e/m = eV/e
    return V

@njit
def calculate_V(test_r, lattice_locations_list, epsilon_0, charge, total_charge):
    V = 0
    
    for atom in lattice_locations_list:
        V += atom_potential(atom[0:3],test_r, atom[3], epsilon_0, charge)
        total_charge += atom[3]

    return V,total_charge

def make_2D_layer(layer_origin, lattice_constants, charge, settings):
    
    """
    2D layer is in XY plane, Z is fixed.
    """
    
    lattice_locations = []
    
    x_points = math.floor((settings["x_max"]-layer_origin[0]-1)/lattice_constants[0])
    y_points = math.floor((settings["y_max"]-layer_origin[1]-1)/lattice_constants[1])
    for x_point in range(x_points+1):
        for y_point in range(y_points+1):
            i = layer_origin[0] + x_point*lattice_constants[0]
            j = layer_origin[1] + y_point*lattice_constants[1]
            k = layer_origin[2]
            lattice_locations.append([i,j,k,charge])
    return lattice_locations

def make_TMD_ML(lattice_locations_list, origin_offset, a, theta, Q_TM,Q_D, settings):
    lattice_constants = [a*2*np.cos(theta),a]
    
    #Mo XY_1:
    charge = Q_TM
    layer_origin = np.array([0,0,0])+origin_offset
    lattice_locations = make_2D_layer(layer_origin, lattice_constants, charge, settings)
    lattice_locations_list.extend(lattice_locations)

    #Mo XY_2:
    charge = Q_TM
    layer_origin = np.array([a*np.cos(theta),a*np.sin(theta),0])+origin_offset
    lattice_locations = make_2D_layer(layer_origin, lattice_constants, charge, settings)
    lattice_locations_list.extend(lattice_locations)

    #Se_bottom_1:
    charge = Q_D
    layer_origin = np.array([lattice_constants[0]/3,0,-a/2])+origin_offset
    lattice_locations = make_2D_layer(layer_origin, lattice_constants, charge, settings)
    lattice_locations_list.extend(lattice_locations)

    #Se_bottom_2:
    charge = Q_D
    layer_origin = np.array([lattice_constants[0]/3+a*np.cos(theta),a/2,-a/2])+origin_offset
    lattice_locations = make_2D_layer(layer_origin, lattice_constants, charge, settings)
    lattice_locations_list.extend(lattice_locations)

    #Se_top_1:
    charge = Q_D
    layer_origin = np.array([lattice_constants[0]/3,0,a/2])+origin_offset
    lattice_locations = make_2D_layer(layer_origin, lattice_constants, charge, settings)
    lattice_locations_list.extend(lattice_locations)

    #Se_top_2:
    charge = Q_D
    layer_origin = np.array([lattice_constants[0]/3+a*np.cos(theta),a/2,a/2])+origin_offset
    lattice_locations = make_2D_layer(layer_origin, lattice_constants, charge, settings)
    lattice_locations_list.extend(lattice_locations)
    
    
    return lattice_locations_list


"""
-----------------------------------------------------------
                        INITIALIZATION
-----------------------------------------------------------
    

"""


# Defining constants and bundling them
epsilon_0 = 55.26349406*10**6 #e^2/eVum #8.8541878118*10**(-12)
charge = 1 #11.60218**(-19)
Angstrom = 10**(-10)


constants = {
    "epsilon_0":epsilon_0,
    "charge":charge,
    "Angstrom":Angstrom
    }

# Setting the size of the domain in Angstroms:
x_min = 0
x_max = 500
y_min = 0
y_max = 500
z_min = 0
z_max = 1000

particle_range = z_max

print("Particle begins " + '%5.1f'%(particle_range/20) + " nm below the junction")
print("Modeling sheet " + '%5.1f'%(x_max/10) + " nm by " + '%5.1f'%(y_max/10) + " nm")


x_middle = (x_max - x_min)/2
y_middle = (y_max - y_min)/2
z_middle = (z_max - z_min)/2

settings = {
                   "x_min":x_min,
                   "x_max":x_max,
                   "y_min":y_min,
                   "y_max":y_max,
                   "z_min":z_min,
                   "z_max":z_max
                   }

"""
-----------------------------------------------------------
                        BUILDING
-----------------------------------------------------------
    

"""

sub_2 = '\u2082'
sub_0 = '\u2080'

lattice_locations_list = []



Q_Mo  =   (42+4)
Q_Se  =   (34-2)
a     =   3.285
theta =   np.pi/6
print("Monolayer of MoSe" + sub_2 +  " at z = 0 Å")
origin_offset = np.array([1,1,z_middle])
lattice_locations_list = make_TMD_ML(lattice_locations_list, origin_offset, a, theta, Q_Mo, Q_Se, settings)

W_layer = False
if W_layer == True:
    Q_W  =   (74+4)
    Q_Se  =   (34-2)
    layer_separation = 6.5
    print("Monolayer of WSe" + sub_2 + " at z = " + '%5.2f Å'%layer_separation)
    origin_offset = origin_offset + np.array([a*np.cos(theta),a/2,layer_separation])
    lattice_locations_list = make_TMD_ML(lattice_locations_list, origin_offset, a, theta, Q_W, Q_Se, settings)

n_atoms = len(lattice_locations_list)
print('Number of atoms = ' + str(n_atoms))


built_time = time.time()



"""
-----------------------------------------------------------
                        PLOTTING
-----------------------------------------------------------
    

"""



# Plotting the locations of each atom in the 2D crystal

if n_atoms < 2000:
    Mo_x_coords = []
    Mo_y_coords = []
    Se_x_coords = []
    Se_y_coords = []
    
    for point in lattice_locations_list:
        if point[3] == 32:
    
            Se_x_coords.append(point[0])
            Se_y_coords.append(point[1])
            
        if point[3] == 46:
            Mo_x_coords.append(point[0])
            Mo_y_coords.append(point[1])
            
            
    plt.scatter(Mo_x_coords,Mo_y_coords,marker='x')
    plt.scatter(Se_x_coords,Se_y_coords,marker='o')
    
    plt.legend(['Mo','Se'])
    plt.title('MoSe2, Top View')
    plt.show()



z_resolution = 201
z = np.linspace(z_min,z_max, z_resolution)
midpoint_i = round(len(z)/2)
V = np.zeros(len(z))
test_x = x_middle
test_y = y_middle
total_charge = 0
for i in range(len(z)):
    test_r = np.array([test_x, test_y, z[i]])
    V[i],total_charge = calculate_V(test_r, np.array(lattice_locations_list), epsilon_0, charge, total_charge)

V_0 = -V[midpoint_i]
plt.plot(z-z_middle,-V)
if z_max>200:
    plt.ylim(1.1*V_0,0)
plt.title('Electric Potential V perpendicular to sheet')
plt.xlabel('z distance (Å)')
plt.ylabel(r'V (eV/$e$)')
plt.show()


U = V*2                         # 2e = charge of a Cooper pair => e*eV/e = eV
U_0 = U[midpoint_i]
plt.plot(z-z_middle,U)



c = 3*10**8
m = 2*0.5109989461 * 10**6/c**2      # eV/c^2
hbar = 6.582119569 * 10**(-16)  # eVs
k_B = 8.61733262*10**(-5)  # eV/K
T = 2 # K
v_thermal = np.sqrt(2*k_B*T/(np.pi*m))  # mean thermal velocity in 1D
v_z = v_thermal
E = (m*v_z**2)/2
print("Typical energy of thermal motion perpendicular to the sheet \nE = " + '%10.3E'%E + " eV at T = " + str(T)+" K")
#E = 0.8*U_0                   # eV


if z_max>200:
    plt.ylim(E*0.1,1.1*U_0)
plt.title('Electric Potential Energy U perpendicular to sheet')
plt.xlabel('z distance (Å)')
plt.ylabel('U (eV)')
plt.hlines(E,-z_max/2,z_max/2,colors='r')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.show()

dz = (z_max-z_min)/(z_resolution-1)*Angstrom # A/# * m/A = m
integral = 0

for i in range(len(z)):
    #if U[i] >= E and abs(z[i]-z_middle) <= particle_range:
    if U[i] >= E:
        integral += np.sqrt(U[i]-E)  # sqrt(eV)

integral = integral
#print('integral/dz = ' + str(integral))
exponential = -2/hbar*np.sqrt(2*m)*integral*dz  # 1/eVs * eVs/m * m = #
#print('exponential = ' + str(exponential))
probability = np.exp(exponential)  

print("Maximum value U"+ sub_0+" = " + '%10.3E'%U_0 + ' eV')
print('Probability that a Cooper pair tunnels through = ' + '%1.3f'%probability)

#print('Building took ' + str(built_time-start)[:5] + ' s')
V_time = time.time()
print('V calculation and plotting took ' + '%5.2f'%(V_time-built_time) + ' s')





