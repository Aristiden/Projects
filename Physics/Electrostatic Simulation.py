# -*- coding: utf-8 -*-
"""
Started on Apr 3, 2024
by Robert Wolle

Phys613 Final Project:  Ion Trap Electrostatic Modeling

The goal of this project was to model a particular ion trap in 2D and calculate the electric field in a vertical line at the center of the trap, thus finding where the particle should be constrained.
It can easily be extended to various geometries since it is implemented by drawing squares and triangles set to a particular voltage, and can even create triangular holes in other shapes.

"""

import time
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200 # setting the dpi of outputted graphs.
from numba import njit
from scipy import sparse

"""
-----------------------------------------------------------
                      FUNCTIONS
-----------------------------------------------------------

List of functions (not the order they appear in code):
    
 -- MAPPING ---------------------------------------------------------------------
      Functions that translate between numberings
        > get_s() -- Maps i,j to s numbering
        > get_i() -- Finds nearest i point corresponding to input x value
        > get_j() -- Finds nearest j point corresponding to input y value
 -- FINDING POINTS --------------------------------------------------------------
      Functions that find (i,j) points matching conditions and store them to or remove them from a list 
        > make_solid_object()    --  Makes list of (i,j) points inside given rectangular boundaries
        > make_triangle_object() --  Makes list of (i,j) points inside given right-triangle boundaries
        > remove_points()        --  Takes two lists of (i,j) points and removes any matching points from the first list
 -- CALCULATING U and A ------------------------------------------------------------
      Functions that take in lists of (i,j) points and set values of U, A at those points
        > border_U()  --  Sets border values of U vector to inputs
        > object_U()  --  Sets internal values of U vector corresponding to input list of points
        > hole_U()    --  Resets internal values of U vector to 0 corresponding to input list of points
        > object_A()  --  Sets the appropriate values of the A matrix for the interior domain
        > Neumann_A() --  Sets the appropriate values of the A matrix for Neumann boundary conditions
 -- APPLYING SHAPES -------------------------------------------------------------
      Front facing functions that get called in the main code.
        > apply_object_U()        --  Runs make_solid_object() and uses the output to run object_U()
        > apply_triangle_U()      --  Runs make_triangle_object() and uses the output to run object_U()
        > apply_triangle_hole_U() --  Runs make_triangle_object() and uses the output to run hole_U()
        > apply_object_A()        --  Runs object_A() which fills the interior structure of the A matrix, and then Neumann_A()
 -- PLOTTING --------------------------------------------------------------------
      Functions for plotting shapes.
        > plot_geometry() -- Plots 2 triangles based on given vertices defining those triangles

"""


"""
 -- MAPPING FUNCTIONS --------------
"""
@njit
def get_s(i,j,nx):
    """
    Maps u_i,u_j --> U_s
    """
    s = i + nx*j
    return s
@njit
def get_i(target,x_min,hx):
    """
    Finds the closest i value that would correspond to the given target value of x
    """
    i = round((target-x_min)/hx)
    return i
@njit
def get_j(target,y_min,hy):
    """
    Finds the closest j value that would correspond to the given target value of y
    """
    j = round((target-y_min)/hy)
    return j

"""
 -- FINDING POINTS ----------------
"""

@njit
def make_solid_object(top_boundary,bottom_boundary,left_boundary,right_boundary, nx,ny, hx,hy, x_min,y_min):
    """
    Makes a list of [i,j] points that are inside the object's boundaries
    """
    solid_ij = []
    
    top_j       = get_j(top_boundary,       y_min,hy)
    bottom_j    = get_j(bottom_boundary,    y_min,hy)
    left_i      = get_i(left_boundary,      x_min,hx)
    right_i     = get_i(right_boundary,     x_min,hx)
    
    # Iterate over all points within the boundary of the shape, and add each point to the list solid_ij:
    for i in range(left_i,right_i+1):
        for j in range(bottom_j,top_j+1):
            solid_ij.append([i,j])
            
    return solid_ij

@njit
def make_triangle_object(top_boundary,bottom_boundary,left_boundary,right_boundary, side, nx,ny, hx,hy, x_min,y_min):
    """
    Makes a list of [i,j] points that are inside the boundaries with the condition that:
        if side = 'left', then the left side is sloped:
            --> inside points are under the line between (left,bottom) and (right,top)
                                                         (x1,      y1)     (x2,    y2)
        if side = 'right', then the right side is sloped:
            --> inside points are under the line between (right,bottom) and (left,top)    
                                                         (x1,       y1)     (x2,   y2)
    """
    solid_ij = []
    
    top_j       = get_j(top_boundary,       y_min,hy)
    bottom_j    = get_j(bottom_boundary,    y_min,hy)
    left_i      = get_i(left_boundary,      x_min,hx)
    right_i     = get_i(right_boundary,     x_min,hx)
    
    if side == 'left':
        # m   = (y2  -  y1)/(x2  -  x1)                  ---->  rise/run
        slope = (top_boundary - bottom_boundary)/(right_boundary - left_boundary)
        for i in range(left_i,right_i+1):
            for j in range(bottom_j,top_j+1):
                # This is just recovering the (x,y) value corresponding to (i,j):
                x = i*hx + x_min
                y = j*hy + y_min
                
                # y  =      m(x - x1)            + y1    ---->  point slope form of a line
                line = slope*(x - left_boundary) + bottom_boundary
                if y <= line:
                    solid_ij.append([i,j])
                    
    # The algorithm is the same when the slope is on the right side of the shape,
    #  except "x1" and "x2" are switched.
    if side == 'right':
        slope = (top_boundary - bottom_boundary)/(left_boundary - right_boundary)
        for i in range(left_i,right_i+1):
            for j in range(bottom_j,top_j+1):
                x = i*hx + x_min
                y = j*hy + y_min
                
                line = slope*(x - right_boundary) + bottom_boundary
                if y <= line:
                    solid_ij.append([i,j])
    
    return solid_ij
        

def remove_points(solid_ij, hole_ij):
    """
    Takes in the points considered to be inside any object so far, solid_ij,
     and removes any points considered to be a hole in an object, hole_ij.
    """
    for point in hole_ij:
        if point in solid_ij:
            solid_ij.remove(point)
            
    # This is probably a horribly inefficient way to do this, and likely accounts for
    #  a large portion of the time spent building shapes.
            
    return solid_ij

"""
 -- CALCULATING U and A ----------------
"""

@njit
def border_U(U, V_top,V_bottom, V_left,V_right, nx,ny):
    """
    Avoids conditional statements by brute force replacing current values in vector U
     with input values of boundary voltages at each i,j along the extreme i,j values.
    """
    
    for i in range(nx):
        j = ny-1
        s = get_s(i,j,nx)
        U[s] = V_top
        j = 0
        s = get_s(i,j,nx)
        U[s] = V_bottom
    for j in range(ny):
        i = nx-1
        s = get_s(i,j,nx)
        U[s] = V_right
        i = 0
        s = get_s(i,j,nx)
        U[s] = V_left
    return U

@njit
def object_U(U, V_object, solid_ij, nx,ny):
    for point in solid_ij:
        # Each point is formatted as the list [i,j]
        # This step is uneccessary, but makes it more obvious what's happening
        i = point[0] 
        j = point[1]
        # Find the corresponding i,j --> s point
        s = get_s(i,j,nx)
        U[s] = V_object
    return U

@njit
def hole_U(U, solid_ij, nx,ny):
    """
    Resets value of U to 0 for a given list of (i,j) points solid_ij
    """
    for point in solid_ij:
        # Each point is formatted as the list [i,j]
        # These steps are uneccessary, but makes it more obvious what's happening
        i = point[0]
        j = point[1]
        # Find the corresponding i,j --> s point
        s = get_s(i,j,nx)
        U[s] = 0
    return U


def object_A(A, solid_ij, nx,ny, hx,hy):
    """
    Applies the 9-point Laplacian stencil for any point not in the list of points inside shapes.
        Input:
            A -- Current A matrix in Ax = U
            solid_ij -- List of points (i,j) inside the any shape's boundaries
            (nx,ny, hx,hy) -- Necessary for get_s() and calculating the Laplacian
            
        Output:
            solid_ij -- List of points (i,j) inside the input boundaries
            U -- Updated U vector with U[s] = V_object for every point (i,j) in solid_ij
    
    """
    s_solid = []
    for point in solid_ij:
        i = point[0]
        j = point[1]
        s_solid.append(get_s(i,j,nx))    
    
    
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            s = get_s(i,j,nx)
            
            # 5-Point stencil, unused
            """if s not in s_solid:
                A[s][s]=-4/hx**2
                A[s][s+nx]=1/hx**2
                A[s][s-nx]=1/hx**2
                A[s][s+1]=1/hx**2
                A[s][s-1]=1/hx**2
            """
            # 9-Point stencil
            if s not in s_solid:
                A[s,s]= -20/hx**2
                A[s,s+nx]= 4/hx**2
                A[s,s-nx]= 4/hx**2
                A[s,s+1] = 4/hx**2
                A[s,s-1] = 4/hx**2
                A[s,s+1+nx] = 1/hx**2
                A[s,s-1+nx] = 1/hx**2
                A[s,s+1-nx] = 1/hx**2
                A[s,s-1-nx] = 1/hx**2
    
    return A, s_solid

def Neumann_A(A, s_solid, nx,ny, hx,hy):
    """
    Creates Neumann boundaries at i = 0,nx-1 and j = 0,ny-1
    """
    
    
    for j in range(ny):
        i = 0
        s = get_s(i,j,nx)
        A[s,s] += 3/(2*hx)
        A[s,s+1] += -2/hx
        A[s,s+2] += 1/(2*hx)
        i = nx-1
        s = get_s(i,j,nx)
        A[s,s] += 3/(2*hx)
        A[s,s-1] += -2/hx
        A[s,s-2] += 1/(2*hx)       
    for i in range(nx):
        j = 0
        s = get_s(i,j,nx)
        A[s,s] += 3/(2*hx)
        A[s,s+nx] += -2/hx
        A[s,s+2*nx] += 1/(2*hx)
        j = ny-1
        s = get_s(i,j,nx)
        A[s,s] += 3/(2*hx)
        A[s,s-nx] += -2/hx
        A[s,s-2*nx] += 1/(2*hx)
    
    
    
        
    return A
    
"""
 -- APPLYING SHAPES --------------------------
"""

@njit
def apply_object_U(U, V_object, top,bottom,left,right, nx,ny, hx,hy, x_min,y_min):
    """
    Makes solid objects with square boundaries.
        Input:
            U -- Current U vector in Ax = U
            V_object -- Voltage of the object
            (top,bottom,left,right) -- Boundaries of the object
            (nx,ny, hx,hy, x_min,y_min) -- Necessary for get_i(), get_j()
            
        Output:
            solid_ij -- List of points (i,j) inside the input boundaries
            U -- Updated U vector with U[s] = V_object for every point (i,j) in solid_ij
    """
    
    solid_ij    =  make_solid_object(top,bottom,left,right, nx,ny, hx,hy, x_min,y_min)
    U = object_U(U, V_object, solid_ij, nx,ny)
    
    return U, solid_ij

@njit
def apply_triangle_U(U, V_object, top,bottom,left,right, side, nx,ny, hx,hy, x_min,y_min):
    """
    Makes solid objects with right-triangle boundaries.
        Input:
            U -- Current U vector in Ax = U
            V_object -- Voltage of the object
            (top,bottom,left,right) -- Boundaries of the object
            side -- String, either 'left' or 'right' indicating which side of the triangle is sloped
            (nx,ny, hx,hy, x_min,y_min) -- Necessary for get_i(), get_j()
            
        Output:
            solid_ij -- List of points (i,j) inside the input triangle boundaries
            U -- Updated U vector with U[s] = V_object for every point (i,j) in solid_ij
    """
    solid_ij    =  make_triangle_object(top,bottom,left,right, side, nx,ny, hx,hy, x_min,y_min)
    U           =  object_U(U, V_object, solid_ij, nx,ny)
    
    return U, solid_ij

@njit
def apply_triangle_hole_U(U, top,bottom,left,right, side, nx,ny, hx,hy, x_min,y_min):
    """
    Makes right-triangle holes in solid objects.
        Input:
            U -- Current U vector in Ax = U
            (top,bottom,left,right) -- Boundaries of the object
            side -- String, either 'left' or 'right' indicating which side of the triangle is sloped
            (nx,ny, hx,hy, x_min,y_min) -- Necessary for get_i(), get_j()
            
        Output:
            hole_ij -- List of points (i,j) inside the input triangle boundaries
            U -- Updated U vector with U[s] = 0 for every point (i,j) in solid_ij
    """
    
    hole_ij     =  make_triangle_object(top,bottom,left,right, side, nx,ny, hx,hy, x_min,y_min)
    U           =  hole_U(U, hole_ij, nx,ny)
    return U, hole_ij




def apply_object_A(A, solid_ij, nx,ny, hx,hy):
    """
    Makes the matrix A.
     Excludes points inside the solid objects given by solid_ij using object_A().
     Creates Neumann boundaries for A with Neumann_A()
    """
    
    A, s_solid = object_A(A, solid_ij, nx,ny, hx,hy)
    A = Neumann_A(A, s_solid, nx,ny, hx,hy)
    
    return A

"""
 -- PLOTTING FUNCTIONS --------------------
"""

def plot_geometry(segments, x_max):
    """
    Plots 2 triangles built out of black lines defining the true edges of the trap geometry.
     The segments are formatted as:
    
    outer_trap_tip          =   segments[0]
    outer_trap_edge_top     =   segments[1]
    outer_trap_edge_bottom  =   segment[2]
    
    inner_trap_tip          =   segments[3]
    inner_trap_tip_bottom   =   segments[4]
    inner_trap_slope_bottom =   segments[5]
    
     with segments[segment] = [x,y] corresponding to that vertex.
    Also mirrors the triangles across the y-axis using at x_max/2.
    """
    linewidth = 0.7
    
    # Outer trap segments
    segment = 0
    x_values = np.array([segments[segment][0], segments[segment+1][0], segments[segment+2][0],segments[segment][0]])
    y_values = np.array([segments[segment][1], segments[segment+1][1], segments[segment+2][1],segments[segment][1]])
    plt.plot(x_values,y_values,'k-', linewidth=linewidth)
    plt.plot(x_max-x_values,y_values,'k-', linewidth=linewidth) # Mirroring the triangle
    
    # Inner trap segments
    segment = 3 # This just redefines which vertex the triangle begins on
    x_values = np.array([segments[segment][0], segments[segment+1][0], segments[segment+2][0],segments[segment][0]])
    y_values = np.array([segments[segment][1], segments[segment+1][1], segments[segment+2][1],segments[segment][1]])
    plt.plot(x_values,y_values,'k-', linewidth=linewidth)
    plt.plot(x_max-x_values,y_values,'k-', linewidth=linewidth)



"""
-----------------------------------------------------------
                      INITIALIZATION
-----------------------------------------------------------

This part of the code defines the problem's domain, grid size, and initial values of 
 A and U in the linear system Ax = U.

"""


#Initialize time:
start = time.time()

# Setting the number of grid points on each side to be the same, nx=ny:
nx = 351
ny = 351

#   Cost of dense (nx*ny x nx*ny) matrices:
#        nx   ny
#       101, 101 costs 0.79 GB of RAM
#       121, 121 costs 1.60 GB of RAM
#       151, 151 costs 3.87 GB of RAM


# Setting the size of the domain (e.g. length in millimeters):
x_min = 0
x_max = 2.2
y_min = 0
y_max = 2.2

# Grid size:
hx = (x_max - x_min)/(nx-1)
hy = (y_max - y_min)/(ny-1)

# Setting the boundary voltages (in Volts):
V_top     =  0
V_bottom  =  0
V_left    =  0
V_right   =  0


U = np.zeros(nx*ny)

# Only sets U values at the boundary if any of the boundary voltages are not 0.
#  This is done to cut down on initialization time.
V_check = V_top**2 + V_bottom**2 + V_left**2 + V_right**2
if V_check != 0:
    # Sets up U at the boundaries with given boundary voltages:
    U = border_U(U, V_top, V_bottom, V_left, V_right, nx,ny)


# We can set up the linear system by beginning with an identity matrix 
#  where every point in the solution is initially equal to the corresponding value of U
A = sparse.lil_matrix((nx*ny,nx*ny))
A.setdiag(1)


border_time = time.time()
print('Initialization took ' + str(border_time-start) + ' s')

"""
-----------------------------------------------------------
                        BUILDING
-----------------------------------------------------------

This part of the code handles the creation of shapes inside the boundary.
It works in steps:
    

    1. Defintions
        a. Voltages of the shapes:
            >  V_trap   -- Voltage of the outer trap
            >  V_ground -- Voltage of the inner trap
        b. Vertices of the shapes
            >  outer_trap_tip -- Position of the upper tip of the outer trap
            >  outer_trap_edge_top -- Position where the top of the outer trap stops
                                      along x_min
            >  outer_trap_edge_bottom -- Position where the bottom of the outer trap stops
                                         along x_min
            >  outer_trap_correction -- Position of the top right corner of a box meant
                                        to fill the space left empty by a slope starting
                                        at x_min.
            >  inner_trap_tip -- Position of the upper tip of the inner trap
            >  inner_trap_tip_bottom -- Position of the upper tip of the inner trap along y_min
            >  inner_trap_slope_bottom -- Position of the inner trap's slope edge along y_min 
        c. Empty python array, to be filled with points inside each shape
        
    2. Building Shapes
        a. Outer trap
            i.   Outer trap slope -- The topmost part of the outer trap
            ii.  Outer trap correction -- The part below the outer trap slope that is
                                          left unfilled when beginning the slope from x_min
            iii. Outer trap slope hole -- The part of the slope that needs to be removed
                                          to get the sloped bottom of the outer trap.
            iv.  Outer trap (mirrored) -- Shapes i-iii. mirrored across x_max/2
        b. Inner trap
            v.   Inner trap slope -- The slope that makes up the inner trap
            vi.  Inner trap slope (mirrored) -- Shape v. mirrored across x_max/2
            
            
    Building shapes requires two parts:
        1. Defining the shape:
            a. The boundaries and voltages of the shape are already defined, but need
                to be set for each new shape.
        2. Applying the shape:
            a. Set U[s] = V_object for every (i,j) --> s in the defined shape.
            b. Add the list of every (i,j) inside the shape to a list of all points
                inside shapes, or remove the list of points inside holes in shapes.
            c. Apply the 9-point Laplacian stencil to create the A matrix, excluding
                points in the list of every point inside shapes. 
                
    The format of this process looks like:
        
        Name of shape being built
        
        V_object        =  Voltage of the shape
        right, top      =  boundaries of the shape
        left, bottom    =  boundaries of the shape
        side            =  if the shape is a triangle, which side is the slope on
        U, list_ij      =  apply_shape_U(...) -- Finds all points in the shape and sets 
                                                 their corresponding U appropriately
        solid_ij.extend(list_ij) -- Stores all points inside of the shape to the list 
                                    of all points inside shapes
           or (if shape is a hole)
           
        solid_ij = remove_points(solid_ij,list_ij) -- Removes any point in the shape
                                                      from the list of all points inside shapes

    Then, after all shapes are built:
        A = apply_object_A(...)


"""



# Defining features of geometry of the trap:
V_trap      =   250
V_ground    =     0

#  Name of vertex                  x           y
outer_trap_tip          =   [x_min + 0.3, y_min + 0.75]
outer_trap_edge_top     =   [x_min,       y_min + 0.41]
outer_trap_edge_bottom  =   [x_min,       y_min + 0.25]
outer_trap_correction   =   [outer_trap_edge_top[0] + 0.09, outer_trap_edge_top[1] - 0.001]

inner_trap_tip          =   [x_min + 0.9, y_min + 0.54]
inner_trap_tip_bottom   =   [inner_trap_tip[0],  y_min]
inner_trap_slope_bottom =   [x_min + 0.66,       y_min]

segments = [outer_trap_tip, outer_trap_edge_top, outer_trap_edge_bottom, inner_trap_tip, inner_trap_tip_bottom, inner_trap_slope_bottom]

# Empty list to be filled with points inside of shapes:
solid_ij = []

"""
Outer trap slope - left  /
"""
V_object        =   V_trap
right, top      =   outer_trap_tip
left, bottom    =   outer_trap_edge_top
side            =   'left'
U, list_ij      =   apply_triangle_U(U, V_object, top,bottom,left,right, side, nx,ny,hx,hy,x_min,y_min)
solid_ij.extend(list_ij)


"""
Outer trap correction - left
"""
V_object        =   V_trap
right, top      =   outer_trap_correction
left, bottom    =   outer_trap_edge_bottom
U, list_ij      =   apply_object_U(U, V_object, top,bottom,left,right, nx,ny,hx,hy,x_min,y_min)
solid_ij.extend(list_ij)


"""
Outer trap slope HOLE - left  /
"""
right, top      =   outer_trap_tip
left, bottom    =   outer_trap_edge_bottom                                                      
side            =   'left'
U, list_ij      =   apply_triangle_hole_U(U, top,bottom,left,right, side, nx,ny,hx,hy,x_min,y_min)
solid_ij = remove_points(solid_ij,list_ij)


"""
Outer trap slope - right  \
"""
V_object        =   V_trap
left, top       =   outer_trap_tip
right, bottom   =   outer_trap_edge_top
left            =   x_max - left
right           =   x_max - right
side            =   'right'
U, list_ij      =   apply_triangle_U(U, V_object, top,bottom,left,right, side, nx,ny,hx,hy,x_min,y_min)
solid_ij.extend(list_ij)


"""
Outer trap correction - right
"""
V_object        =   V_trap
left, top       =   outer_trap_correction
right, bottom   =   outer_trap_edge_bottom
left            =   x_max - left
right           =   x_max - right
U, list_ij      =   apply_object_U(U, V_object, top,bottom,left,right, nx,ny,hx,hy,x_min,y_min)
solid_ij.extend(list_ij)


"""
Outer trap slope HOLE - right  \
"""
left, top       =   outer_trap_tip
right, bottom   =   outer_trap_edge_bottom 
left            =   x_max - left
right           =   x_max - right
side            =   'right'
U, list_ij      =   apply_triangle_hole_U(U, top,bottom,left,right, side, nx,ny,hx,hy,x_min,y_min)
solid_ij = remove_points(solid_ij,list_ij)


"""
Inner trap slope - left  /
"""
V_object        =   V_ground
right, top      =   inner_trap_tip
left, bottom    =   inner_trap_slope_bottom
side            =   'left'
U, list_ij      =   apply_triangle_U(U, V_object, top,bottom,left,right, side, nx,ny,hx,hy,x_min,y_min)
solid_ij.extend(list_ij)


"""
Inner trap slope - right  \
"""
V_object        =   V_ground
left, top       =   inner_trap_tip
right, bottom   =   inner_trap_slope_bottom
left            =   x_max - left
right           =   x_max - right
side            =   'right'
U, list_ij      =   apply_triangle_U(U, V_object, top,bottom,left,right, side, nx,ny,hx,hy,x_min,y_min)
solid_ij.extend(list_ij)




A = apply_object_A(A, solid_ij, nx,ny,hx,hy)

built_time = time.time()
print('Building took ' + str(built_time- border_time) + ' s')



"""
-----------------------------------------------------------
             SOLVING and PLOTTING POTENTIAL
-----------------------------------------------------------

This portion of the code solves the linear system Ax = U

"""




# First, I convert the A matrix into a sparse matrix:
A_sparse = sparse.csr_matrix(A)
# Then, I use a sparse solver from SciPy for a gigantic speedup over NumPy's dense solver
x = sparse.linalg.spsolve(A_sparse,U)


# The solution is in s-ordering, and can be rearranged to display the result in (x_i, y_j) format.
#  This arranged solution has rows corresponding to y points and columns to x points.
V = np.zeros([nx,ny])
for i in range(nx):
    for j in range(ny):
        s = get_s(i, j, nx)
        V[j][i] = x[s]

"""
 -- PLOTTING ---------------------

  1. Plotting contours and geometry of the electric potential.
  2. Plotting electric potential V / V_trap over distance away from trap.
"""
# Plotting line segments of the geometry:
plot_geometry(segments, x_max)

# Plotting the 2D electric potential:
X = np.linspace(0,x_max,nx)
Y = np.linspace(0,y_max,ny)
contours = 300

plt.contourf(X,Y,V, levels=contours)
plt.colorbar()
plt.title('Trap Potential V, n = ' + str(nx))
plt.show()


# Plotting the potential V as a function of vertical distance
midpoint_i = round(nx/2)
plt.plot(Y,V[:,midpoint_i]/V_trap, 'k', linewidth=2)
plt.xlim(0,2.2)

# Zoomed to the potential maximum:
#plt.xlim(0.75,1.75)
#plt.ylim(0.175,0.27)

plt.ylabel(r'V / V$_{out}$')
plt.xlabel(r'vertical distance (mm)')
plt.grid()
plt.show()

# Plotting the potential V as a function of horizontal distance at the potential maximum
trap_point_j = get_j(1.444,y_min,hy)

plt.plot(Y,V[trap_point_j,:]/V_trap, 'k', linewidth=2)
plt.xlim(0,2.2)
plt.ylabel(r'V / V$_{out}$')
plt.xlabel(r'horizontal position (mm)')
plt.grid()
plt.show()



V_time = time.time()
print('V solution and plotting took ' + str(V_time-built_time) + ' s')



"""
-----------------------------------------------------------
             CALCULATING and PLOTTING FIELD
-----------------------------------------------------------

This portion of the code calculates the gradient of the found potential
 using the central-difference approximation. This gradient is the E field.
The norm of the E field is then taken at every point to find the point where
 |E| = 0, which is the location where a charged particle will be trapped.
The value of |E| is then clipped at problematic points so that the minimum 
 value is more visible.
Finally, |E| is plotted.
 
 
"""

# Calculating the electric field as the gradient of the electric potential:
# E = -grad(V)
E = np.zeros((nx-2,ny-2,2))
for i in range(1,nx-1):
    for j in range(1,ny-1):
        # Calculating partial derivatives as central-difference:
        E[j-1][i-1][0] = -(V[j][i+1] - V[j][i-1])/(2*hx)
        E[j-1][i-1][1] = -(V[j+1][i] - V[j-1][i])/(2*hy)

# Calculating the 2-norm of E at every point in space:
normE = np.zeros((nx-2,ny-2))
for i in range(nx-2):
    for j in range(ny-2):
        normE[j][i] = np.linalg.norm(E[j][i])



# Clipping off values at sharp geometry to display the trap location clearly:
halfrange = round(ny/2)
clip_value = 1400
for i in range(nx-2):
    for j in range(halfrange):
        if normE[j][i] > clip_value:
            normE[j][i] = clip_value

        
# Plotting |E|
X = np.linspace(hx,x_max-hx,nx-2)
Y = np.linspace(hy,y_max-hy,ny-2)

plot_geometry(segments, x_max)
contours = 300
plt.contourf(X,Y,normE, levels=contours, cmap='gist_ncar')
plt.colorbar()
plt.xlim(0,2.2)
plt.ylim(0,2.2)
plt.title(r'Trap Field |E|, n = ' + str(nx))
plt.show()

# Finding and printing the distance of the minimum away from the trap tip:
min_j_normE = np.argmin(normE[get_j(outer_trap_tip[1],y_min,hy):get_j(y_max-0.2,y_min,hy),midpoint_i])
min_y_normE = str(min_j_normE*hy)
print('Minimum |E| at ' + min_y_normE[:6] + ' mm away from trap tip')

E_time = time.time()
print('Calculating E and plotting took ' + str(E_time-V_time) + ' s')

print('Total runtime = ' + str(E_time-start) + ' s')




