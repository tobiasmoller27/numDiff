import numpy as np
import scipy.linalg as lg
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import project3Functions as lok
import math


# TASK 1.1

#Constructing the Toeplitzmatrix as in Project 2
#Tdx = np.diag(np.full(N,-2))+np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1)
def g(x):
    return (0.5-abs(x-0.5))
    #return x
    #return (np.cos(math.pi*x)**2)

tstart, tend, M, N = 0, 1, 500, 11
y, x, t = lok.eulerint(g,tstart,tend,M,N)

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(x, t, y)
ax.set_title('Diffusion equation')
ax.set_xlabel('X')
ax.set_ylabel('T')
plt.show()


"""
Dont forget to insert the boundary values

This functions should be used:
- plot_wireframe
- plot_surface
"""