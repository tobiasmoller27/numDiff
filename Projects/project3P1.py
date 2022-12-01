import numpy as np
import scipy.linalg as lg
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import project3Functions as lok


# TASK 1.1

#Constructing the Toeplitzmatrix as in Project 2
#Tdx = np.diag(np.full(N,-2))+np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1)
def g(x):
    return np.sin(x)

tstart, tend, M, N = 0, 1, 100, 100
y, x, t = lok.eulerint(g,tstart,tend,M,N)

fig = plt.figure()
ax = plt.axes(projection='3d')


plt.show()


"""
Dont forget to insert the boundary values

This functions should be used:
- plot_wireframe
- plot_surface
"""