import numpy as np
import scipy.linalg as lg
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import project3Functions as lok
import math

def g(x):
    return np.exp(-100*(x-0.5)**2)
    #return x
    #return x-x**2-x**3
    

tstart = 0 
tend = 5
N = 900
M = 1000
a = 0.2

y1, x1, t1, norm = lok.LaxWenInt(g, a, tstart, tend, M, N)

"""
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = plt.axes(projection='3d')

ax.plot_surface(x1, t1, y1)
ax.set_title('Advection equation')
ax.set_xlabel('X')
ax.set_ylabel('T')
plt.show()
"""
"""
#First subplot
y1, x1, t1, norm = lok.LaxWenInt(g, 1, 0, 1, M, N)
ax = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax.plot_surface(x1,t1,y1)

#Second subplot
y2, x2, t2, norm = lok.LaxWenInt(g, 0.2, 0, 5, M, N)
ax = fig.add_subplot(1, 2, 2, projection='3d')
surf2 = ax.plot_surface(x2,t2,y2)
"""

plt.plot(np.linspace(tstart, tend, M+1),norm)
plt.title('Norm against time')
plt.xlabel('Time')
plt.ylabel('Norm')
plt.show()