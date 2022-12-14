import numpy as np
import scipy.linalg as lin
import matplotlib.pyplot as plt
import project3Functions as lok
from mpl_toolkits.mplot3d import Axes3D

def g(x):
    return np.exp(-100*(x-0.5)**2)
    #return x
    #return np.exp(x)
    

tstart = 0 
tend = 1
N = 500
M = 500
a = 5
d = 0.1

print(a/d)

y, x, t = lok.convdifInt(g, a, d, tstart, tend, M, N)

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(x, t, y)
ax.set_title('wireframe')
ax.set_ylabel('X')
ax.set_xlabel('T')
plt.show()