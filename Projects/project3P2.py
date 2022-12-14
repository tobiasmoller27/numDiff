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
tend = 1
N = 10
M = 1000
a = 1

y, x, t, norm = lok.LaxWenInt(g, a, tstart, tend, M, N)


fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(x, t, y)
ax.set_title('wireframe')
ax.set_ylabel('X')
ax.set_xlabel('T')


#plt.plot(np.linspace(tstart, tend, M+1),norm)
plt.show()