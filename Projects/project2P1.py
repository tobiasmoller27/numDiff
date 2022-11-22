import project2P1Functions as lok
import numpy as np
import math
import matplotlib.pyplot as plt

"""
# TASK 1.1
N = 100
L = 2
Dx = L/(N+1)
#test function y''
def p(x):
    return (2+4*(x**2))*np.exp(x**2)
x = np.linspace(0+Dx,L-Dx,N)
f = p(x)

y = lok.twopBVP(f, 1, math.e**4, L, N)
x = np.linspace(0, L, N+2)
plt.plot(x,y)
plt.show()
"""

# Task 1.2

#Start by calculating M
L = 10
N = 999
q = (-50)*np.ones((N,1))
Dx = L/(N+1)
xInterior = np.linspace(0+Dx, L-Dx, N)
M = lok.twopBVP(q,0,0,L,N)
x = np.linspace(0, L, N+2)
"""plt.title("M(x)")
plt.plot(x,M)
plt.show()"""

#Calculating u
I = 0.001*(3-2*(np.cos(math.pi*xInterior/L))**12)
E = 1.9 * 10**11
MInterior = M[1:N+1].reshape((N,))
u2 = MInterior/(E*I)
u = lok.twopBVP(u2, 0, 0, L, N)
print("u midpoint: "+str(u[501]))
plt.title("u(x)")
plt.plot(x,u)
plt.show()