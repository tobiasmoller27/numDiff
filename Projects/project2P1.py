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





#plt.plot(x,y)
#plt.loglog(x, err)
#plt.show()
Dxs = np.zeros(400-10)
rms = np.zeros(400-10)

for N in range(10, 400):
    Dx = L/(N+1)
    Dxs[N-10] = Dx
    x = np.linspace(0+Dx,L-Dx,N)
    y = lok.twopBVP(p(x), 1, math.e**4, L, N)
    x = np.linspace(0, L, N+2)
    yReal = np.exp(np.square(x))
    err = np.zeros(len(y))
    for i in range(len(y)):
        err[i]=(abs(yReal[i] - y[i])**2*Dx)
    rms[N-10] = math.sqrt(np.sum(err))

plt.title("RMS of Local Error per Dx")
plt.xlabel("Dx")
plt.ylabel("RMS of Local Error")
plt.grid()
plt.loglog(Dxs, Dxs**2, label="Reference")
plt.loglog(Dxs, rms, label = "Error")
plt.legend()
plt.show()
"""


# Task 1.2

#Start by calculating M
L = 10
N = 999
q = (-50000)*np.ones((N,1))
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
plt.title("Solution of the Beam Equation")
plt.ylabel("u(x)")
plt.xlabel("x")
plt.plot(x,u)
plt.show()
