import project1Functions as lok
import matplotlib.pyplot as plt
import math
import numpy as np


# Task 1.1
def p(told, uold):
    return uold*3 
"""
[tgrid, approx, localErr] = lok.RK4walk(p,0,1,1,10)
realY = np.exp(tgrid)
plt.title("Task 1.1")

#plt.plot(tgrid,approx, label = "Approx")
#plt.plot(tgrid, realY, label = "Real")

glb = abs(realY - approx)

plt.loglog(tgrid, glb, label ="Global error")
plt.legend()
plt.show()
print(glb)
"""
"""

# Task 1.2
[unew, rnew] = lok.RK34step(p,0,1,0.01)
print(unew)
print(rnew)
"""
t0 = 0
y0 = 1
tf = 1
tol  = 10**(-6)
[t,y] = lok.adaptiveRK34(p, t0, tf, y0, tol)
plt.plot(t,y)
print(y[-1])
print(t[-1])
plt.show()