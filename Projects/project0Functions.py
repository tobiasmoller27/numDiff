from cmath import exp
import re
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg


def eulerstep(A, uold, h):
    return uold + h*A*uold 

def eulerint(A, y0, t0, tf, N):
    h=(tf-t0)/N
    uold = y0
    tgrid = np.zeros(N)

    for i in range(N):
        unew = eulerstep(A, uold, h)
        uold = unew
        tgrid[i] = unew

    approx=unew
    err= linalg.norm(approx - linalg.expm(tf*A))
    return [tgrid, approx, err]

def errVSh(A, y0, t0, tf):
    NArray = np.arange(10,500)
    errAxis = np.zeros(len(NArray))
    h = (tf-t0)/NArray
    for i in range(len(NArray)):
        [tgrid, approx, err] = eulerint(A, y0, t0, tf, NArray[i])
        errAxis[i] = err
    plt.title("Error per h")
    plt.loglog(h, errAxis)
    plt.show()

"""
N = 10000
[tgrid, approx, err] = eulerint(np.matrix(2), 1, 0, 1, N)
print("Error: " + str(err))

t = np.linspace(0,1,N)
realY = np.exp(2*t)
plt.title("Task 1.2")
plt.plot(t, tgrid, label = "Approximation")
plt.plot(t, realY, label = "Real function")
plt.legend()
plt.show()
"""

errVSh(np.matrix(2), 1, 0, 1)