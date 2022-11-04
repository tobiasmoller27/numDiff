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
    unew = y0
    tgrid = np.arange(t0,tf,1/N)
    err = np.zeros(N)
    approx = np.zeros((2,N))

    for i in range(N):
        r, c = 0, i
        approx[r:r+unew.shape[0], c:c+unew.shape[1]] = unew
        err[i] = linalg.norm(unew - linalg.expm(tgrid[i]*A))
        uold = unew
        unew = eulerstep(A, uold, h)
    return [tgrid, approx, err]

def errVSh(A, y0, t0, tf):
    NArray  = np.arange(10,150)
    endpointErr = np.zeros(len(NArray))
    h = (tf-t0)/NArray
    for i in range(len(NArray)):
        [tgrid, approx, err] = eulerint(A, y0, t0, tf, NArray[i])
        endpointErr[i] = err[NArray[i]-1]
    return h, endpointErr
    """
    plt.title("Error per h")
    plt.loglog(h, endpointErr)
    plt.show()
    """

def ieulerstep(A, uold, h):
    return uold + h*A*uold*+(h**2)*(A**2)*uold

def ieulerint(A, y0, t0, tf, N):
     h=(tf-t0)/N
     uold = y0
     unew = y0
     tgrid = np.arange(t0,tf,1/N)
     err = np.zeros(N)
     approx = np.zeros((2,N))
     for i in range(N):
        r, c = 0, i
        approx[r:r+unew.shape[0], c:c+unew.shape[1]] = unew
        err[i] = linalg.norm(unew - linalg.expm(tgrid[i]*A))
        uold = unew
        unew = eulerstep(A, uold, h)
     return [tgrid, approx, err]


