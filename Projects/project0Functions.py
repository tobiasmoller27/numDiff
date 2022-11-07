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
    tgrid = np.arange(t0,tf+h,h)
    err = np.zeros(N+1)
    approx = np.zeros((2,N+1))

    for i in range(N+1):
        r, c = 0, i
        approx[r:r+unew.shape[0], c:c+unew.shape[1]] = unew
        err[i] = linalg.norm(unew - linalg.expm(tgrid[i]*A)*y0)
        uold = unew
        unew = eulerstep(A, uold, h)
    return [tgrid, approx, err]

def errVSh(A, y0, t0, tf):
    NArray  = np.arange(1,2**9+1)
    endpointErr = np.zeros(len(NArray))
    h = (tf-t0)/NArray
    for i in range(len(NArray)):
        [tgrid, approx, err] = eulerint(A, y0, t0, tf, NArray[i])
        endpointErr[i] = err[-1]
    return h, endpointErr
    """
    plt.title("Error per h")
    plt.loglog(h, endpointErr)
    plt.show()
    """

### IMPLICIT BELOW ###
def ieulerstep(A, uold, h):
    return linalg.inv(np.identity(len(A))-h*A)*uold

def ieulerint(A, y0, t0, tf, N):
    h=(tf-t0)/N
    uold = y0
    unew = y0
    tgrid = np.arange(t0,tf+h,h)
    err = np.zeros(N+1)
    approx = np.zeros((2,N+1))

    for i in range(N+1):
        r, c = 0, i
        approx[r:r+unew.shape[0], c:c+unew.shape[1]] = unew
        err[i] = linalg.norm(unew - linalg.expm(tgrid[i]*A)*y0)
        uold = unew
        unew = ieulerstep(A, uold, h)
    return [tgrid, approx, err]


def ierrVSh(A, y0, t0, tf):
    NArray  = np.arange(1,2**9+1)
    endpointErr = np.zeros(len(NArray))
    h = (tf-t0)/NArray
    for i in range(len(NArray)):
        [tgrid, approx, err] = ieulerint(A, y0, t0, tf, NArray[i])
        endpointErr[i] = err[-1]
    return h, endpointErr


    """
    plt.title("Error per h")
    plt.loglog(h, endpointErr)
    plt.show()
    """
