from cmath import exp
import re
import numpy as np
import matplotlib as plt
import scipy.linalg as linalg


def eulerstep(A, uold, h):
    return uold + h*A*uold 

def eulerint(A, y0, t0, tf, N):
    h=(tf-t0)/N
    uold = y0
    tgrid = np.zeros(shape=(1,N))

    for i in range(N):
        unew = eulerstep(A, uold, h)
        uold = unew
        tgrid[1,i] = unew

    approx=unew
    err= approx - linalg.expm(tf*A)
    return [tgrid, approx, err]

[tgrid, approx, err] = eulerint(2, 1, 0, 1, 100)
