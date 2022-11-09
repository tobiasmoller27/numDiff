from cmath import exp
import re
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as linalg

def RK4step(f, told, uold, h):
    Y1 = f(told,uold)
    Y2 = f(told+h/2,uold+h*Y1/2)
    Y3 = f(told+h/2,uold+h*Y2/2)
    Z3 = f(told+h,uold-h*Y1+2*h*Y2)
    Y4 = f(told+h,uold+h*Y3)

    return [uold+(h/6)*(Y1+2*Y2+2*Y3+Y4), (h/6)*(2*Y2+Z3-2*Y3-Y4)]

def RK4walk(f,t0,y0,tf,N):
    h=(tf-t0)/N
    uold = y0
    told = t0
    unew = uold
    tgrid = np.arange(t0,tf+h,h)
    localErr = np.zeros(N+1)
    approx = np.zeros(N+1)
    l=0

    for i in range(N+1):
        approx[i] = unew
        localErr[i] = l
        [unew,l]=RK4step(f,told,uold,h)
        uold = unew
        

    return [tgrid, approx, localErr]


def RK34step(f,told,uold,h):
    [unew,error] = RK4step(f,told,uold,h)
    rNPlus1 = linalg.norm(error)
    return [unew, rNPlus1]


def newstep(tol, err, errold, hold, k):
    hnew = ((tol/err)**(2/(3*k)))*((tol/errold)**(-1/(3*k)))*hold
    return hnew


def adaptiveRK34(f, t0, tf, y0, tol):
    h = (abs(tf-t0)*tol**(1/4))/(100*(1+linalg.norm(f(t0,y0))))
    k = 4
    errold = tol
    uold = y0
    unew = y0
    tc = t0
    t = []
    y = []

    while(tc <= tf):
        t.append(tc)
        y.append(unew)
        [unew, err] = RK34step(f, tc, uold, h)
        if err > tol:
            h = newstep(tol, err, errold, h, k)
        errold = err
        uold = unew
        tc += h
    hDiff = tf-t[-2]
    [yf, err] = RK34step(f, t[-2], y[-2], hDiff)
    y[-1] = yf
    t[-1] = tf
    
    return [t, y]