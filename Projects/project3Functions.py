import numpy as np
import scipy.linalg as lin
import math
import matplotlib.pyplot as plt
import project3Functions as lok
from scipy.sparse import diags


def eulerstep(Tdx, uold, dt):
    return uold + dt*np.matmul(Tdx,uold)


def eulerint(g, tstart, tend, M, N):
    #Defining timeinterval and timegrid
    dt=tend/(M+1)
    tt = np.linspace(0, tend, M+1)
    
    #Defining x-grid
    dx = 1/(N+1)
    xx = np.linspace(0, 1, N+2)
    xInterior = xx[1:-1]

    #Constructing the meshgrid
    [T,X] = np.meshgrid(tt, xx)
    #Constructing the Toeplitzmatrix as in Project 2
    Tdx = (np.diag(np.full(N,-2))+np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1))/(dx**2)    
    print(Tdx)
    #Constructing solution matrix

    solution = np.zeros((N,M+1))
    solution[:,0] = g(xInterior)
    

    uold = g(xInterior)
    for c in range (M):
        #unew =eulerstep(Tdx, uold, dt)
        unew =TRstep(Tdx, uold, dt)
        solution[:,c+1] = unew
        uold = unew

    solution = np.vstack((np.zeros((1,M+1)),solution,np.zeros((1,M+1))))
    print(dt/(dx**2))
    return solution, X, T
    
def TRstep(Tdx, uold, dt):

    I = np.identity(len(Tdx))
    A = (I - (dt/2) * Tdx)
    B = np.matmul((I+(dt/2)*Tdx) , uold)
    unew = np.linalg.solve(A,B)
    return unew

def LaxWen(u, amu):
    A = 1-amu**2
    B = (amu/2)*(amu + 1)
    C = (amu/2)*(amu - 1)

    N = len(u)

    S = (np.diag(A*np.full(N,1))+np.diag(C*np.ones(N-1),1)+np.diag(B*np.ones(N-1),-1))
    
    S[N-1,0] = C
    S[0,N-1] = B

    
    unew = np.matmul(S,u)

    return unew

def LaxWenInt(g, a, tstart, tend, M, N):
    dx =1/(N+1)
    xx = np.linspace(0, 1, N)
    xInterior = xx

    dt = tend/(M+1)
    tt = np.linspace(tstart,tend,M+1)

    [T,X] = np.meshgrid(tt, xx)

    solution = np.zeros((N,M+1))
    solution[:,0] = g(xInterior)
    uold = g(xInterior)

    amu = a*(dt/dx)
    print(amu)
    
    
    for c in range (M):

        unew =LaxWen(uold, amu)
        solution[:,c+1] = unew
        uold = unew

    #solution = np.vstack((solution, solution[0,:]))

    l2norm = np.zeros((M+1))
    for i in range (M+1):
        l2norm[i]= math.sqrt(dx)*lin.norm(solution[:,i])

    
    return solution, X, T, l2norm

def convdif(uold, a, d, dt, dx):
    N = len(uold)
    M = 1/dt
    dx2 = dx**2
    
    b = d/dx2 + a/(2*dx)
    c = - 2*d/dx2
    e = d/dx2 - a/(2*dx)
    A =  np.diag(c*np.full(N,1))+ np.diag(e*np.ones(N-1) ,1) + np.diag(b*np.ones(N-1),-1)
    A[0,N-1] = e
    A[N-1,0] = b
    unew = TRstep(A, uold, dt)

    return unew

def convdifInt(g, a, d, tstart, tend, M, N):

    dx =1/N
    xx = np.linspace(0, 1, N+1)
    xInterior = xx[:-1]

    dt = (tend-tstart)/M
    tt = np.linspace(tstart,tend,M+1)

    [T,X] = np.meshgrid(tt, xx)

    solution = np.zeros((N,M+1))
    solution[:,0] = g(xInterior)
    uold = g(xInterior)

    
    for c in range (M):

        unew = convdif(uold, a, d, dt, dx)
        solution[:,c+1] = unew
        uold = unew

    solution = np.vstack((solution, solution[0,:]))

    return solution, X, T

def LW(u, dt):
    N = u.size
    dx = 1/(N+1)
    unew = np.zeros(N)
    for i in range(N):
        ux = (u[i-1] - u[(i+1) % N]) / (2*dx)
        uxx = (u[i-1] - 2*u[i] + u[(i+1) % N]) / dx**2
        unew[i] = u[i] - dt*u[i]*ux + dt**2/2*(2 * u[i] * ux**2 + u[i]**2 * uxx)
    return unew
    

def TRLW(Tdx, uold, dt, d):
    T = d * dt/2 * Tdx
    I = np.identity(uold.size)
    return lin.solve(I - T, LW(uold, dt) + T @ uold)

    

def LWInt(g, d, tstart, tend, M, N):

    N = 250
    M = 1000
    d = 0.01
    tend = 1
    dt = tend/(M+1)
    dx = 1 / N
    xgrid = np.linspace(0, 1, N) 

    xx = np.linspace(0, 1, N)
    tt = np.linspace(0, tend, M + 1)
    T, X = np.meshgrid(tt, xx)

    Tdx = diags([1, -2, 1], [-1, 0, 1], shape=(N, N)).toarray()
    Tdx[0][-1] = 1
    Tdx[-1][0] = 1

    U = np.zeros((M+1, N))
    u = g(xgrid)
    for i in range(M + 1):
        U[i] = u
        u = TRLW(1/dx**2 * Tdx, u, dt, d)
        print(i)

    return U, X, T
