import numpy as np
import scipy.linalg as lin


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
    dx =1/N
    xx = np.linspace(0, 1, N+1)
    xInterior = xx[:-1]

    dt = (tend-tstart)/M
    tt = np.linspace(tstart,tend,M+1)

    [T,X] = np.meshgrid(tt, xx)

    solution = np.zeros((N,M+1))
    solution[:,0] = g(xInterior)
    uold = g(xInterior)

    amu = a*(dt/dx)
    
    l2norm = np.zeros((M+1, 1))
    l2norm[0] = lin.norm(g(xInterior))
    print(l2norm[0])
    
    for c in range (M):

        unew =LaxWen(uold, amu)
        solution[:,c+1] = unew
        l2norm[c+1] = lin.norm(unew)
        uold = unew

    solution = np.vstack((solution, solution[0,:]))

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

def LW(uold, N, dx, dt, Tdx):
    Sdx = (np.diag(np.full(N,0))+np.diag((-1)*np.ones(N-1),1)+np.diag(np.ones(N-1),-1))/(2*dx)
    #unew = uold - dt*np.matmul(uold, Sdx) + ((dt**2)/2)*(2*np.matmul(uold, np.matmul(Sdx, Sdx)) + np.matmul(np.square(uold), Tdx))
    unew = uold - dt*np.matmul(uold, Sdx) + ((dt**2)/2)*(2*np.matmul(uold, np.square(Sdx)) + np.matmul(np.square(uold), Tdx))
    return unew

def TRLW(uold, N, dx, dt, d):
    Tdx = (np.diag(np.full(N,-2))+np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1))/(dx**2)
    unew = np.matmul(lin.inv(np.eye(N)-(d*dt/2)*Tdx),(LW(uold, N, dx, dt, Tdx) + (d*dt/2)*np.matmul(Tdx, uold)))
    
    return unew

def LWInt(g, d, tstart, tend, M, N):

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

        unew =TRLW(uold, N, dx, dt, d)
        solution[:,c+1] = unew
        uold = unew

    solution = np.vstack((solution, solution[0,:]))

    return solution, X, T
    

