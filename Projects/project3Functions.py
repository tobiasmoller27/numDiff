import numpy as np


def eulerstep(Tdx, uold, dt):
    return np.reshape(uold + dt*np.matmul(Tdx,uold), uold.shape)


def eulerint(g, tstart, tend, M, N):
    #Defining timeinterval and timegrid
    dt=(tend-tstart)/M
    tt = np.linspace(0, tend, M+1)
    
    #Defining x-grid
    dx = 1/(N+1)
    xx = np.linspace(0, 1, N+2)
    xInterior = xx[1:-1]

    #Constructing the meshgrid
    [X,T] = np.meshgrid(xx, tt)
    #Constructing the Toeplitzmatrix as in Project 2
    Tdx = (np.diag(np.full(N,-2))+np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1))/(dx**2)    

    #Constructing solution matrix

    solution = np.zeros((N,M+1))
    solution[:,0] = g(xInterior)

    uold = g(xInterior)
    for c in range (len(solution)):
        unew =eulerstep(Tdx, uold, dt)
        solution[:,c+1] = unew
        uold = unew

    solution = np.hstack((np.zeros((M,1)),solution,np.zeros((M,1))))
   
    print(xx.shape)
    print(xInterior.shape)
    return solution, X, T
    


