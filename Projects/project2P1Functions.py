import numpy as np
import numpy.linalg as lin

def twopBVP(fvec, alpha, beta, L, N):
    # CREATING THE RIGHT SIDE OF THE EQUATION
    Dx = L/(N+1)
    HL = np.zeros(N)
    for i in range(N):
        HL[i] = fvec[i]
    HL[0] = HL[0] - alpha/((Dx)**2)
    HL[N-1] = HL[N-1] - beta/((Dx)**2)
    HL = HL[np.newaxis]
    HL = HL.T
    
    # MULTIPLYING THE RIGHT SIDE WITH THE INVERSE OF T
    T = np.diag(np.full(N,-2))+np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1)
    Y = np.matmul(lin.inv(T),HL)*((Dx)**2)
    Y = np.vstack([[alpha],Y, [beta]])
    return Y

def statSchrodSolve(V, L, N):
    Dx = L/(N+1)
    xgrid = np.linspace(0,L,N)
    vvec = V(xgrid)

    T = (np.diag(np.full(N,-2))+np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1))/(Dx**2)
    # Add Dx^2*V to the diagonal 
    T = T + np.diagflat(vvec)
    eg, egfRow = lin.eig(T)
    egf = [[] for _ in range(len(T))]
    for i in range(len(T)):
        for v in egfRow:
            egf[i].append(v[i])
    eg = list(eg)
    egf = list(egf)
    sorted_egf = [x for _, x in sorted(zip(eg, egf))]
    sorted_eg = sorted(eg)
    return Dx, sorted_eg, sorted_egf, xgrid