import numpy as np
import scipy.linalg as lin
L = 1
N = 499
Dx = L/(N+1)
def p(x, lmbd):
    return lmbd*x
TDx =  (np.diag(np.full(N,-2))+np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1))/(Dx**2)
print(TDx)
lambdas, eigvectors = lin.eigh(TDx)
#idx = np.argsort(lambdas)
#lambdas = lambdas[idx]
#eigvectors = eigvectors[:, idx]

print(lambdas[-1])
print(lambdas[-2])
print(lambdas[-3])


"""
z
xInterior = np.linspace(0+Dx, L-Dx, N)
M = lok.twopBVP(q,0,0,L,N)
x = np.linspace(0, L, N+2)"""