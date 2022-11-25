import numpy as np
import scipy.linalg as lin
import scipy.integrate as intg
import matplotlib.pyplot as plt
import math
import project2P1Functions as lok
"""
NStart = 5
bigN = 499
NLambda = np.zeros((bigN-NStart,3))

NList = np.arange(NStart, bigN)
#NList = NList.reshape((bigN-NStart,1))


for N in NList:
    L = 1
    Dx = L/(N)
    TDx =  (np.diag(np.full(N+1,-2))+np.diag(np.ones(N),1)+np.diag(np.ones(N),-1))/(Dx**2)
    TDx[N][N] = -1/(Dx**2)
    TDx[N][N-1] = 1/(Dx**2)
    lambdas, eigvectors = lin.eigh(TDx)


    NLambda[N-NStart][0] = lambdas[-1] + math.pi**2/4
    NLambda[N-NStart][1] = lambdas[-2] + (3*math.pi/2)**2
    NLambda[N-NStart][2] = lambdas[-3] + (5*math.pi/2)**2
    eigvectorTest = eigvectors[-3]

NLambda1 = [row[0] for row in NLambda]
NLambda2 = [row[1] for row in NLambda]
NLambda3 = [row[2] for row in NLambda]

plt.loglog(NList, NLambda1,  label = "Lambda1")
plt.loglog(NList, NLambda2, label = "Lambda2")
plt.loglog(NList, NLambda3, label = "Lambda3")
plt.loglog(NList, 1/NList, label = "Referens")
plt.grid(True)
plt.legend()
plt.show()
"""
"""
plt.plot(np.linspace(0,1, 499),eigvectorTest, 'b+')
plt.plot(np.linspace(0,1, 499),np.max(eigvectorTest)*np.sin(np.linspace(0,1, 499)*5*(math.pi/2)), 'r')

plt.show()
"""

def V(x):
    return 80000*x

N = 100
Dx, eg, egf, xgrid = lok.statSchrodSolve(V, 1, N)


plt.plot(xgrid, egf[-3],"+")
plt.show()

# Norm
for i in range(-4,0):
    egfunc = [x**2 for x in (egf[i][:])]
    egint = intg.cumtrapz(egfunc, xgrid, initial=0)
    norm = 1/egint[-1]
    normegfunc = [x*norm for x in egfunc]
    plt.plot(xgrid, normegfunc)
plt.show()

