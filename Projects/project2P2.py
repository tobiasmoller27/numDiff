import numpy as np
import scipy.linalg as lin
import matplotlib.pyplot as plt
import math

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
plt.plot(np.linspace(0,1, 499),eigvectorTest, 'b+')
plt.plot(np.linspace(0,1, 499),np.max(eigvectorTest)*np.sin(np.linspace(0,1, 499)*5*(math.pi/2)), 'r')

plt.show()
"""