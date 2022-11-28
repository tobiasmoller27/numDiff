import numpy as np
import scipy.linalg as lin
import scipy.integrate as intg
import matplotlib.pyplot as plt
import math
import project2P1Functions as lok
"""
NStart = 5
bigN = 500
NLambda = np.zeros((bigN-NStart,3))

NList = np.arange(NStart, bigN)
#NList = NList.reshape((bigN-NStart,1))


for N in NList:
    L = 1
    Dx = L/(N)
    TDx =  (np.diag(np.full(N,-2))+np.diag(np.ones(N-1)*(1),1)+np.diag(np.ones(N-1)*(1),-1))/(Dx**2)
    TDx[N-1][N-2] = 2/(Dx**2)
    
    if (N==10):
        print(TDx)
    
    lambdas, eigvectors = lin.eig(TDx)
    # Sort the arrays
    lambdas = np.sort(lambdas)
    if N==499:
        print(lambdas[-1])
        print(lambdas[-2])
        print(lambdas[-3])
    

    

    NLambda[N-NStart][0] = abs(lambdas[-1] + (math.pi**2)/4)
    NLambda[N-NStart][1] = abs(lambdas[-2] + (3*math.pi/2)**2)
    NLambda[N-NStart][2] = abs(lambdas[-3] + (5*math.pi/2)**2)

NLambda1 = [row[0] for row in NLambda]
NLambda2 = [row[1] for row in NLambda]
NLambda3 = [row[2] for row in NLambda]

plt.loglog(NList, NLambda1,  label = "Lambda1")
plt.loglog(NList, NLambda2, label = "Lambda2")
plt.loglog(NList, NLambda3, label = "Lambda3")
plt.loglog(NList, 1/np.square(NList), label = "Referens")
plt.xlabel("Amount of steps")
plt.ylabel("Error")
plt.title("Error per Amount of Steps")
plt.grid(True)
plt.legend()
plt.show()
"""

"""
# Plotting first three eigenmodes for N = 499
N = 499
L = 1
Dx = L/(N)
TDx =  (np.diag(np.full(N+1,2))+np.diag(np.ones(N)*(-1),1)+np.diag(np.ones(N)*(-1),-1))/(Dx**2)
TDx[N][N-1] = -2/(Dx**2)

eg, egfRow = lin.eig(TDx)
# Sort the arrays
egf = [[] for _ in range(len(TDx))]
for i in range(len(TDx)):
    for v in egfRow:
        egf[i].append(v[i])
eg = list(eg)
egf = list(egf)
sorted_egf = [x for _, x in sorted(zip(eg, egf))]
sorted_eg = sorted(eg)

x = np.linspace(0,1, 500)
for i in range(3):
    plt.plot(x, sorted_egf[i], label="Lambda = "+str(round(np.real(sorted_eg[i]),2)))
plt.grid(True)
plt.xlabel("Gridpoints")
plt.ylabel("Approximation")
plt.title("Approximation of the three first eigenmodes for N=499")
plt.legend()
plt.show()
"""

# Task 2.2 

def V(x):
    #return 700*(0.5 - abs(x-0.5))
    return 800*np.sin(math.pi*x)**2
    #return 0

N = 1000
Dx, eg, egf, xgrid = lok.statSchrodSolve(V, 1, N)

for i in range(1,7):
    plt.plot(xgrid, [300*x for x in egf[-i]]+eg[-i], label="Lambda =  "+ str(round(np.real(eg[-i]),2))) #OBS ANNAN AMPLITUD!!!!
plt.xlabel("Gridpoints")
plt.ylabel("Eigenfunctions")
plt.legend()
plt.show()
# Norm
for i in range(1,7):
    egfSquared = [x**2 for x in egf[-i]]
    egint = abs(intg.cumtrapz(egfSquared, xgrid, initial=0))
    norm = math.sqrt(1/egint[-1])
    possibilityNormed = [norm*x for x in egfSquared]
    plt.plot(xgrid, possibilityNormed,label="Lambda =  "+ str(round(np.real(eg[-i]),2)))
plt.xlabel("Position")
plt.ylabel("Probability density")
plt.legend()
plt.show()

