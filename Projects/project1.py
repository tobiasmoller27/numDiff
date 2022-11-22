import project1Functions as lok
import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.linalg as linalg
import scipy.integrate as scint


"""
# Test of RK4walk
def p(told, uold):
    #return np.matrix('3 9;15 15')*uold
    return 3*uold
"""
"""
[tgrid, approx, localErr] = lok.RK4walk(p,0,1,1,10)
realY = np.exp(3*tgrid)


plt.title("Aprroximated and real curve")
plt.plot()
plt.plot(tgrid,approx, label = "Approx")
plt.plot(tgrid, realY, label = "Real")
plt.legend()

glb = abs(realY - approx)

plt.subplot(1,2,2)
plt.title("Global error in log-scale")
plt.loglog(tgrid, glb, label ="Global error")
plt.legend()

plt.show()
print(glb[-1])

"""
"""
# Task 1.2
[unew, rnew] = lok.RK34step(p,0,1,0.01)
print(unew)
print(rnew)

"""
"""
# Task 1.4
t0 = 0
y0 = 1
tf = 1
tol  = 10**(-6)
[t,y] = lok.adaptiveRK34(p, t0, tf, y0, tol)
plt.plot(t,y)
print(y[-1])
print(t[-1])
plt.show()
"""

"""
t0 = 0
y0 = np.matrix('1;1')
tf = 1
tol  = 10**(-6)
[t,y] = lok.adaptiveRK34m(p, t0, tf, y0, tol)
yClean = []
yClean2 = []
for i in range(len(t)):
    yClean.append(y[0,i])
    yClean2.append(y[1,i])
print("Real endvalues: ")
print(linalg.expm(np.matrix('3 9; 15 15'))*np.matrix('1;1'))
plt.plot(t,yClean, label="Övre")
plt.plot(t,yClean2, label="Undre")
plt.legend()
plt.show()

"""
"""
# Task 2.1

def lotka(t, u):
    a, b, c, d = 3, 9, 15, 15
    x = a*u[0] -b*u[0]*u[1]
    y = c*u[0]*u[1]-d*u[1]
    r = np.zeros((2,1))
    r[0] = x
    r[1] = y
    return r

def H(x,y):
    a, b, c, d = 3, 9, 15, 15
    return c*x + b*y - d*np.log(x)-a*np.log(y)


t0 = 0
tf = 100
tol = 10**(-6)
y0 = np.matrix('60; 70')
[t,y] = lok.adaptiveRK34m(lotka, t0, tf, y0, tol)
HCompare = abs(H(y[0,:], y[1,:])/H(y0[0],y0[1])-1)
yClean = []
yClean2 = []
hPlottable = []
for i in range(len(t)):
    yClean.append(y[0,i])
    yClean2.append(y[1,i])
    hPlottable.append(HCompare[0,i])

plt.plot(t,yClean, label="Rabbit")
plt.plot(t,yClean2, label="Fox")
#plt.plot(yClean2,yClean)
plt.legend()
plt.show()


# Plot H
#plt.semilogy(t,hPlottable)
#plt.show()
"""
"""
# TASK 3.1
def vdP(t,u):
    r = np.zeros((2,1))
    r[0] = u[1]
    r[1] = 200*(1-u[0]**2)*u[1]-u[0]
    return r

    
t0 = 0
tf = 400
tol = 10**(-6)
y0 = np.matrix('3; 0')
[t,y] = lok.adaptiveRK34m(vdP, t0, tf, y0, tol)
yClean = []
yClean2 = []
for i in range(len(t)):
    yClean.append(y[0,i])
    yClean2.append(y[1,i])


#plt.plot(t,yClean2, label="Fox")
plt.plot(yClean,yClean2)
plt.legend()
plt.show()
"""

# Task 3.2


"""
mus = [10,15,22,33,47,68,100,150,220,330,470,680]
stepsMu = []
for mu in mus:
    def vdPwithMu(t,u):
        r = np.zeros((2,1))
        r[0] = u[1]
        r[1] = mu*(1-u[0]**2)*u[1]-u[0]
        return r

    t0 = 0
    tf = 0.7*mu
    tol = 10**(-7)
    y0 = np.matrix('2; 0')
    [t,y] = lok.adaptiveRK34m(vdPwithMu, t0, tf, y0, tol)
    stepsMu.append(len(t))
    print(len(t))

plt.title("Steps per mu^2")
plt.loglog(np.square(mus),stepsMu)
plt.show()
"""

# Task 3.3

def vdP(t,u):
    r = [0,0]
    r[0] = u[1]
    r[1] = 100*(1-u[0]**2)*u[1]-u[0]
    return r

    
t0 = 0

tol = 10**(-6)
y0 = [2,0]
#resultIVP= scint.solve_ivp(vdP, (t0,tf), y0, method="BDF")

"""
print(resultIVP.t)
print(resultIVP.y)
plt.plot(resultIVP.t,resultIVP.y[0,:])
plt.plot(resultIVP.t, resultIVP.y[1,:])
plt.show()
"""


mus = [10,15,22,33,47,68,100,150,220,330,470,680,1000]
nSteps = []
for mu in mus:
    tf=0.7*mu
    
    def vdPwithMu(t,u):
        r = [0,0]
        r[0] = u[1]
        r[1] = mu*(1-u[0]**2)*u[1]-u[0]
        return r
    resultIVP= scint.solve_ivp(vdPwithMu, (t0,tf), y0, method="BDF")
    nSteps.append(len(resultIVP.t))
    print("Mu: " +str(mu) + "   Nr of steps: " +str(len(resultIVP.t)))
plt.loglog(mus,nSteps)
plt.show()


