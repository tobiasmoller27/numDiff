import project0Functions as lok
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as linalg


#TASK 1.2 
N = 100
[tgrid, approx, err] = lok.eulerint(np.matrix('-1 10; 0 -3'), np.matrix('1;1'), 0, 1, N)
#[tgrid, approx,err] = lok.eulerint(np.matrix('3'),np.matrix('1'),0,1,20)
print("Error: " + str(err))
print(len(tgrid))
realY = [0]*tgrid
for i in range(len(tgrid)):
    A = linalg.expm(tgrid[i]*np.matrix('-1 10; 0 -3'))*np.matrix('1;1')
    realY[i] = A[0][0]
plt.title("Task 1.2")
#plt.plot(tgrid, approx, label = "Approximation")
plt.plot(tgrid, realY, label = "Real function")
plt.plot(tgrid, err, label = "Error")
plt.plot(tgrid, approx[0][:], label = "Approximation")
plt.legend()
plt.show()


"""
#TASK 1.3

h1, e1 = lok.errVSh(np.matrix(1), 1, 0, 1)
h2, e2 = lok.errVSh(np.matrix(2), 1, 0, 1)
h3, e3 = lok.errVSh(np.matrix(3), 1, 0, 1)
h4, e4 = lok.errVSh(np.matrix(4), 1, 0, 1)
h5, e5 = lok.errVSh(np.matrix(5), 1, 0, 1)

plt.title("Error per h")
plt.loglog(h1, e1, label ="A = 1")
plt.loglog(h2, e2, label ="A = 2")
plt.loglog(h3, e3, label ="A = 3")
plt.loglog(h4, e4, label ="A = 4")
plt.loglog(h5,e5, label = "A = 5")
plt.legend()
plt.show()

# Felet ökar avtagande 

"""
"""
# TASK 1.4
[tgrid, approx, err] = lok.eulerint(np.matrix(3), 1, 0, 1, 1000)
plt.title("Error per timepoint")
plt.xlabel("Time")
plt.ylabel("log(err)")
plt.semilogy(tgrid, err)
plt.show()
"""
"""
# TASK 1.5

A = np.matrix('-1 10; 0 -3')
B = np.matrix('1 3 9; 2 3 8; 0 12 7')
y0 = np.matrix('1;1')
h1, e1 = lok.errVSh(A, y0, 0, 10)
h2, e2 = lok.errVSh(A, y0, 0, 10)
#h3, e3 = lok.errVSh(np.matrix(3), 1, 0, 1)
#h4, e4 = lok.errVSh(np.matrix(4), 1, 0, 1)
#h5, e5 = lok.errVSh(np.matrix(5), 1, 0, 1)

plt.title("Error per h")
plt.loglog(h1, e1, label ="A = 1")
plt.loglog(h2, e2, label ="A = 2")
#plt.loglog(h3, e3, label ="A = 3")
#plt.loglog(h4, e4, label ="A = 4")
#plt.loglog(h5,e5, label = "A = 5")
plt.legend()
plt.show()
"""

#Task 1.6  Inexplicit euler
"""
N = 100
[tgrid, approx, err] = lok.ieulerint(np.matrix('-1 100; 0 -30'), np.matrix('1;1'), 0, 1, N)
print("err: ")
print(err)





plt.title("Task 1.2")
plt.plot(tgrid, approx, label = "Approximation")
plt.plot(tgrid, realY, label = "Real function")
plt.legend()
plt.show()

"""
"""

A = np.matrix('-1 10; 0 -3')
#B = np.matrix('1 3 9; 2 3 8; 0 12 7')
y0 = np.matrix('1;1')
h1, e1 = lok.errVSh(A, y0, 0, 10)
#h2, e2 = lok.errVSh(A, y0, 0, 10)
#h3, e3 = lok.errVSh(np.matrix(3), 1, 0, 1)
#h4, e4 = lok.errVSh(np.matrix(4), 1, 0, 1)
#h5, e5 = lok.errVSh(np.matrix(5), 1, 0, 1)
print("e1 size: " + str(len(e1)) + "h1 size: " + str(len(h1)))
print(e1)
print(h1)
plt.title("Error per h")

plt.loglog(h1, e1, label ="A")
#plt.loglog(h2, e2, label ="A = 2")
#plt.loglog(h3, e3, label ="A = 3")
#plt.loglog(h4, e4, label ="A = 4")
#plt.loglog(h5,e5, label = "A = 5")
plt.legend()
plt.show()
"""
"""
A = np.matrix('-1 10; 0 -3')
B = np.matrix('1 3 9; 2 3 8; 0 12 7')
y0 = np.matrix('1;1')
h1, e1 = lok.ierrVSh(A, y0, 0, 10)
#h2, e2 = lok.ierrVSh(A, y0, 0, 10)
#h3, e3 = lok.errVSh(np.matrix(3), 1, 0, 1)
#h4, e4 = lok.errVSh(np.matrix(4), 1, 0, 1)
#h5, e5 = lok.errVSh(np.matrix(5), 1, 0, 1)

plt.title("Error per h")
plt.loglog(h1, e1, label ="A = 1")
#plt.loglog(h2, e2, label ="A = 2")
#plt.loglog(h3, e3, label ="A = 3")
#plt.loglog(h4, e4, label ="A = 4")
#plt.loglog(h5,e5, label = "A = 5")
plt.legend()
plt.show()
"""
"""
N = 100
[tgrid, approx, err] = lok.eulerint(np.matrix('3'), np.matrix('1'), 0, 1, N)
print("tgrid: " + str(tgrid))

realY = []


for i in range(len(tgrid)):
    realY.append(linalg.expm(tgrid[i]*np.matrix('3'))*np.matrix('1'))
#upperApprox = approx[0, :]
#lowerApprox = approx[1, :]

plt.title("Task 1.2")
#plt.plot(tgrid, approx, label = "Approximation")
#plt.plot(tgrid, realY, label = "Real function")
plt.plot(tgrid, approx)
plt.plot(tgrid, realY)


plt.legend()
plt.show()
"""