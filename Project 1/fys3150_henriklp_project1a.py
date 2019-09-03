from numpy import *
from matplotlib.pyplot import *
import os
import sys
import time
from scipy.linalg import lu_solve, lu_factor
print("This program solves a linear second order-differential equation by using a tridiagonal matrix. The calculation speeds and relative errors will be determined in the case of a general algorithm, a special algorithm and by the use of LU-decomposition.")
n = [pow(10,i) for i in range(1,int(input("Choose exponent of 10 (mesh grid): "))+1)] # Choose the number of iterations from power of 10

def f(x):
    return 100*exp(-10*x)

def analyticSolution(x):
    return 1-(1-exp(-10))*x-exp(-10*x)

def TDMAGeneralSolution(a,d,c,b,n):
    start_time = time.time()
    for i in range(1,n): #Forward substitution on the interval [1,n]
        ad = a[i-1]/d[i-1]
        d[i] = d[i] - c[i-1]*(ad)
        b[i+1] = b[i+1] - b[i]*(ad)

    u = zeros(len(b))
    u[n] = b[n]/d[n-1]
    for i in range(n-2,-1,-1): #Backwards substitution on the interval [n-1,1]
        u[i+1] = (b[i+1] - c[i]*u[i+2])/d[i]
    elapsed_time = time.time()-start_time
    print("[TDMA General (n = "+str(n)+")], calculation time: "+str(elapsed_time)+" s")
    return u

def TDMASpecalSolution(n):
    start_time = time.time()
    for i in range(1,n): #Forward substitution on the interval [1,n]
        ad = -1/d[i-1]
        d[i] = d[i] + ad
        b[i+1] = b[i+1] + ad

    u = zeros(len(b))
    u[n] = b[n]/d[n-1]
    for i in range(n-2,-1,-1): #Backwards substitution on the interval [n-1,1]
        u[i+1] = (b[i+1] + u[i+2])/d[i]
    elapsed_time = time.time()-start_time
    print("[TDMA Special (n = "+str(n)+")], calculation time: "+str(elapsed_time)+" s")
    return u

def LUdecomp(n):
    h = 1/n #Step size
    n = n-1 #Exclude endpoints
    A = zeros((n,n))
    b = zeros(n)
    x = zeros(n)
    A[0,0] = 2; A[0,1] = -1; #Set first row of A
    x = linspace(h,1-h,n); b[0] = h**2*f(x[0])
    x[n-1] = n*h; b[n-1] = h**2*f(x[n-1])
    for i in range(1,n-1): #Set middle of A
        x[i] = x[i] + h
        b[i] = h**2*f(x[i])
        A[i,i-1] = -1 #Lower diagona
        A[i,i] = 2 #Main diagonal
        A[i,i+1] = -1 #Upper diagonal
    A[n-1,n-1] = 2; A[n-1,n-2] = -1 #Set last row of A
    lu, piv = lu_factor(A)
    start_time = float(time.time())
    u = lu_solve((lu,piv),b)
    elapsed_time = time.time() - start_time
    print("[LU (n = "+str(n+1)+")], calculation time: "+str(elapsed_time)+" s")
    plot(x,u,label="LU-Decomp, n = "+ str(n+1))
    return u

#General Solution

for i in n:
    h = 1/(i+2) #Step size
    x = linspace(0,1,i+2) #x-array
    b_arr = h**2*f(x) #Source term
    a_arr = zeros(i)
    a_arr += -1 #Numbers on lower diagonal
    d_arr = zeros(i)
    d_arr += 2 #Numbers on main diagonal
    c_arr = zeros(i)
    c_arr += -1 #Numbers on upper diagonal
    u_arr = TDMAGeneralSolution(a_arr,d_arr,c_arr,b_arr,i)
    plot(x,u_arr, label = "n = " + str(i))
    if i == n[-1]:
        plot(x,analyticSolution(x), label="Analytical solution")

#Special solution

#LU-decomposition

for i in n:
    u = LUdecomp(i)
    if i == n[-1]:
        x = linspace(0,1,i)
        plot(x,analyticSolution(x), label="Analytical solution")

#Show plot
'''
legend()
xlabel("x")
ylabel("u")
show()
'''
