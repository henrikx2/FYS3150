from numpy import *
from matplotlib.pyplot import *
import os
import sys
import time
from scipy.linalg import solve, lu
print("This program solves a linear second order-differential equation by using a tridiagonal matrix. The calculation speeds and relative errors will be determined in the case of a general algorithm, a special algorithm and by the use of LU-decomposition.")
n = [pow(10,i) for i in range(1,int(input("Choose exponent of 10 (# of mesh points): "))+1)] # Choose the number of iterations from power of 10

def f(x):
    return 100*exp(-10*x)

def analyticSolution(x):
    return 1-(1-exp(-10))*x-exp(-10*x)

def TDMAGeneral(a,d,c,b,v,n):
    start_time = time.perf_counter()
    for i in range(1,n): #Forward substitution on the interval [1,n]
        ad = a[i-1]/d[i-1]
        d[i] = d[i] - c[i-1]*(ad)
        b[i+1] = b[i+1] - b[i]*(ad)

    v[n] = b[n]/d[n-1]
    for i in range(n-2,-1,-1): #Backwards substitution on the interval [n-1,1]
        v[i+1] = (b[i+1] - c[i]*v[i+2])/d[i]
    elapsed_time = time.perf_counter()-start_time
    print("(n = "+str(n)+")[TDMA General], CPU Time: "+str(elapsed_time))
    return v

def TDMASpecial(d,b,v,n):
    start_time = time.perf_counter()
    for i in range(1,n):
        ad = -1/d[i-1]
        d[i] = d[i] + ad
        b[i+1] = b[i+1] - b[i]*ad

    v[n] = b[n]/d[n-1]
    for i in range(n-2,-1,-1):
        v[i+1] = (b[i+1] + v[i+2])/d[i]
    elapsed_time = time.perf_counter() - start_time
    print("(n = "+str(n)+")[TDMA Special], CPU Time: "+str(elapsed_time))
    return v

def LUdecomp(n):
    h = 1/(n+1) #Step size
    A = zeros((n,n))
    b = zeros(n)
    x = zeros(n)
    A[0,0] = 2; A[0,1] = -1; #Set first row of A
    x = linspace(h,1-h,n); b[0] = h**2*f(x[0])
    x[n-1] = n*h; b[n-1] = h**2*f(x[n-1])
    for i in range(1,n-1): #Set middle of A
        x[i] = x[i] + h
        b[i] = h**2*f(x[i])
        A[i,i-1] = -1 #Lower diagonal
        A[i,i] = 2 #Main diagonal
        A[i,i+1] = -1 #Upper diagonal
    A[n-1,n-1] = 2; A[n-1,n-2] = -1 #Set last row of A
    L, U, piv = lu(A)
    start_time = float(time.perf_counter())
    z = solve(U,b)
    v_lu = solve(L,z)
    elapsed_time = time.perf_counter() - start_time
    print("(n = "+str(n)+")[LU-Decomp], CPU Time: "+str(elapsed_time))
    return v_lu

for i in n:
    h = 1/(i+2) #Step size
    x = linspace(0,1,i+2) #x-array
    exact_arr = analyticSolution(x)

    b_arr = h**2*f(x) #Source term
    a_arr = zeros(i)
    a_arr += -1 #Numbers on lower diagonal
    d_arr = zeros(i)
    d_arr += 2 #Numbers on main diagonal
    c_arr = zeros(i)
    c_arr += -1 #Numbers on upper diagonal
    v_g = zeros(len(b_arr))
    v_s = zeros(len(b_arr))

    #General Solution
    v_general = TDMAGeneral(a_arr,d_arr,c_arr,b_arr,v_g,i)
    relative_error_general = abs(exact_arr[1:-2] - (v_general[1:-2]/exact_arr[1:-2]))

#Special solution
    v_special = TDMASpecial(d_arr,b_arr,v_s,i)
    relative_error_special = abs(exact_arr[1:-2] - (v_special[1:-2]/exact_arr[1:-2]))

#LU-decomposition
    v_lu = append([0],append(LUdecomp(i),[0]))
    relative_error_lu = abs(exact_arr[1:-2] - (v_lu[1:-2]/exact_arr[1:-2]))

#Plots
    #SOLUTIONS
    plot(x,v_general, label = "n = " + str(i))
    #plot(x,v_special, label = "n = " + str(i))
    #plot(x,v_lu,label="LU-Decomp, n = "+ str(i))

    if i == n[-1]:
        plot(x,exact_arr, label="Analytical solution")

    #RELATIVE ERROR
    #plot(x[1:-2],relative_error_general[1:-2], label = "n = " + str(i))
    #plot(x,relative_error_special, label = "n = " + str(i))
    #plot(x,relative_error_lu, label = "n = " + str(i))
    title("Relative error")
    xlabel("x")
    ylabel("Log($\epsilon$)")


#Show plot

legend()
show()
