from numpy import *
from matplotlib.pyplot import *
import os
import sys

def f(x):
    return 100*exp(-10*x)

def analyticSolution(x):
    return 1-(1-exp(-10))*x-exp(-10*x)

def generalSolution(a,d,c,b,n):
    for i in range(1,n): #Forward substitution on the interval [1,ng-1]
        ad = a[i-1]/d[i-1]
        d[i] = d[i] - c[i-1]*(ad)
        b[i+1] = b[i+1] - b[i]*(ad)

    u = zeros(len(b))
    u[n] = b[n]/d[n-1]
    for i in range(n-2,-1,-1): #Backwards substitution on the interval [n-1,0]
        u[i+1] = (b[i+1] - c[i]*u[i+2])/d[i]
    return u

def specalSolution():
    return 0
def LUdecomp():
    return 0

#General Solution
n = [pow(10,i) for i in range(1,int(input("Type power of 10: "))+1)] # Choose the number of iterations from power of 10
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
    u_arr = generalSolution(a_arr,d_arr,c_arr,b_arr,i)
    plot(x,u_arr, label = "n = " + str(i))
    if i == n[-1]:
        plot(x,analyticSolution(x), label="Analytical solution")
legend()
xlabel("x")
ylabel("u")
show()

#Special solution

#LU-decomposition
