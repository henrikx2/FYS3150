from numpy import *
from matplotlib.pyplot import *
import os
import sys

def f(x):
    return 100*exp(-10*x)

def analyticSolution(x):
    return 1-(1-exp(-10))*x-exp(-10*x)

def generalSolution(a,b,c,d):
    a1, b1, c1, d1 = map(array,(a,b,c,d))
    ng = len(d)
    for i in range(1,ng): #Forward substitution on the interval [1,ng-1]
        w1 = zeros(ng-1)
        w1[i] = a1[i]/b1[i-1]
        b1[i] = b1[i] - c1[i-1]*(w1[i])
        d1[i] = d1[i] - d1[i-1]*(w1[i])
        a1[i] = a1[i] - b1[i-1]*(w1[i]) #Equals zero and can be skipped

    v1 = a1
    v1[n] = d1[n]/b1[n]
    for i in range(ng-2,-1,-1): #Backwards substitution on the interval [ng-2,0]
        v1[i] = (d1[i] - c1[i]*v1[i+1])/b1[i]

    return v1

def specalSolution():
    return 0
def LUdecomp():
    return 0

#Analytic Solution

#General Solution
n = 10 # Choose the approperiate n-value

a_arr = open(os.path.join(sys.path[0], "a"+str(n)), "r")
b_arr = open(os.path.join(sys.path[0], "b"+str(n)), "r")
c_arr = open(os.path.join(sys.path[0], "c"+str(n)), "r")

h = 1/(n+1)
x = linspace(0,1,n)
'''
u = generalSolution(a_arr,b_arr,c_arr,d)

plot(x,h**2*u, label="General Algorithm")
legend()
xlabel("x")
ylabel("u")
'''
