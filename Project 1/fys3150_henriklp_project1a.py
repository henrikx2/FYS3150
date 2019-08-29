from numpy import *
from matplotlib.pyplot import *
import os
import sys

def f(x):
    return 100*exp(-10*x)

def analyticSolution(x):
    return 1-(1-exp(-10))*x-exp(-10*x)

def generalSolution(a,b,c,d,n):
    a1, b1, c1, d1 = map(array,(a,b,c,d))

    for i in range(1,n): #Forward substitution on the interval [1,ng-1]
        b1[i] = b1[i] - c1[i-1]*(a1[i-1]/b1[i-1])
        d1[i] = d1[i] - d1[i-1]*(a1[i-1]/b1[i-1])

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
n = 10 # Choose the number of iterations
h = 1/(n+2) #Step size
x = linspace(0,1,n+2) #x-array
d = analyticSolution(x)

a_File = open(os.path.join(sys.path[0], "a"+str(n)),"r") #Opens file of context "an"
b_File = open(os.path.join(sys.path[0], "b"+str(n)),"r")
c_File = open(os.path.join(sys.path[0], "c"+str(n)),"r")
a_arr = zeros(n)
b_arr = zeros(n)
c_arr = zeros(n)

for i in range(0,n):
    '''
     a_arr[i] = eval(a_File.read().split(",")[i])
     b_arr[i] = eval(b_File.read().split(",")[i])
     b_arr[i] = eval(b_File.read().split(",")[i])
     '''
     print(a_File.read().split(","))
     print(type(a_File.read().split(",")))

u = generalSolution(a_arr,b_arr,c_arr,d,n)
'''
plot(x,h**2*u, label="General Algorithm")
legend()
xlabel("x")
ylabel("u")
'''
