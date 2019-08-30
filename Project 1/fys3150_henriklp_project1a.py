from numpy import *
from matplotlib.pyplot import *
import os
import sys

def f(x):
    return 100*exp(-10*x)

def analyticSolution(x):
    return 1-(1-exp(-10))*x-exp(-10*x)

def generalSolution(a,d,c,g,n):
    d_tilde = zeros(n)
    d_tilde[0] = d[0]
    g_tilde = zeros(n)
    for i in range(1,n): #Forward substitution on the interval [1,ng-1]
        d_tilde[i] = d[i] - c[i-1]*(a[i-1]/d_tilde[i-1])
        g_tilde[i] = g[i] - g_tilde[i-1]*(a[i-1]/d[i-1])

    u = a
    u[n-1] = g_tilde[n-1]/d_tilde[n-1]
    for i in range(n-2,-1,-1): #Backwards substitution on the interval [n-1,0]
        u[i] = (g[i] - c[i]*u[i+1])/d_tilde[i]
    return u

def specalSolution():
    return 0
def LUdecomp():
    return 0

#Analytic Solution

#General Solution
n = 100 # Choose the number of iterations
h = 1/(n) #Step size
x = linspace(0,1,n+2) #x-array
g_arr = analyticSolution(x)

a_arr = zeros(n)
a_arr += -1
d_arr = zeros(n)
d_arr += 2
c_arr = zeros(n)
c_arr += -1

'''
for i in n: #Use for all graphs in same plot
    u_arr = append([0],append(generalSolution(a_arr,d_arr,c_arr,g_arr,i),[0]))
    plot(x,u_arr, label="n = "+str(i))
'''

u_arr = append([0],append(generalSolution(a_arr,d_arr,c_arr,g_arr,n),[0])) #Just one plot
plot(x,u_arr, label="General Algorithm")

legend()
xlabel("x")
ylabel("u")
show()
