from numpy import *
from matplotlib.pyplot import *

def f(x):
    return 100*exp(-10*x)

def analyticSolution(x):
    return 1-(1-exp(-10))*x-exp(-10x)

def generalSolution(a,b,c,d):
    for i in range(2,n):
        b[i] = b[i] - c[i-1]*(a[i]/b[i-1])
        d[i] = d[i] - d[i-1]*(a[i]/b[i-1])
        a[i] = a[i] - b[i-1]*(a[i]/b[i-1])

    v[n] = s[n]/b[n]
    for i in range(n-1, 1):
        


def specalSolution():

def LUdecomp():

n = [10,100,1000]

for i in n:
    h = 1/(n+1)
    x = linspace(0,1,n)
