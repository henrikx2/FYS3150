from numpy import *
from matplotlib.pyplot import *
import sys
from tabulate import tabulate

n = [pow(10,i) for i in range(1,int(input("Choose exponent of 10 to define numbers of mesh points: "))+1)] # Choose the power of 10 for the number of iterations

def f(x):
    return 100*exp(-10*x)

def analyticSolution(x):
    return 1-(1-exp(-10))*x-exp(-10*x)

def TDMASpecial(d,b,v,n):
    for i in range(1,n):
        b[i+1] = b[i+1] + b[i]/d[i-1]

    v[n] = b[n]/d[n-1]
    for i in range(n-2,-1,-1):
        v[i+1] = (b[i+1] + v[i+2])/d[i]
    return v

error_arr = zeros(len(n)) #Initialize relative error arroy

for i in n:
    h = 1/(i+2) #Step size
    x = linspace(0,1,i+2) #x-array
    exact_arr = analyticSolution(x) #Analytical solution

    v_s = zeros(i+2)
    b_s = h**2*f(x) #Source term
    d_s = zeros(i)
    d_s[0] = 2 #First number in diagonal-array
    for k in range(1,i): #Setting rest of diagonal-array before calculation
        d_s[k] = (k+2)/(k+1)
    v_s = TDMASpecial(d_s,b_s,v_s,i) #Activation special algorithm
    error = max((exact_arr[1:-1]-v_s[1:-1])/exact_arr[1:-1])
    error_arr[int(log10(i))-1] = error

table_list = [[str(points),str(err)] for points,err in zip(n,error_arr)]
print(tabulate(table_list, headers=["n","Relative error"]))

plot_error_ask = input("Plot relative error?\n(y/n): ")
if plot_error_ask == "y":
    plot(log10(n),error_arr)
    title("Relative error of $Log_{10}$(n)")
    xlabel("$Log_{10}$(n)")
    ylabel("$\epsilon(n)$")
    show()
