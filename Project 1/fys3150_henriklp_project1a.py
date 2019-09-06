from numpy import *
from matplotlib.pyplot import *
import sys
import time
from scipy.linalg import solve, lu

#User interface Q&A
print("This program solves a linear second order-differential equation by using a tridiagonal matrix derived from an approximation of the 2nd derivative. The calculation speeds and relative error will be determined in the case of a general algorithm, a special algorithm and by the use of LU-decomposition. The relative error will be shown as a function of mesh points.")
n = [pow(10,i) for i in range(1,int(input("Choose exponent of 10 to define numbers of mesh points: "))+1)] # Choose the power of 10 for the number of iterations
g_ask = input("Calculate with General algorithm?\n(y/n): ")
s_ask = input("Calculate with Special algorithm?\n(y/n): ")
lu_ask = input("Calculate with LU-Decomposition?\n(y/n): ")
if g_ask != "y":
    print("Skipping General algorithm.")
if s_ask != "y":
    print("Skipping Special algorithm.")
if lu_ask != "y":
    print("Skipping LU-Decomposition.")

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
    for i in range(n-2,-1,-1): #Backwards substitution
        v[i+1] = (b[i+1] - c[i]*v[i+2])/d[i]
    elapsed_time = time.perf_counter()-start_time
    print("(n = "+str(n)+")[TDMA General], CPU Time: "+str(elapsed_time))
    return v

def TDMASpecial(d,b,v,n):
    start_time = time.perf_counter()
    for i in range(1,n):
        b[i+1] = b[i+1] + b[i]/d[i-1]

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
    v = solve(L,z)
    elapsed_time = time.perf_counter() - start_time
    print("(n = "+str(n)+")[LU-Decomp], CPU Time: "+str(elapsed_time))
    return v

relative_error = zeros(len(n)) #Initialize relative error arroy

for i in n:
    h = 1/(i+2) #Step size
    x = linspace(0,1,i+2) #x-array
    exact_arr = analyticSolution(x) #Analytical solution

    #General Solution
    v_g = zeros(i+2)
    if g_ask == "y":
        b_g = h**2*f(x) #Source term
        a_g = zeros(i)
        a_g += -1 #Numbers on lower diagonal
        d_g = zeros(i)
        d_g += 2 #Numbers on main diagonal
        c_g = zeros(i)
        c_g += -1 #Numbers on upper diagonal
        v_g = TDMAGeneral(a_g,d_g,c_g,b_g,v_g,i) #Activation general algorithm

    #Special solution
    v_s = zeros(i+2)
    if s_ask == "y":
        b_s = h**2*f(x) #Source term
        d_s = zeros(i)
        d_s[0] = 2 #First number in diagonal-array
        for k in range(1,i): #Setting rest of diagonal-array before calculation
            d_s[k] = (k+2)/(k+1)
        v_s = TDMASpecial(d_s,b_s,v_s,i) #Activation special algorithm


    #LU-decomposition
    v_lu =zeros(i)
    if lu_ask == "y" and i <= pow(10,4):
        v_lu = append([0],append(LUdecomp(i),[0])) #Activation LU-Decomposition
    if i > pow(10,4): #Trying LU-Decomposition with more than 10000x10000 matrix elements
        try:
            v_lu = append([0],append(LUdecomp(i),[0])) #Activation LU-Decomposition

        except:
            print("Unable to allocate array with shape ("+str(i)+", "+str(i)+").")


    #Relative error
    relative_error[int(log10(i))-1] = max((exact_arr[1:-1] - v_s[1:-1])/exact_arr[1:-1])

#Plots
    #SOLUTIONS
    #plot(x,v_g, label = "n = " + str(i)))
    #plot(x,v_s, label = "n = " + str(i))
    #plot(x,v_lu,label="LU-Decomp, n = "+ str(i))
    '''
    if i == n[-1]:
        plot(x,exact_arr, label="Analytical solution")
    '''
    #title("Relative error")
    #xlabel("x")
    #ylabel("Log($\epsilon$)")

#Plots

plot_error_ask = input("Plot relative error?\n(y/n): ")
if plot_error_ask == "y":
    plot(log10(n),relative_error)
    title("Relative error of $Log_{10}$(n)")
    xlabel("$Log_{10}$(n)")
    ylabel("\epsilon(n)")
    show()
