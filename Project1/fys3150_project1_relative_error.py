from numpy import *
from matplotlib.pyplot import *
import sys
import os
from tabulate import tabulate

#User interface Q&A
print("This program computes the relative error as a function of iterations for a linear 2nd order differental equation solved with the specialized algorithm.")
n = [pow(10,i) for i in range(1,int(input("Choose exponent of 10 to define numbers of mesh points: "))+1)] # Choose the power of 10 for the number of iterations

#Make directory to save data (if ot already created)
path = os.getcwd()+"\Data"
try:
    os.mkdir(path)
except OSError:
    print ("Creation of the directory %s failed. It may already exist." % path)
else:
    print ("Successfully created the directory %s " % path)

def f(x): #Source term
    return 100*exp(-10*x)

def analyticSolution(x): #Analytical solution
    return 1-(1-exp(-10))*x-exp(-10*x)

def TDMASpecial(d,b,v,n):
    for i in range(1,n): #Forward substitution
        b[i+1] = b[i+1] + b[i]/d[i-1]

    v[n] = b[n]/d[n-1]
    for i in range(n-2,-1,-1): #Backwards substitution
        v[i+1] = (b[i+1] + v[i+2])/d[i]
    return v

error_arr = zeros(len(n)) #Initialize relative error arroy
h_arr = zeros(len(n))

#Loop over different n-values
for i in n:
    h = 1/(i+1) #Step size
    x = linspace(0,1,i+2) #x-array for plotting
    exact_arr = analyticSolution(x) #Analytical solution

    v_s = zeros(i+2) #Initialize solution array
    b_s = h**2*f(x) #Source term
    d_s = zeros(i) #Iniialize diagonal array
    d_s[0] = 2 #First number in diagonal-array
    for k in range(1,i): #Setting rest of diagonal-array before calculation
        d_s[k] = (k+2)/(k+1)
    v_s = TDMASpecial(d_s,b_s,v_s,i) #Activation of special algorithm
    error = max((exact_arr[1:-1]-v_s[1:-1])/exact_arr[1:-1]) #Computes largest relative error
    error_arr[int(log10(i))-1] = error
    h_arr[int(log10(i))-1] = h

#Make a table with relative error
table_list = [[str(points),str(err)] for points,err in zip(log10(h_arr),log10(error_arr))]
error_table = tabulate(table_list, headers=["log10(h)","log10(Relative error)"])
#Save table to file
error_file = open(path+"\Relative_Error.txt", "w")
error_file.write(tabulate(error_table))
error_file.close()
print(error_table)
print("Table saved as Relative_Error.dat")

#Plot relative error
plot_error_ask = input("Plot relative error?\n(y/n): ")
if plot_error_ask == "y":
    loglog(h_arr,error_arr)
    title("Loglog plot of relative error as function of h.")
    xlabel("h (Step size)")
    ylabel("$\epsilon$ (Relative error)")
    savefig(path+"\Relative_Error_Plot.png")
    slope, i = polyfit(log10(h_arr[:4]),log10(error_arr[:4]),1)
    print("Slope: "+str(slope))
    print("Saved Relative_Error_Plot.png")
    show()
