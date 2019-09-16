from numpy import *
from matplotlib.pyplot import*
from tabulate import tabulate
import fys3150_project2_hamiltonian_matrix
import fys3150_project2_jacobi_method
import os
import sys

path_main = os.getcwd()
path = os.path.normpath(path_main+os.sep+os.pardir+"\Data")
try:
    os.mkdir(path)
except OSError:
    print ("Creation of the directory %s failed. It may already exist." % path)
else:
    print ("Successfully created the directory %s " % path)

#Define equation variables
rho_min = 10E-10
rho_max = 3
n = 325 #Number of integration points/size of hamiltonian

#Build tridiagonal Hamiltonian for 1 electron in HO-potential
H = fys3150_project2_hamiltonian_matrix.hamiltonian(rho_min,rho_max,n)

omega_r = 0.01
#H = fys3150_project2_hamiltonian_matrix.hamiltonian_two_electrons(rho_min,rho_max,omega_r,n,potential=True)

eigenvects = zeros((n,n)) #Matrix to hold eigenvectors
lamdas = zeros(n) #Array to hold eigenvalues

lamdas_numpy = sort(linalg.eigvals(H))[::-1]

lamdas, eigenvects = fys3150_project2_jacobi_method.Jacobi_Method(H,eigenvects,lamdas)

#eigen_set = zip(lamdas,[i[k] for i in eigenvects)

index = lamdas.argmin()

plot(linspace(rho_min,rho_max,n),[i[index] for i in eigenvects])
show()



'''
#Print eigenvalues in table
eigvals_zip = [[jac,num] for jac,num in zip(lamdas,lamdas_numpy)]
eigvals_table = tabulate(eigvals_zip, headers=["Jacobi eigvals:","Numpy eigvals:"])
print(eigvals_table)

#Save eigenvals to table
eigvals_file = open(path+"\Jacobi_Eigvals_vs_Numpy_Eigvals.txt", "w")
eigvals_file.write("Jacobi eigvals:          Numpy eigvals:\n")
for i in eigvals_zip:
    eigvals_file.write("%s%25s\n" %(round(i[0],5),round(i[1],5)))
eigvals_file.close()
'''
