from numpy import *
import fys3150_project2_hamiltonian_matrix
import fys3150_project2_jacobi_method

#Define equation variables
rho_min = 0
rho_max = 1
n = 10

#Build tridiagonal Hamiltonian
H = fys3150_project2_hamiltonian_matrix.hamiltonian(rho_min,rho_max,n)

eigenvects = zeros((n,n))
lamdas = zeros(n)

print("Numpy eigvals:"+str(sort(linalg.eigvals(H))[::-1]))

lamdas = fys3150_project2_jacobi_method.Jacobi_Method(H,eigenvects,lamdas)
print("Jacobi eigvals:"+str(lamdas))
