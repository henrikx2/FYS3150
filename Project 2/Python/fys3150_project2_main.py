from numpy import *
import fys3150_project2_jacobi_method

#Define equation variables
n = 50
h = 1/n

#Build matrix A
A = zeros((n,n))
A[0][0] = 2; A[0][1] = -1
for i in range(1,n):
    A[i][i] = 2
    A[i][i-1] = -1
    A[i-1][i] = -1
A[n-1][n-1] = 2; A[n-1][n-2]
#A *= 1/h**2
eigenvects = zeros((n,n))
lamdas = zeros(n)

print("Numpy eigvals:"+str(linalg.eigvals(A)))

lamdas = fys3150_project2_jacobi_method.Jacobi_Method(A,eigenvects,lamdas)
print("Jacobi eigvals:"+str(lamdas))
