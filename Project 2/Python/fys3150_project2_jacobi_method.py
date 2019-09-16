from numpy import *
import fys3150_project2_max_offdiag
import fys3150_project2_jacobi_rotate

def Jacobi_Method(A,eigenvects,lamdas):
    n = len(A[0])
    tol = 1E-10 #Limit which gives off-diagonal elements zero
    iter = 0
    iter_max = 10E+8 #Number of iterations
    offdiag_max = 1000.0 #Just a number bigger that tol

    A_new = A.copy() #Copying A, as so not to overwrite it
    eigenvects = eye(n) #Identity matrix for assigning eigenvalues

    while (fabs(offdiag_max) > tol and iter <= iter_max):
        l = 0
        k = 0
        l,k = fys3150_project2_max_offdiag.max_offdiag(A_new,l,k,n)
        offdiag_max = A_new[l][k] #Updating max off-diagonal element of A
        fys3150_project2_jacobi_rotate.JacobiRotate(A_new,eigenvects,l,k,n)
        iter += 1
    print(iter)
    lamdas = diagonal(A_new)

    return lamdas, eigenvects
