from numpy import *
import fys3150_project2_jacobi_rot_max

def Jacobi_Method(A):
    n = len(A[0])
    tol = 1E-10 #Limit which gives off-diagonal elements zero
    iter = 0
    iter_max = 10E+8 #Number of iterations
    offdiag_max = 1000.0 #Just a number bigger that tol

    A_new = A.copy() #Copying A, as so not to overwrite it
    R = eye(n) #Identity matrix for assigning eigenvalues

    while (fabs(offdiag_max) > tol and iter <= iter_max):
        #l = 0
        #k = 0
        l,k = fys3150_project2_jacobi_rot_max.max_offdiag(A_new,n)
        offdiag_max = A_new[l][k] #Updating max off-diagonal element of A
        fys3150_project2_jacobi_rot_max.JacobiRotate(A_new,R,l,k,n)
        iter += 1
    lamdas = diagonal(A_new)

    return lamdas, [R[:,i] for i in range(0,n)], iter
