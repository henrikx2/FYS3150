from numpy import *
from numba import jit

#@jit

def Jacobi_Method(lamdas,R,A):
    n = len(A[0])
    tol = 1E-10 #Limit which gives off-diagonal elements zero
    iter = 0
    iter_max = 1E+6 #Number of iterations
    offdiag_max = 1000.0 #Just a number bigger that tol

    A_new = A.copy() #Copying A, as so not to overwrite it
    R = eye(n) #Identity matrix for eigenvalues

    while (fabs(offdiag_max) > tol and iter <= iter_max):
        l = 0.0
        k = 0.0
        max_offdiag(A_new,l,k,n)
        JacobiRotate(A_new,R,l,k,n)
        offdiag_max = A_new[l][k] #Writing new max element to A
        iter += 1

    lamdas = array([i for i in A_new[i][i]])

def max_offdiag(A,l,k,n): #Finding max off-diagonal element and choosing rotation indexes
    max = 0.0
    for i in range(0,n):
        for j in range(i+1,n):
            aij = fabs(A[i][j])
            if aij > max: #If new element is bigger than last, overwrite
                max = aij
                l = i
                k = j
    return l,k

def JacobiRotate(A,R,l,k,n):
    s = 0
    c = 0
    if A[k][l] != 0.0:
        t = 0
        tau = (A[l][l]-A[k][k])/(2*A[k][l])
        if tau >= 0:
            t = 1.0/(-tau+sqrt(1.0+tau*tau))
        else:
            t = -1.0/(-tau+sqrt(1.0+tau*tau))
        c = 1/sqrt(1+t**2)
        s = c*t
    else:
        c = 1.0
        s = 0.0

    a_kk = A[k][k]
    a_ll = A[l][l]
    A[k][k] = c**2*a_kk - 2.0*c*s*A[k][l] + s**2*a_ll
    A[l][l] = s**2*a_kk + 2.0*c*s*A[k][l] + c**2*a_ll
    A[k][l] = 0.0
    A[l][k] = 0.0
    for i in range(0,n):
        if (i != k and i != l):
            a_ik = A[i][k]
            a_il = A[i][l]
            A[i][k] = c*a_ik - s*a_il
            A[k][i] = A[i][k]
            A[i][l] = c*a_il + s*a_ik
            A[l][i] = A[i][l]
        r_ik = R[i][k]
        r_il = R[i][l]
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    return 0
