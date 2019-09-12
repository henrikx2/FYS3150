from numpy import *


def similarityTransform():
    return 0

def offdiag(A,p,q,n):
    max = 0
    for i in range(0,n):
        for j in range(i+1,n):
            aij = fabs(A[i][j])
            if aij > max:
                max = aij
                p = i
                q = j

def JacobiRotate(A,R,k,l,n):
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


tol = 1E-10
iterations = 0
while (maxnondiag > tol and iterations <= maxiter):
    #p =
    #q =
    offdiag(A,p,q,n)
    JacobiRotate(A,R,p,q,n)
    iterations += 1
