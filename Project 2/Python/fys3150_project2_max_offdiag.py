from numpy import *
from numba import jit

@jit

#Arguments(A:matrix of size nxn,l:index integer,k:index integer,n: size of matrix)
def max_offdiag(A,l,k,n): #Finding max off-diagonal element and choosing its indexes for rotation
    max = 0.0
    for i in range(0,n):
        for j in range(i+1,n):
            aij = fabs(A[i][j])
            if aij > max: #If new element is bigger than last, overwrite
                max = aij
                l = i
                k = j
    return l,k
