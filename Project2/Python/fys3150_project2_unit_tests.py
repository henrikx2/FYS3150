from numpy import *
import fys3150_project2_jacobi_rot_max

def test_eigenvalues(A):
    analytical = [3,7,11]
    calculated = sort(linalg.eigvals(A))

def test_max_offdiag():
    A = array([
    [1,0,2,5],
    [0,1,2,0],
    [2,2,1,2],
    [5,0,2,1]])
    l,k = fys3150_project2_jacobi_rot_max.max_offdiag(A,4)
    if A[l][k] != 5:
        print("max_offdiag() algorithm not working properly!")

def test_JacobiRotate(tol):
    A = array([
    [1,0,2,5],
    [0,1,2,0],
    [2,2,1,2],
    [5,0,2,1]])
    n = len(A[0])
    l = 0
    k = 3
    R = eye(n)
    fys3150_project2_jacobi_rot_max.JacobiRotate(A,R,l,k,n)
    R_T = transpose(R)
    dotproduct = dot(R,R_T)
    for i in range(0,n):
        for m in range(0,n):
            if i != k and dotproduct[i][k] > tol:
                print("JacobiRotate() does not preserve orthogonality!")
                break
