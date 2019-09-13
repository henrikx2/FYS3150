from numpy import *
import fys3150_project2_jacobi_method

def hamiltonian(rho_min,rho_max,n):
    h = (rho_max-rho_min)/n
    H = zeros((n,n))

    V = zeros(n)
    V[0] = rho_min**2
    for i in range(1,n):
        rho = rho_min + i*h
        V[i] = rho**2

    e = -1/h**2 #Constant on lower/upper diagonal
    H[0][0] = 2/h**2 + V[0]; #Initialize first diagonal element in Hamiltonian
    for i in range(1,n):
        H[i][i] = 2/h**2 + V[i] #Initialize main diagonal in Hamiltonian
        H[i][i-1] = e #Initialize lower diagonal
        H[i-1][i] = e #Initialize upper diagonal
    H[n-1][n-1] = 2/h**2 + V[n-1] #Initialze last element on tridiagonal
    return H
