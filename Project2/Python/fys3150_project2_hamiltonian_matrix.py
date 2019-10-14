from numpy import *
from numba import jit

@jit

def hamiltonian(rho_min,rho_max,n):
    h = (rho_max-rho_min)/n
    rho_int = linspace(rho_min,rho_max,n)
    H = zeros((n,n))

    V = rho_int**2

    e = -1/h**2 #Constant on lower/upper diagonal
    d1 = 2/h**2
    H[0][0] = d1 + V[0]; #Initialize first diagonal element in Hamiltonian
    for i in range(1,n):
        H[i][i] = d1 + V[i] #Initialize main diagonal in Hamiltonian
        H[i][i-1] = e #Initialize lower diagonal
        H[i-1][i] = e #Initialize upper diagonal
    return H

@jit

def hamiltonian_two_electrons(rho_min,rho_max,omega_r,n,potential):
    h = (rho_max-rho_min)/n
    rho_int = linspace(rho_min,rho_max,n)
    H = zeros((n,n))
    V = zeros(n)

    if potential == True:
        V = omega_r**2*rho_int**2 + 1/rho_int
    else:
        V = omega_r**2*rho_int**2

    e = -1/h**2 #Constant on lower/upper diagonal
    d1 = 2/h**2
    H[0][0] = d1 + V[0]; #Initialize first diagonal element in Hamiltonian
    for i in range(1,n):
        H[i][i] = d1 + V[i] #Initialize main diagonal in Hamiltonian
        H[i][i-1] = e #Initialize lower diagonal
        H[i-1][i] = e #Initialize upper diagonal
    return H
