from numpy import *
from matplotlib.pyplot import*
import fys3150_project2_hamiltonian_matrix
import fys3150_project2_jacobi_method
import fys3150_project2_analyze
import fys3150_project2_plot_data
import os
import sys
import time

analyze = input("Do full analysis? (approx. 2 min)\n(y/n): ")

plot_data = 0
if analyze != "y":
    plot_data = input("Just plot precalculated data?\n(y/n): ")

if analyze != "y" and plot_data != "y":
    #Define equation variables
    rho_min = 10E-10
    rho_max = float(input("Choose rho_max: "))
    n = int(input("Choose number of integration points (int): ")) #Number of integration points/size of hamiltonian
    rho = linspace(rho_min,rho_max,n)
    H = zeros((n,n))

    single_double = float(input("One electron or two electrons in HO?\n(1/2): "))

    if single_double == 1:
        #Build tridiagonal Hamiltonian for 1 electron in HO-potential
        H = fys3150_project2_hamiltonian_matrix.hamiltonian(rho_min,rho_max,n)

    if single_double == 2:
        pot = input("With coulomb interactions?\n(y/n): ")
        omega_r = float(input("Choose frequency of HO: "))
        potential = False
        if pot == "y":
            potential = True
        H = fys3150_project2_hamiltonian_matrix.hamiltonian_two_electrons(rho_min,rho_max,omega_r,n,potential)

    lamdas, eigenvects, iter = fys3150_project2_jacobi_method.Jacobi_Method(H)

    print("The lowest eigenvalues of Hamiltonian:")
    sorted = sort(lamdas)
    print(sorted[0])
    print(sorted[1])
    print(sorted[2])
    print("Number of orthogonal transformations:")
    print(iter)


'''Full analysis and save data to files'''
if analyze == "y":
    start_time = float(time.time())
    #Make direcory for data files
    path_main = os.getcwd()
    path = os.path.normpath(path_main+os.sep+os.pardir+"\Data")
    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed. It may already exist." % path)
    else:
        print ("Successfully created the directory %s " % path)

    #Number of orthogonal transformations and eigenvalue as function of n
    fys3150_project2_analyze.numberOfTransformationsAndEigenvaluesOfN(path)

    #Three lowest state single electron 4 decimal precision
    fys3150_project2_analyze.threeLowestStates4Decimal(path)

    #Lowest state two-electron system no interaction
    fys3150_project2_analyze.lowestStateTwoElectronsNoInteraction(path)

    #Lowest state two-electron system with interaction
    fys3150_project2_analyze.lowestStateTwoElectronsWithInteraction(path)

    elapsed_time = float(time.time()-start_time)
    print("All data saved! Time of full analysis: %i seconds" %elapsed_time)

if plot_data == "y":
    fys3150_project2_plot_data.plot_data()
