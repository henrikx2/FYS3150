from numpy import *
from matplotlib.pyplot import *
import os

def readFileTwoColoumns(inpFile):
    infile = open(inpFile,"r")
    first = []
    second = []
    for line in infile:
        split = line.split()
        first.append(float(split[0]))
        second.append(float(split[1]))
    return first,second

def readFileThreeColoumns(inpFile):
    infile = open(inpFile,"r")
    first = []
    second = []
    third = []
    for line in infile:
        split = line.split()
        first.append(float(split[0]))
        second.append(float(split[1]))
        third.append(float(split[2]))
    return first,second,third

def plot_data():

    #make path to store figures and fecht data
    main_path = os.getcwd()
    data_path = os.path.normpath(main_path+os.sep+os.pardir+"\Data")
    figure_path = os.path.normpath(main_path+os.sep+os.pardir+"\Figures")
    try:
        os.mkdir(figure_path)
    except OSError:
        print ("Creation of the directory %s failed. It may already exist." % figure_path)
    else:
        print ("Successfully created the directory %s " % figure_path)

    #Three lowest states in single electron system
    rho0_single0,psi0_single = readFileTwoColoumns(data_path+"\E_0.dat")
    rho1_single1,psi1_single = readFileTwoColoumns(data_path+"\E_1.dat")
    rho2_single2,psi2_single = readFileTwoColoumns(data_path+"\E_2.dat")

    plot(rho0_single0,psi0_single)
    title("GS, single electron")
    xlabel(r"$\rho$")
    ylabel("$|\psi(\\rho)|^{2}$")
    savefig(figure_path+"\E_0.png")
    print("Saved: E_0.png")
    clf()

    plot(rho1_single1,psi1_single)
    title("1st exited state, single electron")
    xlabel(r"$\rho$")
    ylabel("$|\psi(\\rho)|^{2}$")
    savefig(figure_path+"\E_1.png")
    print("Saved: E_1.png")
    clf()

    plot(rho2_single2,psi2_single)
    title("2nd exited state, single electron")
    xlabel(r"$\rho$")
    ylabel("$|\psi(\\rho)|^{2}$")
    savefig(figure_path+"\E_2.png")
    print("Saved: E_2.png")
    clf()

    #Lowest state in two-electron system without potential and 4 different frequencies
    rho0_double_001, psi_001 = readFileTwoColoumns(data_path+"\E_GS_0_0_1_No_Interaction.dat")
    rho0_double_05, psi_05 = readFileTwoColoumns(data_path+"\E_GS_0_5_No_Interaction.dat")
    rho0_double_1, psi_1 = readFileTwoColoumns(data_path+"\E_GS_1_No_Interaction.dat")
    rho0_double_5, psi_5 = readFileTwoColoumns(data_path+"\E_GS_5_No_Interaction.dat")

    plot(rho0_double_001,psi_001, label = "$\\omega_{r} = 0.01$")
    plot(rho0_double_05,psi_05, label = "$\\omega_{r} = 0.5$")
    plot(rho0_double_1,psi_1, label = "$\\omega_{r} = 1$")
    plot(rho0_double_5,psi_5, label = "$\\omega_{r} = 5$")
    title("GS, two electrons, varying HO-frequency, no interaction.")
    xlabel(r"$\rho$")
    ylabel("$|\psi(\\rho)|^{2}$")
    legend()
    savefig(figure_path+"\E_GS_No_Interaction.png")
    print("Saved: E_GS_No_Interaction.png")
    clf()


    #Lowest state in two-electron system with potential and 4 different frequencies
    rho0_double_001_int, psi_001_int = readFileTwoColoumns(data_path+"\E_GS_0_0_1_With_Interaction.dat")
    rho0_double_05_int, psi_05_int = readFileTwoColoumns(data_path+"\E_GS_0_5_With_Interaction.dat")
    rho0_double_1_int, psi_1_int = readFileTwoColoumns(data_path+"\E_GS_1_With_Interaction.dat")
    rho0_double_5_int, psi_5_int = readFileTwoColoumns(data_path+"\E_GS_5_With_Interaction.dat")

    plot(rho0_double_001_int,psi_001_int, label = "$\\omega_{r} = 0.01$")
    plot(rho0_double_05_int,psi_05_int, label = "$\\omega_{r} = 0.5$")
    plot(rho0_double_1_int,psi_1_int, label = "$\\omega_{r} = 1$")
    plot(rho0_double_5_int,psi_5_int, label = "$\\omega_{r} = 5$")
    title("GS, two electrons, varying HO-frequency, with interaction.")
    xlabel(r"$\rho$")
    ylabel("$|\psi(\\rho)|^{2}$")
    legend()
    savefig(figure_path+"\E_GS_With_Interaction.png")
    print("Saved: E_GS_With_Interaction.png")
    clf()

    #Number of iterations as function of steps and value of eigenvalue approximation as function of steps
    n_step, n_iter, lamda = readFileThreeColoumns(data_path+"\iterations_transformations_eigenvalues_rho_max_3.dat")

    plot(n_step,n_iter, label="$\\rho = 3.0$")
    title("Number of transformations as function of integration points n.")
    xlabel("n")
    ylabel("Number of orthogonal transformations")
    legend()
    savefig(figure_path+"\\Number_Of_Transformations.png")
    print("Saved: Number_Of_Transformations.png")
    clf()

    plot(n_step,lamda, label="$\\rho = 3.0$")
    title("Eigenvalue approximation of GS as function of integration points n.")
    xlabel("n")
    ylabel("Eigenvalue approximation")
    legend()
    savefig(figure_path+"\Eigenvalues_Approx.png")
    print("Saved: Eigenvalues_Approx.png")
    clf()
