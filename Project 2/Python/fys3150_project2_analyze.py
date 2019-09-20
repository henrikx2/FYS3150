from numpy import *
import fys3150_project2_hamiltonian_matrix
import fys3150_project2_jacobi_method

rho_min = 10E-10

#Calculates number of transformations and eigenvalues in ground state as function of iterations n
def numberOfTransformationsAndEigenvaluesOfN(path):
    n_max = 200
    n = linspace(10,n_max,n_max/10)
    rho_max = 3

    file = open(path+"\iterations_transformations_eigenvalues_rho_max_%i.dat"%rho_max, "w")
    for i in n:
        i = int(i)
        H = fys3150_project2_hamiltonian_matrix.hamiltonian(rho_min,rho_max,i)
        lamdas, eigenvects, iter = fys3150_project2_jacobi_method.Jacobi_Method(H)
        lamda_GS = min(lamdas)
        file.write("%i   %i    %f\n"%(i,iter,lamda_GS))
    file.close()
    print("Saved: iterations_transformations_eigenvalues_rho_max_%i.dat"%rho_max)

def threeLowestStates4Decimal(path):
    print("One electron, three lowest states:")
    rho_max = [3.0,3.7,4.2] #Experimentally obtained values by manuelly optimizing the eigenvalues
    n = [325,336,254]  #Experimentally obtained values by manuelly optimizing the eigenvalues

    state = 0
    for i,k in zip(n,rho_max):
        rho_int = linspace(rho_min,k,i)
        H = fys3150_project2_hamiltonian_matrix.hamiltonian(rho_min,rho_max[state],i)
        lamdas, eigenvects, iter = fys3150_project2_jacobi_method.Jacobi_Method(H)
        index = argsort(lamdas)

        print("Eigenvalue %i-state: %f"%(state,lamdas[index[state]]))

        lamda_state = lamdas[index[state]]
        file = open(path+"\E_%i.dat"%state, "w")
        for dex in range(0,i):
            file.write("%f    %f\n"%(rho_int[dex],eigenvects[index[state]][dex]**2))
        file.close()
        print("Saved: E_%i.dat"%state)
        state += 1

def lowestStateTwoElectronsNoInteraction(path):
    print("Two electrons, no interaction:")
    potential = False
    rho_max = 30
    n = 300
    rho_int = linspace(rho_min,rho_max,n)

    omega_r = array([0.01,0.5,1,5])
    omega_name = array(["0_0_1","0_5","1","5"])

    for i in range(0,len(omega_r)):
            H = fys3150_project2_hamiltonian_matrix.hamiltonian_two_electrons(rho_min,rho_max,omega_r[i],n,potential)
            lamdas, eigenvects, iter = fys3150_project2_jacobi_method.Jacobi_Method(H)
            index = argsort(lamdas)

            print("Omega = %.2f   Eigenvalue GS: %f"%(omega_r[i],lamdas[index[0]]))

            file = open(path+"\E_GS_%s_No_Interaction.dat"%omega_name[i], "w")
            for k in range(0,n):
                file.write("%f    %f\n"%(rho_int[k],eigenvects[index[0]][k]**2))
            file.close()
            print("Saved: E_GS_%s_No_Interaction.dat"%omega_name[i])

def lowestStateTwoElectronsWithInteraction(path):
    print("Two electrons, with interaction:")
    potential = True
    rho_max = 40
    n = 300
    rho_int = linspace(rho_min,rho_max,n)
    omega_r = array([0.01,0.5,1,5])
    omega_name = array(["0_0_1","0_5","1","5"])

    for i in range(0,len(omega_r)):
            H = fys3150_project2_hamiltonian_matrix.hamiltonian_two_electrons(rho_min,rho_max,omega_r[i],n,potential)
            lamdas, eigenvects, iter = fys3150_project2_jacobi_method.Jacobi_Method(H)
            index = argsort(lamdas)

            print("Omega = %.2f   Eigenvalue GS: %f"%(omega_r[i],lamdas[index[0]]))

            file = open(path+"\E_GS_%s_With_Interaction.dat"%omega_name[i], "w")
            for k in range(0,n):
                file.write("%f    %f\n"%(rho_int[k],eigenvects[index[0]][k]**2))
            file.close()
            print("Saved: E_GS_%s_With_Interaction.dat"%omega_name[i])
