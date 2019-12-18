from numpy import *
from matplotlib.pyplot import *
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D

def plot_minima(variation,energies,variances,trail,parameter):

    fig1, axs = subplots(2,1,sharex=True,gridspec_kw={'hspace': 0})
    fig1.text(0.5, 0.94, "Mean energy and energy variance of $\\Psi_{T"+str(trail)+"}$ as function of "+parameter, ha='center', va='center')
    axs[0].plot(variation,energies,label="Mean energy")
    axs[0].legend()
    axs[0].set_ylabel("Mean energy, $E$")
    axs[1].plot(variation,variances,label="Energy variance")
    axs[1].legend()
    axs[1].set_ylabel("Variance, $\sigma^{2}$")
    axs[1].set_xlabel("Variational parameter, "+parameter)
    savefig("../figures/plot_minima_trail_"+trail+".pdf")
    show()

def plot_stability(n_vals,energies,variances,trail):

    fig1, axs = subplots(2,1,sharex=True,gridspec_kw={'hspace': 0})
    fig1.text(0.5, 0.94, str(trail)+". trail wave function: Energy and variance as function of MC-cycles", ha='center', va='center')
    axs[0].plot(log10(n_vals),energies,label="Mean energy")
    axs[0].legend()
    axs[0].set_ylabel("Mean energy, $E$")
    axs[1].plot(log10(n_vals),variances,label="Energy variance")
    axs[1].legend()
    axs[1].set_ylabel("Variance, $\sigma^{2}$")
    axs[1].set_xlabel("$log_{10}$ of # of MC-cycles, n")
    savefig("../figures/plot_stability_trail_"+str(trail)+".pdf")
    show()

def plot_virial(omegas,Ek,Ep_HO,Ep_electron,trail):

    omegas = array(omegas)
    Ek = array(Ek)
    Ep_HO = array(Ep_HO)
    Ep_electron = array(Ep_electron)

    figure(1)
    plot(omegas,Ek/(Ep_HO+Ep_electron),label="With repulsive interaction")
    plot(omegas,Ek/Ep_HO,label="Without repulsive interaction")
    legend()
    ylabel("$\\langle T\\rangle$/$\\langle V\\rangle$")
    xlabel("$\\omega$")
    title(str(trail)+". Trail wave function: $\\langle T\\rangle$/$\\langle V\\rangle$ as function of $\\omega$")
    savefig("../figures/plot_virial_trail_"+str(trail)+".pdf")
    show()



def plot_best_steplength(values):

    omega = array([i[0] for i in values[0]])
    alpha = array([i[0] for i in values[1]])
    stepsize = array(values[2])

    style.use('seaborn-darkgrid')

    fig = figure()
    ax = fig.add_subplot(111, projection='3d')
    fig.text(0.5, 0.94, "Stepsize as function of $\\omega$ and $\\alpha$" , ha='center', va='center')

    X, Y = meshgrid(omega, alpha)

    ax.plot_surface(X, Y, stepsize)

    ax.set_xlabel("$\\omega$")
    ax.set_ylabel("$\\alpha$")
    ax.set_zlabel("$\\delta$")
    tight_layout()
    ax.view_init(30, 50)
    savefig("../figures/best_steplength_3D.pdf")
    show()


