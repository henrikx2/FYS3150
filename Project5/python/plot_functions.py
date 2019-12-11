from numpy import *
from matplotlib.pyplot import *

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
    axs[1].set_xlabel("# of MC-cycles, n")
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