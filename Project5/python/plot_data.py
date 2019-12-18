from numpy import *
from matplotlib.pyplot import *
from scipy.optimize import curve_fit
import py_functions
import plot_functions

# Uncomment the values to plot

# 1 Plot energy and variance as functions of alpha and beta

values_1 = py_functions.readarrays("../Data/Results_5_plot_minima_trail_1.txt")
values_2 = py_functions.readarrays("../Data/Results_5_plot_minima_trail_2.txt")
plot_functions.plot_minima(values_1[0],values_1[1],values_1[2],"1","$\\alpha$")
plot_functions.plot_minima(values_2[0],values_2[1],values_2[2],"2","$\\beta$")

# 2 Plot energy stability as function of MC-cycles

mc_max = 100000
values_stability_1 = py_functions.readarrays("../Data/Results_5_plot_stability_trail_1_MC_"+str(mc_max)+".txt")
values_stability_2 = py_functions.readarrays("../Data/Results_5_plot_stability_trail_2_MC_"+str(mc_max)+".txt")


plot_functions.plot_stability(values_stability_1[0],values_stability_1[1],values_stability_1[2],1)
plot_functions.plot_stability(values_stability_2[0],values_stability_2[1],values_stability_2[2],2)


# 3 Plot virial theorem

mc_max_virial = 100000
values_virial_1 = py_functions.readarrays("../Data/Results_5_plot_virial_trail_1_MC_"+str(mc_max_virial)+".txt")
values_virial_2 = py_functions.readarrays("../Data/Results_5_plot_virial_trail_2_MC_"+str(mc_max_virial)+".txt")


plot_functions.plot_virial(values_virial_1[0],values_virial_1[1],values_virial_1[2],values_virial_1[3],1)
plot_functions.plot_virial(values_virial_2[0],values_virial_2[1],values_virial_2[2],values_virial_2[3],2)

# 4 Calculate function for best step size

def func(x,a,b,c):
    return a/(x+b)+c

values_step = py_functions.readarrays("../data/Results_5_best_steplength_MC_100000.txt")

plot(values_step[0],values_step[2], label="Data")

coeffs = polyfit(values_step[0],values_step[2],2)

popt, pcov = curve_fit(func,values_step[0],values_step[2])

x = linspace(0,1.6,100)
y = func(x,popt[0],popt[1],popt[2])
plot(x,y,label="Optimal fit")

xlabel("$\\alpha$")
ylabel("Stepsize")
title("Optimized stepsize as function of alpha")
text(1.0, 7.0, "$f(x) = \\frac{1.11587718}{x+0.11152204}+1.00514885$", ha='center', va='center')
legend()
savefig("../figures/best_steplength.pdf")
show()
