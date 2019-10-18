import numpy as np
import matplotlib.pyplot as plt

r = np.linspace(-5,5,1000)
psi = np.exp(-2*r)

plt.plot(r,psi,label="$\\psi(\lambda = 3)$ = "+str(round(np.exp(-2*3),6)))
plt.xlabel("r")
plt.ylabel("$\\psi$")
plt.legend()
plt.savefig("../Plots/singleParticle.pdf")
plt.show()