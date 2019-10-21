import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join

onlyFiles = [f for f in listdir("../Data/") if isfile(join("../Data/", f))]

gaussIter = 11
mcIter = 7

legN = []
legInt = []
legError = []
legTime = []

lagN = []
lagInt = []
lagError = []
lagTime = []

mcbfN = []
mcbfInt = []
mcbfError = []
mcbfSTD = []
mcbfVar= []
mcbfTime = []

mcisN = []
mcisInt = []
mcisError = []
mcisSTD = []
mcisVar= []
mcisTime = []

for file in onlyFiles[0:gaussIter]:
    lagFile = open("../Data/"+str(file), "r")
    lagN.append(float(lagFile.readline().split()[-1]))
    lagInt.append(float(lagFile.readline().split()[-1]))
    lagError.append(float(lagFile.readline().split()[-1]))
    lagTime.append(float(lagFile.readline().split()[-2]))
    lagFile.close()

for file in onlyFiles[gaussIter:2*gaussIter]:
    legFile = open("../Data/"+str(file), "r")
    legN.append(float(legFile.readline().split()[-1]))
    legInt.append(float(legFile.readline().split()[-1]))
    legError.append(float(legFile.readline().split()[-1]))
    legTime.append(float(legFile.readline().split()[-2]))
    legFile.close()

for file in onlyFiles[2*gaussIter:2*gaussIter+mcIter]:
    mcbfFile = open("../Data/"+str(file), "r")
    mcbfN.append(float(mcbfFile.readline().split()[-1]))
    mcbfInt.append(float(mcbfFile.readline().split()[-1]))
    mcbfError.append(float(mcbfFile.readline().split()[-1]))
    mcbfSTD.append(float(mcbfFile.readline().split()[-1]))
    mcbfVar.append(float(mcbfFile.readline().split()[-1]))
    mcbfTime.append(float(mcbfFile.readline().split()[-2]))
    mcbfFile.close()

for file in onlyFiles[2*gaussIter+mcIter:2*gaussIter+2*mcIter]:
    mcisFile = open("../Data/"+str(file), "r")
    mcisN.append(float(mcisFile.readline().split()[-1]))
    mcisInt.append(float(mcisFile.readline().split()[-1]))
    mcisError.append(float(mcisFile.readline().split()[-1]))
    mcisSTD.append(float(mcisFile.readline().split()[-1]))
    mcisVar.append(float(mcisFile.readline().split()[-1]))
    mcisTime.append(float(mcisFile.readline().split()[-2]))
    mcisFile.close()

exact = (5*np.pi**2)/(16**2)

plt.plot(np.linspace(0,450,25),np.array([0 for i in range(25)]),label="ZERO")
plt.plot(legTime,legError,label="guLeg")
plt.plot(lagTime,lagError,label="guLag")
plt.plot(mcbfTime,mcbfVar,label="MCBF")
plt.plot(mcisTime,mcisVar,label="MCIS")
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Error")
plt.savefig("../Plots/errorConvergance.pdf")
#plt.show()
plt.clf()
plt.cla()
plt.loglog(mcbfN,mcbfVar,label="MCBF")
plt.loglog(mcisN,mcisVar,label="MCIS")
plt.legend()
plt.xlabel("$log_{10}(N)$")
plt.ylabel("$log_{10}(\\sigma^{2})$")
plt.savefig("../Plots/varianceOfN.pdf")
plt.show()
