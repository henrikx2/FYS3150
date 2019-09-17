from numpy import *
from matplotlib.pyplot import *

#make path to store figures and fecht data
main_path = os.getcwd()
data_path = os.path.normpath(path_main+os.sep+os.pardir+"\Data")
figure_path = os.path.normpath(path_main+os.sep+os.pardir+"\Figures")
try:
    os.mkdir(figure_path)
except OSError:
    print ("Creation of the directory %s failed. It may already exist." % figure_path)
else:
    print ("Successfully created the directory %s " % figure_path)

#Three lowest states in single electron system

n0,psi0 = readFileTwoColoumns("data_path\E_0.dat")
n1,psi1 = readFileTwoColoumns("data_path\E_1.dat")
n2,psi2 = readFileTwoColoumns("data_path\E_2.dat")

n0

#Lowest state in two-electron system without potential and 4 different frequencies

#Lowest state in two-electron system with potential and 4 different frequencies


def readFileTwoColoumns(inpFile):
    infile = open(inpFile,"r")
    first = []
    second = []
    for line in infile:
        first.append(line.split()[0])
        second.append(line.split()[1])
    return first,second

def readFileThreeColoumns(inpFile):
    infile = open(inpFile,"r")
    first = []
    second = []
    third = []
    for line in infile:
        first.append(line.split()[0])
        second.append(line.split()[1])
        third.append(line.split()[2])
    return first,second,third
