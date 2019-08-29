from numpy import *
import os
import sys

n = [10,100,1000]

for i in n:
    aFile = open(os.path.join(sys.path[0], "a"+str(i)),"w")
    bFile = open(os.path.join(sys.path[0], "b"+str(i)),"w")
    cFile = open(os.path.join(sys.path[0], "c"+str(i)),"w")

    a = full(i,-1,int)
    b = full(i,2,int)
    c = a

    for i in range(0,len(a)):
        aFile.write(str(a[i])+",")
        bFile.write(str(a[i])+",")
        cFile.write(str(a[i])+",")

    aFile.close()
    bFile.close()
    cFile.close()
