CC = g++
CFLAGS = -O3 -Wall -std=c++11
LFLAGS = -fopenmp

OBJS = functions.cpp
TARGETS = gaussLeg gaussLag monteCarloBF monteCarloIS monteCarloISPar

all: $(TARGETS)

gaussLeg: gaussLeg.cpp
	$(CC) $(OBJS) gaussLeg.cpp -o gaussLeg

gaussLag: gaussLag.cpp
	$(CC) $(OBJS) gaussLag.cpp -o gaussLag

monteCarloBF: monteCarloBF.cpp
	$(CC) $(OBJS) monteCarloBF.cpp -o monteCarloBF

monteCarloIS: monteCarloIS.cpp
	$(CC) $(OBJS) monteCarloIS.cpp -o monteCarloIS

monteCarloISPar: monteCarloISPar.cpp
	$(CC) -o monteCarloISPar $(CFLAGS) $(OBJS) monteCarloISPar.cpp  $(LFLAGS)
