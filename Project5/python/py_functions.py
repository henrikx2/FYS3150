from numpy import array, zeros

def readarrays(filename):
	values = []
	with open(filename) as file:
		arr = []
		for line in file:
			if line != "\n":
				arr.append(float(line))
			if line == "\n":
				values.append(arr)
				arr = []
	return values

def readarrays_sideways(filename):
	values = open(filename, "r")
	lines = values.readlines()

	#Counting
	C = 0
	D = 0
	Dims = []
	A = []

	for i in lines:
		if i != "\n":
			D += 1
		if i == "\n":
			C += 1
			Dims.append(D)
			A.append(zeros(D))
			D = 0

	#Filling
	F = 0
	G = 0
	for i in lines:
		if i != "\n":
			A[F][G] = i
			G += 1
		if i == "\n":
			F += 1
			G = 0
	values.close()
	return A,len(A)

def readmatrices(filename):
	values = open(filename, "r")
	#print values.read()
	lines = values.readlines()

	#Counting
	C = 0
	D = 0
	Dims = []
	A = []

	for i in range(len(lines)):
		if lines[i] != "\n":
			D += 1
		if lines[i] == "\n":
			C += 1
			Dims.append(D)
			A.append(zeros(shape=(D,len(lines[i-1].split()))))
			D = 0

	#Filling
	F = 0
	G = 0
	for i in lines:
		if i != "\n":
			for j in range(len(i.split())):
				A[F][G][j] = i.split()[j]
			G += 1
		if i == "\n":
			F += 1
			G = 0
	values.close()

	return A,len(A)
