from gurobipy import *
import numpy as np
from numpy import genfromtxt
import math
n = 5000
model = 'borda'

#Read observed matrix based on estimated permutations and true matrix
sample = genfromtxt(str(n) + model + '.csv', delimiter=',')
sample = np.array(sample)


base = genfromtxt(str(n) + model + '_truth.csv', delimiter = ',')
base = (np.array(base)/8.0)


#Set up convex optimization problem for isotonic regression
#Scale reduces the number of variables by taking the mean of contiguous scale * scale submatrices
m = Model("ir1")
scale = 20
a = int(n/scale)
compressed = np.zeros((a,a))
for i in range(a):
    for j in range(a):
        compressed[i][j] = np.mean(sample[i*scale : (i+1) * scale, j*scale : (j+1) * scale])

print('preprocessing_done')

mat = m.addVars(a * a, lb = 0.0, ub = 1.0, vtype = GRB.CONTINUOUS, name = "M")
for i in range(a):
    m.addConstrs(mat[i*a + j] <= mat[i*a + j+1] for j in range(a-1))
for j in range(a):
    m.addConstrs(mat[i*a + j] <= mat[(i + 1)*a + j] for i in range(a-1))
obj = QuadExpr()
obj = quicksum((mat[i*a+ j] - compressed[i][j])*(mat[i*a + j] - compressed[i][j]) for i in range(a) for j in range(a))


#Solve the regression
m.setObjective(obj, GRB.MINIMIZE)
print("Solving")
m.optimize()


output = []
error = 0
error = sum((mat[int(i/scale)*a + int(j/scale)].x - base[i][j])* (mat[int(i/scale)*a + int(j/scale)].x - base[i][j]) for i in range(n) for j in range(n))
print(error)



f= open(str(n) + model + ".txt","w+")
f.write(str(n))
f.write('\n')
f.write(str(error))
f.close()