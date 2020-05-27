#! /usr/bin/env python3

from pyfinite import ffield
import numpy as np
from numpy import random

q = 5
F = ffield.FField(q) # finite field GF(2^5)
N = 2**q





#alpha^j_l = i^{k[j,l]}
K = np.zeros((N,N), dtype=np.uint8)

for l in range(N):
    print(l)
    for m in range(N):
        for n in range(N):
            v = F.Multiply(((l >> m)%2) << m, ((l >> n)%2) << n)
            for j in range(N):
                K[j,l] += F.Multiply(j, v)

for l in range(N):
    for j in range(N):
        K[j,l] = K[j,l] % 4    

basisVec = np.zeros((N+1,N,N), dtype=np.complex)        
for j in range(N):
    basisVec[N, j, j] = 1

for j in range(N):
    for k in range(N):
        basisVec[0,j,k] = (-1)**F.Multiply(F.Subtract(0, k), j)/ np.sqrt(N)


for i in range(1,N):
    for k in range(N):
        for l in range(N):
            basisVec[i,k,l] = (-1)**F.Multiply(k,l) *  (1j)**(-K[i,l]) / np.sqrt(N)
        
def delta(a,b):
    if a == b:
        return 1
    return 0
        
for _ in range(1000):
    i, j, = np.random.randint(0,N+1, 2)
    k, l, = np.random.randint(0,N, 2)

    d = 1/N + delta(i,j)*(delta(k,l) - 1/N)
    a = abs(np.vdot(basisVec[i,k], basisVec[j,l]))**2
    if abs(a - d)  > 1e-10:
        print("error", a, d)
    

                     



