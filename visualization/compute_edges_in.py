import scipy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math as m
import sys
import random as r
name = sys.argv[1]
edges = sys.argv[2]
draw_ch = (len(sys.argv) > 3)

file1 = open(name, "r")
file2 = open(edges, "r")

X = []
S = []

partition = []
sets=0
for line in file1.readlines():
    if line == "sets\n":
        sets = 1
    else:
        temp = []
        for p in line.split(","):
            if p != "\n":
                temp.append(float(p))
        if sets == 0:
                X.append(temp)
        elif sets == 1:
                S.append(temp)

part = [0]*len(X)
for i in range(len(X)):
    for j in range(len(S)):
        if S[j][i]:
            part[i] = j

c = 0
total = 0

for line in file2.readlines():
    total += 1
    temp = line.split(" ")
    if part[int(temp[0])-1] != part[int(temp[1])-1]:
        c += 1

print(c, total, c/total*100)

file1.close()
file2.close()
