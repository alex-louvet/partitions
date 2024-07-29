import scipy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math as m
import sys
import random as r
name = sys.argv[1]

file1 = open(name, "r")

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

for s in S:
    temp = []
    for i in range(len(X)):
        if s[i]:
            temp.append(X[i])
    partition.append(temp)

n = sum([len(part) for part in partition])
t = len(partition)

legend = name[:-4]
        
if len(partition[0][0]) != 2:
    raise ValueError('Wrong dimension to plot: required 2, got ', len(partition[0][0]))

colors = cm.viridis(np.linspace(0, 1, len(partition)))
for p,c in zip(partition,colors):
    #ch = scipy.spatial.ConvexHull(p)
    plt.scatter([x[0] for x in p[:1]], [x[1] for x in p[:1]], marker='.',s=3, label='partition' + str(i+1),zorder=100,color='blue')
    plt.scatter([x[0] for x in p[:1]], [x[1] for x in p[:1]], marker='.',s=3, label='partition' + str(i+1),zorder=100,color='blue')
    #for simplex in ch.simplices:
        #plt.plot([p[simplex[0]][0],p[simplex[1]][0]], [p[simplex[0]][1], p[simplex[1]][1]], 'k-',zorder=101)
plt.axis('equal')
plt.title('Sample size $\\frac{n}{4}$', fontsize=25)
ax = plt.gca()
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
plt.axis('off')
plt.savefig(name[:-4]+'.png',dpi=300)
plt.clf()

file1.close()

sam = r.sample(X,t)
plt.scatter([x[0] for x in sam], [x[1] for x in sam], marker='.',s=3, color='blue')
plt.axis('equal')
ax = plt.gca()
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
plt.axis('off')
plt.savefig(name[:-4]+'approx.png',dpi=300)