import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np
import math as m
import matplotlib as mpl
name = sys.argv[1]

dlist = [2,3,4,5,10]
tlist = [16,32,64,128,256]

var = "d"

label_table = ["","minWeight","greedy","sampling","greedy","atOnce","atOnce","atOnce","atOnce","atOnce","atOnce"]
color_table = ["", "red","green","orange","green","purple","purple","purple","purple","purple","purple"]
'''
for x in dlist:

    df = pd.read_csv("test.csv", sep=";")
    toplot = df.drop(df[(df["n"] != int(sys.argv[1])) | (df["d"] != x) |(df["t"] != int(sys.argv[3])) | (df["ss_type"] != sys.argv[4])  | (df["algo"] >= 12)].index)

    cross = toplot.groupby('algo').mean("max_crossing")
    print(cross)
    res_1.append(cross.iat[0,10])
    res_10.append(cross.iat[1,10])
    res_approx.append(1/2 * (cross.iat[1,11] + cross.iat[0,11]))

plt.plot(dlist, res_1, c='red', label="min algo")
plt.plot(dlist, res_10, c='green', label="parallel algo")
plt.plot(dlist, res_approx, c='blue', label="uniform sample")
plt.plot(dlist, [1/m.sqrt(int(sys.argv[3])) for d in dlist], c='grey')
plt.plot(dlist, [1/int(sys.argv[3])**((d+1)/(d)) for d in dlist], c='grey')
'''

fig, ax = plt.subplots()

df = pd.read_csv(name, sep=";")
toplot = df.drop(df[(df["n"] != int(sys.argv[1])) | (df["t"] != int(sys.argv[3])) | (df["ss_type"] != sys.argv[4])  | (df["algo"] >= 12)].index)
for key, grp in toplot.groupby('algo'):
    temp = grp.groupby([var],as_index=False).agg(min=pd.NamedAgg(column="approx_part", aggfunc="min"), max=pd.NamedAgg(column="approx_part", aggfunc="max"))
    print(temp)
    grp.groupby([var]).mean("approx_part").plot(ax=ax,y="approx_part", label=label_table[int(key)], color=color_table[int(key)])

    ax.fill_between(x=var,y1='min',y2='max',data=temp,color=mpl.colors.to_rgba(color_table[int(key)], 0.05))

temp = toplot.groupby([var],as_index=False).agg(min=pd.NamedAgg(column="approx_random", aggfunc="min"), max=pd.NamedAgg(column="approx_random", aggfunc="max"))
print(temp)
grp.groupby([var]).mean("approx_random").plot(ax=ax,y="approx_random", label="random", color='blue')

ax.fill_between(x=var,y1='min',y2='max',data=temp,color=mpl.colors.to_rgba('blue', 0.05))

#plt.plot(tlist, [1/m.sqrt(t) for t in tlist], c='grey', linestyle='dotted', label="$\\frac{1}{\\sqrt{t}}$")
#plt.plot(tlist, [1/m.sqrt(t)**((int(sys.argv[2])+1)/(int(sys.argv[2]))) for t in tlist], c='grey', linestyle='dashdot', label="$\\frac{1}{\\sqrt{t^{\\frac{d+1}{d}}}}$")
plt.plot(dlist, [1/m.sqrt(int(sys.argv[3])) for d in dlist], c='grey', linestyle='dotted', label="$\\frac{1}{\\sqrt{t}}$")
plt.plot(dlist, [1/m.sqrt(int(sys.argv[3]))**((d+1)/d) for d in dlist], c='grey', linestyle='dashdot', label="$\\frac{1}{\\sqrt{t^{\\frac{d+1}{d}}}}$")

ax.set_xlabel(var, fontsize=25)
ax.set_ylabel("epsilon",fontsize=25)

plt.legend(fontsize=20)
plt.show()