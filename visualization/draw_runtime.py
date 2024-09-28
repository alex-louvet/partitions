import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd
import numpy as np
import math as m

label_table = ["MP-Mat","min","potential","sampling","deque","atOnce","atOnce","atOnce","atOnce","atOnce","atOnce"]
color_table = ["blue", "red","green","orange","blue","purple","purple","purple","purple","purple","purple"]
join_table = {'d': 2, 't': 512, 'n':8192}

df = pd.read_csv("results.csv", sep=";")
toplot = df.drop(df[(df["t"] != join_table['t']) | (df["d"] != join_table['d']) | (df["ss_type"] != sys.argv[2]) | (df["ss_type"] != sys.argv[2])| (df["algo"] == 2)| (df["algo"] == 4)].index).sort_values(sys.argv[1])
if sys.argv[1] == 't':
    toplot = df.drop(df[(df["n"] != join_table['n']) | (df["d"] != join_table['d']) | (df["ss_type"] != sys.argv[2]) | (df["ss_type"] != sys.argv[2])| (df["algo"] == 2)| (df["algo"] == 4)].index).sort_values(sys.argv[1])
if sys.argv[1] == 'd':
    toplot = df.drop(df[(df["t"] != join_table['t']) | (df["n"] != join_table['n']) | (df["ss_type"] != sys.argv[2]) | (df["ss_type"] != sys.argv[2])| (df["algo"] == 2)| (df["algo"] == 4)].index).sort_values(sys.argv[1])

fig, ax = plt.subplots()

for key, grp in toplot.groupby('algo'):
    temp = grp.groupby([sys.argv[1]],as_index=False).agg(min=pd.NamedAgg(column="time", aggfunc="min"), max=pd.NamedAgg(column="time", aggfunc="max"))
    grp.groupby([sys.argv[1]]).mean("time").plot(ax=ax,y="time", label=label_table[int(key)], color=color_table[int(key)])
    ax.fill_between(x=sys.argv[1],y1='min',y2='max',data=temp,color=mpl.colors.to_rgba(color_table[int(key)], 0.15))
'''
if sys.argv[1] == 'n':
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[(x*m.sqrt(x)*512) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='$nmt$')
if sys.argv[1] == 't':
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[(8192*m.sqrt(8192)*x) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='$nmt$')
if sys.argv[1] == 'd':
    plt.plot(np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10),[(8192*m.sqrt(8192)*512) for x in np.linspace(toplot.min(axis=0)[sys.argv[1]],toplot.max(axis=0)[sys.argv[1]],10)],color='black',label='$nmt$')
'''
ax.set_xlabel(sys.argv[1], fontsize=25)
ax.set_ylabel("runtime",fontsize=25)

plt.legend(loc='best', fontsize="large")
plt.savefig(sys.argv[3])
